# Inventory of the 16symbols-cxy-hmd tree: one row per file, plus one per directory.
# Read-only. Writes only into _reliability_audit\data\.
$ErrorActionPreference = 'Stop'

$root  = 'E:\nas-shared\Endmembers\16symbols-cxy-hmd'
$out   = Join-Path $root '_reliability_audit\data'
$mine  = Join-Path $root '_reliability_audit'

Write-Host "root = $root"
$sw = [System.Diagnostics.Stopwatch]::StartNew()

$files = New-Object System.Collections.Generic.List[object]
$dirs  = New-Object System.Collections.Generic.List[object]

Get-ChildItem -LiteralPath $root -Force | ForEach-Object {
    if ($_.PSIsContainer) {
        # one row per directory, with the count of files directly inside it
        $dirs.Add([pscustomobject]@{
            name      = $_.Name
            mtime_utc = $_.LastWriteTimeUtc.ToString('s')
        })
    } else {
        $files.Add([pscustomobject]@{
            rel       = $_.Name
            size      = $_.Length
            mtime_utc = $_.LastWriteTimeUtc.ToString('s')
        })
    }
}

# Now recurse per endpoint, skipping our own working directory.
$endpoints = $dirs | Where-Object { $_.name -ne '_reliability_audit' }
$perEndpoint = New-Object System.Collections.Generic.List[object]

foreach ($ep in $endpoints) {
    $epPath = Join-Path $root $ep.name
    $n = 0; $bytes = 0
    $top = @()
    Get-ChildItem -LiteralPath $epPath -Force | ForEach-Object {
        if ($_.PSIsContainer) {
            $sub = Get-ChildItem -LiteralPath $_.FullName -Recurse -File -Force -ErrorAction SilentlyContinue
            $n += $sub.Count
            $bytes += ($sub | Measure-Object -Property Length -Sum).Sum
            $top += ('DIR:' + $_.Name)
        } else {
            $n += 1
            $bytes += $_.Length
            $top += $_.Name
        }
    }
    $perEndpoint.Add([pscustomobject]@{
        endpoint  = $ep.name
        n_files   = $n
        n_bytes   = [int64]$bytes
        top_level = ($top -join '|')
    })
    Write-Host ("  {0,-28} {1,7} files  {2,12:N0} B" -f $ep.name, $n, $bytes)
}

$files      | Export-Csv -NoTypeInformation -Encoding UTF8 (Join-Path $out 'root_files.csv')
$endpoints  | Export-Csv -NoTypeInformation -Encoding UTF8 (Join-Path $out 'endpoints.csv')
$perEndpoint| Export-Csv -NoTypeInformation -Encoding UTF8 (Join-Path $out 'endpoint_inventory.csv')

$sw.Stop()
Write-Host ("done in {0:N1}s; endpoints={1}" -f $sw.Elapsed.TotalSeconds, $perEndpoint.Count)
