#!/usr/bin/env python3
"""Prepare an independent W90 -> TB2J rerun for all Phase-1 magnetic states."""
from __future__ import annotations
import argparse, gzip, hashlib, json, re, shutil
from pathlib import Path
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar

STATES = {
    "fm": [3.0, 3.0],
    "afm_ti_pos_nb_neg": [3.0, -3.0],
    "afm_ti_neg_nb_pos": [-3.0, 3.0],
    "weak_fm": [0.5, 0.5],
    "weak_afm": [0.5, -0.5],
    "nm": [0.0, 0.0],
}
BASE = {
    "SYSTEM":"BCC-Nb-Ti independent W90/TB2J rerun",
    "ENCUT":520,"PREC":"Accurate","ALGO":"Normal","ISTART":0,"ICHARG":2,
    "ISPIN":2,"ISMEAR":1,"SIGMA":0.2,"EDIFF":1e-6,"NELM":200,"NELMIN":6,"NELMDL":-6,
    "IBRION":-1,"ISIF":2,"NSW":0,"ISYM":0,"LREAL":False,"LASPH":True,"LMAXMIX":4,
    "ENAUG":1360,"LORBIT":11,"LWAVE":False,"LCHARG":False,"KPAR":1,"NCORE":1,
    "NBANDS":64,"NUM_WANN":18,"LWANNIER90":True,"LWANNIER90_RUN":True,
    "LWRITE_MMN_AMN":True,"GGA":"PE",
}

def sha(p):
 h=hashlib.sha256(); h.update(p.read_bytes()); return h.hexdigest()

def titles(p): return re.findall(r"TITEL\s*=\s*([^\n]+)",p.read_text(errors="replace"))

def wrapper(root):
 states=list(STATES)
 lines='\n'.join(f'  "{x}"' for x in states)
 return f'''#!/usr/bin/env bash
#SBATCH --job-name=nbti-w90j
#SBATCH --partition=partCPU
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output={root}/logs/w90j-%A_%a.out
#SBATCH --error={root}/logs/w90j-%A_%a.err
set -euo pipefail
ROOT={root!s}
PYTHON="$ROOT/../../../../.conda/bin/python"
STATES=(
{lines}
)
INDEX="${{SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}}"
STATE="${{STATES[$INDEX]:?invalid array index}}"
RUN_DIR="$ROOT/states/$STATE"
cd "$RUN_DIR"
for required in INCAR POSCAR KPOINTS POTCAR wannier90.win; do test -s "$required" || {{ echo "missing $required" >&2; exit 2; }}; done
. /etc/profile.d/modules.sh
ulimit -s unlimited
module purge
module load vasp-cpu
command -v vasp_std
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
printf 'stage=w90_tb2j\\nstate=%s\\njob_id=%s\\narray_id=%s\\nhost=%s\\nstarted=%s\\n' "$STATE" "${{SLURM_JOB_ID:-unknown}}" "${{SLURM_ARRAY_TASK_ID:-unknown}}" "$(hostname)" "$(date --iso-8601=seconds)" > run-environment.txt
srun --ntasks="$SLURM_NTASKS" --cpu-bind=cores vasp_std > vasp.out 2>&1
printf 'vasp_finished=%s\\n' "$(date --iso-8601=seconds)" >> run-environment.txt
EFERMI="$($PYTHON -c 'from pymatgen.io.vasp.outputs import Vasprun; v=Vasprun("vasprun.xml",parse_dos=True,parse_eigen=False,parse_potcar_file=False); print(v.efermi)')"
test "$EFERMI" != "None"
printf 'efermi_eV=%s\\n' "$EFERMI" >> run-environment.txt
printf '%s\\n' "$PYTHON -m TB2J.scripts.wann2J --path $RUN_DIR --posfile $RUN_DIR/POSCAR --efermi $EFERMI --kmesh 9 9 9 --elements Ti Nb --prefix_up wannier90.1 --prefix_down wannier90.2 --output_path $RUN_DIR/tb2j_solve" > tb2j-command.txt
mkdir -p tb2j_solve
"$PYTHON" -m TB2J.scripts.wann2J --path "$RUN_DIR" --posfile "$RUN_DIR/POSCAR" --efermi "$EFERMI" --kmesh 9 9 9 --elements Ti Nb --prefix_up wannier90.1 --prefix_down wannier90.2 --output_path "$RUN_DIR/tb2j_solve" > tb2j.out 2>&1
printf 'tb2j_finished=%s\\n' "$(date --iso-8601=seconds)" >> run-environment.txt
'''

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--phase1-root',type=Path,default=Path('data/endmembers/BCC-Nb-Ti/magnetic-search-20260911')); args=ap.parse_args()
 out=args.output.resolve(); p1=args.phase1_root.resolve()
 if out.exists() and any(out.iterdir()): raise SystemExit(f'output must be new or empty: {out}')
 out.mkdir(parents=True); (out/'states').mkdir(); (out/'logs').mkdir()
 source=p1/'phase1'/'fm'
 srcs={x:source/x for x in ['POSCAR','POTCAR','KPOINTS']}
 for x,p in srcs.items():
  if not p.is_file(): raise SystemExit(f'missing source {p}')
 species=list(Poscar.from_file(srcs['POSCAR']).site_symbols)
 pot_titles=titles(srcs['POTCAR'])
 if species != ['Ti','Nb'] or pot_titles != ['PAW_PBE Ti_sv 07Sep2000','PAW_PBE Nb_sv 17Jan2003']:
  raise SystemExit(f'input identity mismatch species={species} titles={pot_titles}')
 kpoints=Kpoints.gamma_automatic(kpts=(8,8,8)); kp_path=out/'KPOINTS'; kpoints.write_file(kp_path)
 # Preserve exact phase1 POSCAR/POTCAR bytes; use fresh explicit KPOINTS.
 records=[]
 for state,mom in STATES.items():
  d=out/'states'/state; d.mkdir()
  for x,p in srcs.items(): shutil.copy2(p,d/x)
  kpoints.write_file(d/'KPOINTS')
  incar=dict(BASE); incar['MAGMOM']=mom; Incar(incar).write_file(d/'INCAR')
  win='''num_wann = 18\nnum_bands = 64\ndis_win_min = -50\ndis_win_max = 50\ndis_froz_min = -10\ndis_froz_max = 10\nbegin projections\nTi: s,p,d\nNb: s,p,d\nend projections\nwrite_hr = .true.\nwrite_xyz = .true.\n'''
  (d/'wannier90.win').write_text(win)
  records.append({'state':state,'site_order':species,'initial_moments_muB':mom,'input_sha256':{x:sha(d/x) for x in ['INCAR','POSCAR','KPOINTS','POTCAR','wannier90.win']}})
 w=out/'run_w90_tb2j_array.sbatch'; w.write_text(wrapper(out)); w.chmod(0o755)
 sub=out/'submit.sh'; sub.write_text(f'#!/usr/bin/env bash\nset -euo pipefail\nROOT="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"\nmkdir -p "$ROOT/logs"\nexec sbatch --parsable --array=0-5 "$ROOT/run_w90_tb2j_array.sbatch"\n'); sub.chmod(0o755)
 manifest={'schema_version':'nbti-w90-tb2j-rerun/v1','date':'2026-09-11','source_phase1':p1.as_posix(),'site_order':species,'potcar_titles':pot_titles,'recipe':BASE,'w90_kpoints':[8,8,8],'tb2j_kmesh':[9,9,9],'states':records,'safety':{'formal_tb2jflow_overwritten':False,'source_phase1_overwritten':False}}
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 (out/'README.md').write_text('Independent W90->TB2J rerun. Source formal tb2jflow is untouched. All six Phase-1 magnetic initial states are rerun with ICHARG=2 fresh SCF; TB2J is launched after VASP in the same Slurm array task.\n')
 print(json.dumps({'output':str(out),'states':list(STATES),'site_order':species,'potcar_titles':pot_titles,'wrapper':str(w)},indent=2))
if __name__=='__main__': main()
