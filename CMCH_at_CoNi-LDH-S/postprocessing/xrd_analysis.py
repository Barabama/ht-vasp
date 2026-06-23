"""
XRD 分析 — 从实验数据文件中提取峰位
"""
import argparse, math, warnings
import numpy as np
from pathlib import Path

warnings.filterwarnings("ignore")

def read_xrd(path):
    """读取 XRD 数据文件，返回 (2theta, intensity) 数组"""
    with open(path) as f:
        lines = f.readlines()
    # 跳过 header
    data = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) == 2:
            try:
                t = float(parts[0])
                i = float(parts[1])
                data.append((t, i))
            except:
                continue
    arr = np.array(data)
    return arr[:, 0], arr[:, 1]


def find_peaks(two_theta, intensity, min_prominence=5, min_height=10):
    """用滑动窗口找峰"""
    peaks = []
    n = len(intensity)
    # 窗口宽度 ~0.5°
    window = max(3, int(0.5 / (two_theta[1] - two_theta[0])))

    for i in range(window, n - window):
        if intensity[i] < min_height:
            continue
        left = intensity[i-window:i].max()
        right = intensity[i+1:i+window+1].max()
        if intensity[i] > left and intensity[i] > right:
            prom = intensity[i] - max(left, right)
            if prom >= min_prominence:
                peaks.append((two_theta[i], intensity[i], prom))

    peaks.sort(key=lambda x: x[1], reverse=True)
    return peaks


def _d_spacing(two_theta, wavelength=1.5406):
    """Bragg 方程: d = λ/(2sinθ)"""
    theta_rad = np.deg2rad(two_theta / 2)
    return wavelength / (2 * np.sin(theta_rad))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-prominence", type=float, default=5)
    parser.add_argument("--min-height", type=float, default=10)
    args = parser.parse_args()

    d1 = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S")
    for fname in ["s1-PH.txt", "s3-PH.txt"]:
        tt, inten = read_xrd(d1 / fname)
        peaks = find_peaks(tt, inten, args.min_prominence, args.min_height)

        print(f"\n{'='*60}")
        print(f"  File: {fname}  ({len(tt)} points)")
        print(f"  Peaks found: {len(peaks)} (prominence ≥ {args.min_prominence})")
        print(f"{'='*60}")
        print(f"  {'#':>3s}  {'2θ (°)':>8s}  {'Intensity':>10s}  {'Prom':>6s}  {'d (Å)':>8s}")
        print(f"  {'─'*3}  {'─'*8}  {'─'*10}  {'─'*6}  {'─'*8}")
        for idx, (tth, iten, prom) in enumerate(peaks[:25]):
            d = _d_spacing(tth)
            print(f"  {idx+1:>3d}  {tth:>8.3f}  {iten:>10.1f}  {prom:>6.1f}  {d:>8.4f}")
