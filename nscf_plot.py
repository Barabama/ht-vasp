import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from monty.json import MontyDecoder
from pymatgen.electronic_structure.dos import CompleteDos, Spin, Dos
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter, BSDOSPlotter

dos_path = Path("/nfs_hdd/2025/gaominliang/ht-vasp/data/CoNiOH/dos_out.json")
band_path = Path("/nfs_hdd/2025/gaominliang/ht-vasp/data/CoNiOH/band_out.json")
output_dir = Path("/nfs_hdd/2025/gaominliang/ht-vasp/data/CoNiOH")
output_dir.mkdir(exist_ok=True)

print("=" * 60)
print("1. DOS 数据分析")
print("=" * 60)

with open(dos_path, "r", encoding="utf-8") as f:
    dos_data = json.load(f, cls=MontyDecoder)

dos: CompleteDos = dos_data.vasp_objects["dos"]

print(f"\n[DOS 基本信息]")
print(f"  费米能级: {dos.efermi:.4f} eV")
print(f"  能量范围: {dos.energies.min():.4f} ~ {dos.energies.max():.4f} eV")
print(f"  能量点数: {len(dos.energies)}")
is_spin_polarized = Spin.down in dos.densities
print(f"  自旋极化: {is_spin_polarized}")

print(f"\n[结构信息]")
print(f"  化学式: {dos.structure.formula}")
print(f"  元素: {[str(el) for el in dos.structure.elements]}")

print(f"\n[DOS 数据类型]")
print(f"  energies: {type(dos.energies)}, shape: {dos.energies.shape}")
print(f"  densities: {type(dos.densities)}")
print(f"  pdos: {type(dos.pdos)}")

print(f"\n[可用方法]")
dos_methods = [m for m in dir(dos) if not m.startswith('_') and callable(getattr(dos, m))]
print(f"  {', '.join(dos_methods[:15])} ...")

print(f"\n[带隙信息]")
try:
    bandgap = dos.get_gap()
    print(f"  带隙: {bandgap:.4f} eV")
except Exception as e:
    print(f"  无法计算带隙: {e}")

print(f"\n[元素投影 DOS]")
for el in dos.structure.elements:
    el_dos = dos.get_element_dos()[el]
    print(f"  {el}: 可用")

print(f"\n[轨道投影 DOS (s, p, d)]")
spd_dos = dos.get_spd_dos()
for orb, orb_dos in spd_dos.items():
    print(f"  {orb.name}: 可用")

print("\n" + "=" * 60)
print("2. Band 数据分析")
print("=" * 60)

with open(band_path, "r", encoding="utf-8") as f:
    band_data = json.load(f, cls=MontyDecoder)

band: BandStructureSymmLine = band_data.vasp_objects["bandstructure"]

print(f"\n[Band 基本信息]")
print(f"  费米能级: {band.efermi:.4f} eV")
print(f"  自旋极化: {band.is_spin_polarized}")
print(f"  是否金属: {band.is_metal()}")

print(f"\n[能带数据]")
print(f"  k点数: {len(band.kpoints)}")
print(f"  能带数: {band.nb_bands}")

print(f"\n[带隙信息]")
bg = band.get_band_gap()
print(f"  带隙: {bg['energy']:.4f} eV")
print(f"  直接带隙: {bg['direct']}")
print(f"  跃迁: {bg['transition']}")

print(f"\n[VBM/CBM]")
vbm = band.get_vbm()
cbm = band.get_cbm()
print(f"  VBM 能量: {vbm['energy']:.4f} eV")
print(f"  VBM k点: {vbm['kpoint_index']}")
print(f"  CBM 能量: {cbm['energy']:.4f} eV")
print(f"  CBM k点: {cbm['kpoint_index']}")

print(f"\n[高对称点标签]")
for label, kpoint in band.labels_dict.items():
    print(f"  {label}: {kpoint}")

print(f"\n[可用方法]")
band_methods = [m for m in dir(band) if not m.startswith('_') and callable(getattr(band, m))]
print(f"  {', '.join(band_methods[:15])} ...")

print("\n" + "=" * 60)
print("3. 绘制 DOS 图 (使用 DosPlotter)")
print("=" * 60)

plotter = DosPlotter()
plotter.add_dos("Total", dos)
ax = plotter.get_plot()
ax.set_title("Total DOS", fontsize=14)
ax.axvline(dos.efermi, color='k', linestyle='--', linewidth=0.8, label='E_F')
ax.legend()
ax.grid(True, alpha=0.3)

dos_plot_path = output_dir / "dos_total.png"
plt.savefig(dos_plot_path, dpi=300, bbox_inches='tight')
print(f"\n总 DOS 图已保存: {dos_plot_path}")
plt.close()

plotter2 = DosPlotter()
for el in dos.structure.elements:
    el_dos = dos.get_element_dos()[el]
    plotter2.add_dos(str(el), el_dos)
ax2 = plotter2.get_plot()
ax2.set_title("Element-projected DOS", fontsize=14)
ax2.axvline(dos.efermi, color='k', linestyle='--', linewidth=0.8, label='E_F')
ax2.legend()
ax2.grid(True, alpha=0.3)

dos_plot_path2 = output_dir / "dos_element.png"
plt.savefig(dos_plot_path2, dpi=300, bbox_inches='tight')
print(f"元素投影 DOS 图已保存: {dos_plot_path2}")
plt.close()

plotter3 = DosPlotter()
spd_dos = dos.get_spd_dos()
for orb, orb_dos in spd_dos.items():
    plotter3.add_dos(orb.name, orb_dos)
ax3 = plotter3.get_plot()
ax3.set_title("Orbital-projected DOS (s, p, d)", fontsize=14)
ax3.axvline(dos.efermi, color='k', linestyle='--', linewidth=0.8, label='E_F')
ax3.legend()
ax3.grid(True, alpha=0.3)

dos_plot_path3 = output_dir / "dos_orbital.png"
plt.savefig(dos_plot_path3, dpi=300, bbox_inches='tight')
print(f"轨道投影 DOS 图已保存: {dos_plot_path3}")
plt.close()

if is_spin_polarized:
    plotter4 = DosPlotter()
    spin_up_dos = Dos(dos.efermi, dos.energies, {Spin.up: dos.densities[Spin.up]})
    spin_down_dos = Dos(dos.efermi, dos.energies, {Spin.up: dos.densities[Spin.down]})
    plotter4.add_dos("Spin Up", spin_up_dos)
    plotter4.add_dos("Spin Down", spin_down_dos)
    ax4 = plotter4.get_plot()
    ax4.set_title("Spin-polarized DOS", fontsize=14)
    ax4.axvline(dos.efermi, color='k', linestyle='--', linewidth=0.8, label='E_F')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    dos_plot_path4 = output_dir / "dos_spin.png"
    plt.savefig(dos_plot_path4, dpi=300, bbox_inches='tight')
    print(f"自旋极化 DOS 图已保存: {dos_plot_path4}")
    plt.close()

print("\n" + "=" * 60)
print("4. 绘制能带图")
print("=" * 60)

from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine

if isinstance(band, BandStructureSymmLine) and len(band.labels_dict) > 0:
    bs_plotter = BSPlotter(band)
    plt_bs = bs_plotter.get_plot()
    plt_bs.gca().set_title("Band Structure", fontsize=14)

    band_plot_path = output_dir / "band_plot.png"
    plt_bs.savefig(band_plot_path, dpi=300, bbox_inches='tight')
    print(f"\n能带图已保存: {band_plot_path}")
    plt.close()
else:
    print("\n注意: band_out.json 是均匀网格上的能带数据 (NSCF Uniform)")
    print("      需要沿高对称线的能带数据 (NSCF Line) 才能绘制能带图")
    print("      跳过能带图绘制...")

    fig, ax = plt.subplots(figsize=(10, 6))
    for i in range(min(10, band.nb_bands)):
        if band.is_spin_polarized:
            energies = band.bands[Spin.up][i] - band.efermi
        else:
            energies = band.bands[Spin.up][i] - band.efermi
        ax.plot(range(len(energies)), energies, 'b-', linewidth=0.5, alpha=0.7)
    ax.set_xlabel("K-point index", fontsize=12)
    ax.set_ylabel("Energy - E_F (eV)", fontsize=12)
    ax.set_title("Band Structure (Uniform Grid)", fontsize=14)
    ax.axhline(0, color='k', linestyle='--', linewidth=0.8, label='E_F')
    ax.legend()
    ax.grid(True, alpha=0.3)

    band_plot_path = output_dir / "band_uniform.png"
    plt.savefig(band_plot_path, dpi=300, bbox_inches='tight')
    print(f"\n均匀网格能带图已保存: {band_plot_path}")
    plt.close()

print("\n" + "=" * 60)
print("5. 绘制 DOS + Band 联合图")
print("=" * 60)

if isinstance(band, BandStructureSymmLine) and len(band.labels_dict) > 0:
    bsdos_plotter = BSDOSPlotter()
    plt_combined = bsdos_plotter.get_plot(band, dos)

    combined_plot_path = output_dir / "band_dos_combined.png"
    plt_combined.savefig(combined_plot_path, dpi=300, bbox_inches='tight')
    print(f"\n联合图已保存: {combined_plot_path}")
    plt.close()
else:
    print("\n注意: 需要 BandStructureSymmLine 才能绘制联合图")
    print("      跳过联合图绘制...")

print("\n" + "=" * 60)
print("完成！所有图像已保存到:")
print(f"  {output_dir}")
print("=" * 60)
