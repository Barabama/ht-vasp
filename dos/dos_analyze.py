import json
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.plotter import DosPlotter
import matplotlib
matplotlib.use("Agg")  # 使用非交互式后端

def load_dos_data_from_json(json_file):
    """
    从JSON文件加载DOS数据
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # 提取DOS相关信息
    dos_info = data.get("vasp_objects", {}).get("dos", {})
    
    if not dos_info:
        print("在JSON文件中未找到DOS数据")
        return None
    
    # 从JSON数据重建CompleteDos对象
    # 由于CompleteDos对象不能直接从JSON重建，我们需要提取关键数据
    return data


def analyze_dos_data(json_file):
    """
    分析DOS数据并生成可视化
    """
    print("开始分析DOS数据...")
    
    # 加载JSON数据
    data = load_dos_data_from_json(json_file)
    if not data:
        return
    
    print(f"分析材料: {data.get('name', 'Unknown')}")
    print(f"化学系统: {data.get('chemsys', 'Unknown')}")
    print(f"体积: {data.get('volume', 'Unknown')} Å³")
    print(f"密度: {data.get('density', 'Unknown')} g/cm³")
    print(f"计算类型: {data.get('calc_type', 'Unknown')}")
    print(f"完成时间: {data.get('completed_at', 'Unknown')}")
    
    # 输出结构信息
    structure_info = data.get("structure", {})
    if structure_info:
        lattice_info = structure_info.get("lattice", {})
        if lattice_info:
            print(f"晶格参数: a={lattice_info.get('a', 'Unknown'):.3f} Å, "
                  f"b={lattice_info.get('b', 'Unknown'):.3f} Å, "
                  f"c={lattice_info.get('c', 'Unknown'):.3f} Å")
            print(f"角度: α={lattice_info.get('alpha', 'Unknown'):.2f}°, "
                  f"β={lattice_info.get('beta', 'Unknown'):.2f}°, "
                  f"γ={lattice_info.get('gamma', 'Unknown'):.2f}°")
    
    # 输出DOS信息
    output_info = data.get("output", {})
    if output_info:
        print(f"能量: {output_info.get('energy', 'Unknown')} eV")
        print(f"每原子能量: {output_info.get('energy_per_atom', 'Unknown')} eV")
        print(f"带隙: {output_info.get('bandgap', 'Unknown')} eV")
    
    # 从JSON中提取DOS数据
    vasp_objects = data.get("vasp_objects", {})
    dos_info = vasp_objects.get("dos", {})
    
    if dos_info:
        energies = dos_info.get("energies", [])
        total_dos = dos_info.get("densities", {}).get("1", [])  # 自旋向上
        spin_down_dos = dos_info.get("densities", {}).get("-1", [])  # 自旋向下
        
        if energies and total_dos:
            print(f"能量点数: {len(energies)}")
            print(f"Fermi能级: {dos_info.get('efermi', 'Unknown')} eV")
            
            # 计算费米能级附近的DOS值
            efermi = dos_info.get('efermi', 0)
            print(f"费米能级处的DOS: {np.interp(efermi, energies, total_dos):.4f}")
            
            # 绘制总态密度图
            plt.figure(figsize=(10, 6))
            plt.plot(energies, total_dos, label='Spin Up', color='blue')
            
            # 如果有自旋向下数据，也绘制它
            if spin_down_dos:
                plt.plot(energies, [-d for d in spin_down_dos], label='Spin Down', color='red')
            
            plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
            plt.axvline(x=efermi, color='r', linestyle='--', label=f'Fermi Level (E={efermi:.3f} eV)')
            plt.xlabel('Energy (eV)')
            plt.ylabel('Density of States')
            plt.title(f'Total Density of States - {data.get("name", "Unknown")}')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig('total_dos.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("总态密度图已保存为 total_dos.png")
            
            # 分析磁性性质
            if spin_down_dos:
                integrated_magnetization = np.trapz(
                    np.array(total_dos) - np.array(spin_down_dos), 
                    x=np.array(energies)
                )
                print(f"积分磁化强度 (费米能级以下): {integrated_magnetization:.4f} μB/atom")
    
    # 分析元素的DOS贡献（如果可用）
    pdos_info = dos_info.get("pdos", [])
    if pdos_info:
        print(f"找到 {len(pdos_info)} 个原子的PDOS数据")
        
        # 获取元素信息
        structure = dos_info.get("structure", {})
        if structure:
            sites = structure.get("sites", [])
            elements = [site.get("label", "Unknown") for site in sites]
            print(f"结构中的元素: {elements}")
    
    print("\n分析完成！")


def analyze_elements_dos(json_file):
    """
    分析各元素对DOS的贡献
    """
    print("\n开始分析各元素的DOS贡献...")
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    vasp_objects = data.get("vasp_objects", {})
    dos_info = vasp_objects.get("dos", {})
    
    if not dos_info:
        print("在JSON文件中未找到DOS数据")
        return
    
    energies = dos_info.get("energies", [])
    if not energies:
        print("未找到能量数据")
        return
    
    efermi = dos_info.get('efermi', 0)
    energies_relative = [e - efermi for e in energies]  # 相对于费米能级的能量
    
    pdos_data = dos_info.get("pdos", [])
    if not pdos_data:
        print("未找到分波态密度(PDOS)数据")
        return
    
    # 获取结构信息以确定元素
    structure = dos_info.get("structure", {})
    sites = structure.get("sites", [])
    elements = [site.get("label", "Unknown") for site in sites]
    
    # 检查第一个原子的PDOS结构
    if pdos_data and len(pdos_data) > 0:
        first_atom_pdos = pdos_data[0]
        print(f"第一个原子PDOS数据类型: {type(first_atom_pdos)}")
        if isinstance(first_atom_pdos, dict):
            print(f"第一个原子PDOS键: {list(first_atom_pdos.keys())}")
        elif isinstance(first_atom_pdos, list):
            print(f"第一个原子PDOS是列表，长度: {len(first_atom_pdos)}")
            if len(first_atom_pdos) > 0:
                print(f"第一个元素类型: {type(first_atom_pdos[0])}")
                if isinstance(first_atom_pdos[0], dict):
                    print(f"第一个元素键: {list(first_atom_pdos[0].keys())}")
    
    # 重新实现PDOS分析
    if pdos_data and elements:
        plt.figure(figsize=(12, 8))
        
        colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
        
        for idx, (atom_pdos, element) in enumerate(zip(pdos_data, elements)):
            color = colors[idx % len(colors)]
            
            # 如果是字典格式，直接使用
            if isinstance(atom_pdos, dict):
                for orbital_type, orbital_data in atom_pdos.items():
                    if isinstance(orbital_data, dict) and 'densities' in orbital_data:
                        densities = orbital_data['densities']
                        # 获取自旋向上数据
                        if '1' in densities:
                            orbital_dos = densities['1']
                            plt.plot(energies_relative, orbital_dos, 
                                   label=f'{element} {orbital_type}-orbital', 
                                   linestyle='-', color=color, alpha=0.7)
            # 如果是列表格式，按轨道类型处理
            elif isinstance(atom_pdos, list) and len(atom_pdos) > 0:
                # 尝试获取s, p, d轨道的DOS贡献
                orbital_names = ['s', 'p', 'd', 'f']
                for i, orbital_name in enumerate(orbital_names):
                    if i < len(atom_pdos):
                        orbital_data = atom_pdos[i]
                        if isinstance(orbital_data, dict) and 'densities' in orbital_data:
                            densities = orbital_data['densities']
                            if '1' in densities:
                                orbital_dos = densities['1']
                                plt.plot(energies_relative, orbital_dos, 
                                       label=f'{element} {orbital_name}-orbital', 
                                       linestyle='-', color=color, alpha=0.7)
        
        plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        plt.axvline(x=0, color='r', linestyle='--', label='Fermi Level')
        plt.xlabel('Energy - E_Fermi (eV)')
        plt.ylabel('Density of States')
        plt.title('Projected Density of States by Orbital')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('pdos_orbitals.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("分波态密度图已保存为 pdos_orbitals.png")


def main():
    """
    主函数，执行DOS数据分析
    """
    json_file = "data/Co-Cr/dos_result.json"
    
    # 执行基本DOS分析
    analyze_dos_data(json_file)
    
    # 执行元素DOS分析
    analyze_elements_dos(json_file)


if __name__ == "__main__":
    main()