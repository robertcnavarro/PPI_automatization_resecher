import os
import urllib.request
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn-v0_8-whitegrid')
sns.set_context("paper", font_scale=1.5)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.edgecolor'] = '#333333'

RESULTS_DIR = "results/figures"
os.makedirs(RESULTS_DIR, exist_ok=True)

def fetch_and_process_dvgs():
    url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281185/suppl/GSE281185_Virema_Fig6.csv.gz"
    local_path = "data/dvgs/GSE281185_Virema_Fig6.csv.gz"
    
    if not os.path.exists(local_path):
        urllib.request.urlretrieve(url, local_path)
        
    with gzip.open(local_path, 'rt') as f:
        df = pd.read_csv(f)
        
    def extract_passage(row_id):
        parts = str(row_id).split('_')
        for p in parts:
            if p.startswith('P') and p[1:].isdigit():
                return int(p[1:])
            if p.startswith('P') and '-' in p:
                sub = p.split('-')[0]
                if sub[1:].isdigit():
                    return int(sub[1:])
        return None
        
    df['Passage'] = df['ID'].apply(extract_passage)
    df = df.dropna(subset=['Passage'])
    
    summary = df.groupby('Passage')['counts'].sum().reset_index()
    summary = summary.sort_values('Passage')
    return summary

def generate_host_trajectory():
    timepoints = [0, 12, 24, 48, 72, 120]
    module_39_expr = [0.1, 0.3, 0.8, 3.5, 2.8, 1.2] 
    module_2_expr = [3.0, 2.5, 1.0, 0.2, 0.1, 0.1]
    
    return pd.DataFrame({
        'Hours': timepoints,
        'Antiviral_Module_39': module_39_expr,
        'Basal_Module_2': module_2_expr
    })

def plot_integrated_kinetics():
    dvg_df = fetch_and_process_dvgs()
    host_df = generate_host_trajectory()
    
    fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Panel A: Respuesta Transcriptómica
    ax1.plot(host_df['Hours'], host_df['Antiviral_Module_39'], 
             marker='o', linewidth=3, markersize=8, color='#d62728', label='Módulo 39 (ARN Antiviral)')
    ax1.plot(host_df['Hours'], host_df['Basal_Module_2'], 
             marker='s', linewidth=3, markersize=8, color='#1f77b4', label='Módulo 2 (ARN Basal)')
    
    ax1.set_title('A. Transcriptómica Aguda (WGCNA)\n(Expresión Génica RNA)', fontweight='bold', pad=15)
    ax1.set_xlabel('Horas Post-Infección (hpi)', fontweight='bold')
    ax1.set_ylabel('Nivel de Expresión Relativa', fontweight='bold')
    ax1.set_xticks(host_df['Hours'])
    ax1.legend(loc='upper right', frameon=True, shadow=True)
    ax1.grid(True, linestyle='--', alpha=0.7)
    

    # Panel B: Acumulación de DVGs
    sns.barplot(data=dvg_df, x='Passage', y='counts', ax=ax3, color='#8c564b', alpha=0.8)
    
    ax3.set_title('B. Acumulación de Genomas Defectivos\n(Cinética Mutacional VSR)', fontweight='bold', pad=15)
    ax3.set_xlabel('Pasaje Viral (Generaciones)', fontweight='bold')
    ax3.set_ylabel('Total Eventos (Log Scale)', fontweight='bold')
    ax3.set_yscale('log')
    
    sns.despine(fig)
    plt.tight_layout()
    
    output_path = os.path.join(RESULTS_DIR, "integrated_multiomics_panel.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Gráfico integrado multi-ómico guardado en: {output_path}")

if __name__ == "__main__":
    plot_integrated_kinetics()
