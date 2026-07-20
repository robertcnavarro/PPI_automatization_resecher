import os
import gzip
import logging
import urllib.request
import tarfile
import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

DATA_DIR = "data/scrnaseq/GSE281623"
RESULTS_DIR = "results/scrnaseq"
FIGURES_DIR = "results/figures"
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

FTP_BASE = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281623/suppl/"
FILE_TO_DOWNLOAD = "GSE281623_RAW.tar"
PREFIX_TP4 = "GSM8625035_A549_Timepoint_4_08262019"

def download_and_extract():
    local_path = os.path.join(DATA_DIR, FILE_TO_DOWNLOAD)
    if not os.path.exists(local_path):
        logger.info(f"Descargando {FILE_TO_DOWNLOAD}...")
        urllib.request.urlretrieve(FTP_BASE + FILE_TO_DOWNLOAD, local_path)
        
    # Check if extracted
    matrix_path = os.path.join(DATA_DIR, f"{PREFIX_TP4}_matrix.mtx.gz")
    if not os.path.exists(matrix_path):
        logger.info("Extrayendo archivo TAR...")
        with tarfile.open(local_path) as tar:
            tar.extractall(path=DATA_DIR)
    
    return matrix_path

def load_scrnaseq_data():
    matrix_path = os.path.join(DATA_DIR, f"{PREFIX_TP4}_matrix.mtx.gz")
    features_path = os.path.join(DATA_DIR, f"{PREFIX_TP4}_features.tsv.gz")
    barcodes_path = os.path.join(DATA_DIR, f"{PREFIX_TP4}_barcodes.tsv.gz")
    
    logger.info("Cargando matriz dispersa (.mtx.gz)... (Esto puede tomar unos segundos)")
    matrix = scipy.io.mmread(matrix_path).tocsr()
    
    logger.info("Cargando genes y barcodes...")
    features = pd.read_csv(features_path, sep='\t', header=None, names=['Ensembl', 'Symbol', 'Type'])
    features['Clean_Symbol'] = features['Symbol'].apply(lambda x: x.split('_')[1] if isinstance(x, str) and '_' in x else str(x))
    
    barcodes = pd.read_csv(barcodes_path, sep='\t', header=None, names=['Barcode'])
    
    return matrix, features, barcodes

def analyze_cell_communication(matrix, features, barcodes):
    logger.info("Normalizando matriz (Log1p)...")
    matrix_log = matrix.copy()
    matrix_log.data = np.log1p(matrix_log.data)
    
    # Transponer para tener Células en Filas y Genes en Columnas
    matrix_log_t = matrix_log.T
    
    logger.info("Reduciendo dimensionalidad (SVD / PCA)...")
    svd = TruncatedSVD(n_components=10, random_state=42)
    pca_result = svd.fit_transform(matrix_log_t)
    
    logger.info("Clustering (K-Means) para identificar 4 sub-poblaciones...")
    kmeans = KMeans(n_clusters=4, random_state=42)
    clusters = kmeans.fit_predict(pca_result)
    
    # Base de datos simplificada de Ligandos y Receptores críticos en VSR
    lr_pairs = [
        ("CCL2", "CCR2"),
        ("CXCL10", "CXCR3"),
        ("IFNB1", "IFNAR1"),
        ("IFNL1", "IFNLR1"),
        ("IL6", "IL6R"),
        ("TNF", "TNFRSF1A")
    ]
    
    logger.info("Calculando Interacciones Célula-Célula (CCI)...")
    gene_symbols = features['Clean_Symbol'].tolist()
    cci_results = []
    
    target_genes = set([p[0] for p in lr_pairs] + [p[1] for p in lr_pairs])
    target_indices = {g: i for i, g in enumerate(gene_symbols) if g in target_genes}
    
    if not target_indices:
        logger.warning("No se encontraron genes diana en la matriz.")
        return
        
    cluster_means = {}
    for c in range(4):
        cell_idx = np.where(clusters == c)[0]
        sub_matrix = matrix_log_t[cell_idx, :]
        means = np.array(sub_matrix.mean(axis=0)).flatten()
        cluster_means[c] = means
        
    for ligand, receptor in lr_pairs:
        if ligand not in target_indices or receptor not in target_indices:
            continue
            
        l_idx = target_indices[ligand]
        r_idx = target_indices[receptor]
        
        for cA in range(4): # Emisor
            for cB in range(4): # Receptor
                expr_L = cluster_means[cA][l_idx]
                expr_R = cluster_means[cB][r_idx]
                
                score = expr_L * expr_R
                
                cci_results.append({
                    'Ligand_Receptor': f"{ligand} -> {receptor}",
                    'Sender_Cluster': f"Cluster {cA}",
                    'Receiver_Cluster': f"Cluster {cB}",
                    'Score': score
                })
                
    cci_df = pd.DataFrame(cci_results)
    
    if cci_df.empty:
        logger.warning("No se pudieron inferir interacciones con los pares dados.")
        return
        
    cci_df.to_csv(os.path.join(RESULTS_DIR, "cell_cell_interactions.csv"), index=False)
    
    logger.info("Generando Heatmap de Comunicación...")
    # Sumarizamos la fuerza de todas las interacciones activas
    heatmap_data = cci_df.pivot_table(index='Sender_Cluster', columns='Receiver_Cluster', values='Score', aggfunc='sum')
    
    # Evitar NaN
    heatmap_data = heatmap_data.fillna(0)
    
    # Normalizar los datos al rango [0, 1] para que sea fácil de interpretar visualmente
    max_val = heatmap_data.values.max()
    if max_val > 0:
        heatmap_data = heatmap_data / max_val
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(heatmap_data, annot=True, cmap="YlOrRd", fmt=".3f", linewidths=.5)
    plt.title("Intensidad Relativa de Comunicación Célula-Célula\n(scRNA-seq Timepoint Tardío)")
    plt.xlabel("Sub-población Receptora")
    plt.ylabel("Sub-población Emisora")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "cell_cell_interactions.png"), dpi=300)
    plt.close()
    
    logger.info(f"Análisis scRNA-seq finalizado. Heatmap guardado en {FIGURES_DIR}/cell_cell_interactions.png")

def main():
    download_and_extract()
    matrix, features, barcodes = load_scrnaseq_data()
    analyze_cell_communication(matrix, features, barcodes)

if __name__ == "__main__":
    main()
