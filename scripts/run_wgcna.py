import os
import glob
import logging
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.decomposition import PCA

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

DATA_DIR = "data/transcriptomics/GSE247298"
RESULTS_DIR = "results/wgcna"
os.makedirs(RESULTS_DIR, exist_ok=True)

TIME_POINTS = ["0hrs", "12hrs", "24hrs", "48hrs", "72hrs", "120hrs"]

def load_preprocessed_data():
    """Load the preprocessed expression data we created earlier or re-load the raw data"""
    logger.info("Cargando datos de expresión para WGCNA...")
    files = glob.glob(os.path.join(DATA_DIR, "*_NCF*.tsv.gz"))
    
    data_dict = {}
    for f in files:
        basename = os.path.basename(f)
        parts = basename.replace(".tsv.gz", "").split("_")
        subject = parts[2]
        timepoint = parts[3]
        
        df = pd.read_csv(f, sep='\t')
        df = df[['target_id', 'tpm']]
        df = df.rename(columns={'tpm': f"{subject}_{timepoint}"})
        
        if len(data_dict) == 0:
            data_dict['target_id'] = df['target_id']
            
        data_dict[f"{subject}_{timepoint}"] = df[f"{subject}_{timepoint}"]
        
    merged_df = pd.DataFrame(data_dict)
    
    # Filter for top 5000 most variable genes to make WGCNA computationally feasible in pure python
    numeric_cols = [c for c in merged_df.columns if c != 'target_id']
    variances = merged_df[numeric_cols].var(axis=1)
    
    # Get top 2000 most variable genes
    top_indices = variances.nlargest(2000).index
    filtered_df = merged_df.loc[top_indices]
    
    # Log2 transform
    for col in numeric_cols:
        filtered_df[col] = np.log2(filtered_df[col] + 1)
        
    filtered_df.set_index('target_id', inplace=True)
    return filtered_df

def construct_coexpression_network(df):
    logger.info("Calculando matriz de correlación (WGCNA step 1)...")
    # Compute correlation matrix (Genes x Genes)
    corr_matrix = df.T.corr(method='pearson')
    
    # Convert correlation to adjacency (soft thresholding approximation: power of 6)
    # Using absolute correlation to group both positive and negative correlated genes,
    # or just positive. Standard WGCNA uses signed or unsigned. We'll use unsigned:
    adjacency = np.power(np.abs(corr_matrix), 6)
    
    # Convert to Topological Overlap Matrix (TOM) equivalent - Distance matrix
    # TOM is complex, so we use a simplified correlation distance
    logger.info("Calculando matriz de distancias y clustering jerárquico...")
    distance_matrix = 1 - corr_matrix
    
    # Hierarchical clustering
    # We use condensed distance matrix for scipy linkage
    import scipy.spatial.distance as ssd
    condensed_dist = ssd.squareform(distance_matrix)
    Z = linkage(condensed_dist, method='average')
    
    # Cut tree to form modules (height parameter)
    # 0.5 is an arbitrary height for demonstration
    labels = fcluster(Z, t=0.5, criterion='distance')
    
    module_assignments = pd.DataFrame({
        'Gene': df.index,
        'Module': labels
    })
    
    unique_modules = module_assignments['Module'].nunique()
    logger.info(f"Se identificaron {unique_modules} módulos de coexpresión.")
    
    module_assignments.to_csv(os.path.join(RESULTS_DIR, "module_assignments.csv"), index=False)
    return module_assignments, Z

def calculate_eigengenes(df, module_assignments):
    logger.info("Calculando Module Eigengenes (1er Componente Principal por módulo)...")
    modules = module_assignments['Module'].unique()
    
    eigengenes = {}
    
    for mod in modules:
        genes_in_mod = module_assignments[module_assignments['Module'] == mod]['Gene']
        mod_data = df.loc[genes_in_mod].T # Samples x Genes
        
        # PCA to get first component
        if mod_data.shape[1] > 1:
            pca = PCA(n_components=1)
            first_pc = pca.fit_transform(mod_data)
            eigengenes[f"Module_{mod}"] = first_pc[:, 0]
        else:
            eigengenes[f"Module_{mod}"] = mod_data.values[:, 0]
            
    eigengenes_df = pd.DataFrame(eigengenes, index=df.columns)
    
    # Create trait matrix (timepoints)
    # Convert timepoints to numerical values for correlation (0, 12, 24, 48, 72, 120)
    time_values = []
    for col in eigengenes_df.index:
        tp_str = col.split('_')[1].replace('hrs', '')
        time_values.append(int(tp_str))
        
    eigengenes_df['Time'] = time_values
    
    # Correlate modules with Time
    correlations = []
    for mod in modules:
        mod_name = f"Module_{mod}"
        r, p = pearsonr(eigengenes_df[mod_name], eigengenes_df['Time'])
        correlations.append({
            'Module': mod_name,
            'Correlation_with_Time': r,
            'P_value': p,
            'Size': len(module_assignments[module_assignments['Module'] == mod])
        })
        
    corr_df = pd.DataFrame(correlations)
    corr_df = corr_df.sort_values(by='Correlation_with_Time', key=abs, ascending=False)
    
    corr_df.to_csv(os.path.join(RESULTS_DIR, "module_time_correlations.csv"), index=False)
    logger.info("Correlación de módulos con el tiempo de infección completada.")
    
    # Print top 5 modules
    print("\nTop 5 Módulos correlacionados con el avance temporal de la infección:")
    print(corr_df.head(5).to_string(index=False))

def main():
    df = load_preprocessed_data()
    module_assignments, Z = construct_coexpression_network(df)
    calculate_eigengenes(df, module_assignments)
    
if __name__ == "__main__":
    main()
