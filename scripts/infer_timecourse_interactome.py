import os
import glob
import time
import logging
import requests
import numpy as np
import pandas as pd
from scipy import stats
import networkx as nx

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Constants
DATA_DIR = "data/transcriptomics/GSE247298"
RESULTS_DIR = "results/interactomes"
CYTOSCAPE_DIR = os.path.join(RESULTS_DIR, "cytoscape")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(CYTOSCAPE_DIR, exist_ok=True)

TIME_POINTS = ["12hrs", "24hrs", "48hrs", "72hrs", "120hrs"]
LOG2FC_THRESHOLD = 1.0
PVAL_THRESHOLD = 0.05
STRING_SCORE_THRESHOLD = 400
SPECIES = 9606 # Homo sapiens

# Interacciones Viral-Huésped Curadas por Literatura (Separación Temporal)
EARLY_RSV_INTERACTIONS = {
    "RSV-P": ["MAPK1", "MAPK3"],
    "RSV-NS1": ["STAT1", "STAT2", "MAVS"],
    "RSV-NS2": ["STAT2", "RIGI"]
}

LATE_RSV_INTERACTIONS = {
    "RSV-F": ["EGFR", "TLR4", "IGF1R"],
    "RSV-N": ["AKT1", "PIK3CA"],
    "RSV-M": ["NFKB1", "RELA"],
    "RSV-P": ["MAPK1", "MAPK3"] # P es esencial en ambas fases
}

def load_data():
    logger.info("Loading transcriptomics data...")
    # Load all NCF (healthy) samples
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
    merged_df['target_id_clean'] = merged_df['target_id'].apply(lambda x: x.split('.')[0])
    
    logger.info(f"Mapping {len(merged_df)} RefSeq IDs to Gene Symbols...")
    
    tpm_cols = [c for c in merged_df.columns if c not in ['target_id', 'target_id_clean']]
    expressed = merged_df[tpm_cols].max(axis=1) >= 1.0
    expressed_df = merged_df[expressed].copy()
    expressed_ids = expressed_df['target_id_clean'].tolist()
    
    gene_map = {}
    batch_size = 1000
    for i in range(0, len(expressed_ids), batch_size):
        batch = expressed_ids[i:i+batch_size]
        try:
            res = requests.post('https://mygene.info/v3/query', data={'q': ','.join(batch), 'scopes': 'refseq', 'fields': 'symbol', 'species': 'human'})
            if res.ok:
                results = res.json()
                for item in results:
                    if 'symbol' in item and 'query' in item:
                        gene_map[item['query']] = item['symbol']
        except Exception as e:
            logger.warning(f"Error mapping batch: {e}")
        time.sleep(0.5)
        
    expressed_df['gene_symbol'] = expressed_df['target_id_clean'].map(gene_map)
    expressed_df = expressed_df.dropna(subset=['gene_symbol'])
    
    grouped_df = expressed_df.groupby('gene_symbol')[tpm_cols].max().reset_index()
    logger.info(f"Loaded data for {len(grouped_df)} unique expressed genes.")
    
    return grouped_df

def calculate_differential_expression(df):
    logger.info("Calculating differential expression for all timepoints vs 0h...")
    control_cols = [c for c in df.columns if "0hrs" in c]
    
    for col in [c for c in df.columns if c != 'gene_symbol']:
        df[col] = np.log2(df[col] + 1)
        
    control_data = df[control_cols].values
    results = {}
    
    for tp in TIME_POINTS:
        tp_cols = [c for c in df.columns if tp in c]
        if not tp_cols:
            continue
            
        tp_data = df[tp_cols].values
        
        mean_control = np.mean(control_data, axis=1)
        mean_tp = np.mean(tp_data, axis=1)
        log2fc = mean_tp - mean_control
        
        with np.errstate(divide='ignore', invalid='ignore'):
            t_stat, p_val = stats.ttest_ind(tp_data, control_data, axis=1, equal_var=False)
            
        p_val = np.nan_to_num(p_val, nan=1.0)
        
        res_df = pd.DataFrame({
            'gene': df['gene_symbol'],
            'log2fc': log2fc,
            'pval': p_val
        })
        
        degs = res_df[(abs(res_df['log2fc']) >= LOG2FC_THRESHOLD) & (res_df['pval'] <= PVAL_THRESHOLD)]
        results[tp] = degs
        degs.to_csv(os.path.join(RESULTS_DIR, f"degs_{tp}.csv"), index=False)
        
    return results

def get_string_network(genes):
    if len(genes) < 2:
        return pd.DataFrame()
        
    genes = list(genes)
    if len(genes) > 2000:
        genes = genes[:2000]
        
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(genes),
        "species": SPECIES,
        "required_score": STRING_SCORE_THRESHOLD,
        "caller_identity": "RSV_Interactome_Script"
    }
    
    try:
        res = requests.post(url, data=params)
        if res.ok and res.text.strip():
            df = pd.read_json(res.text)
            if not df.empty:
                return df[['preferredName_A', 'preferredName_B', 'score']]
    except Exception as e:
        logger.warning(f"Error fetching STRING network: {e}")
        
    return pd.DataFrame()

def analyze_topology(network_df, timepoint):
    if network_df.empty:
        return pd.DataFrame()
        
    G = nx.from_pandas_edgelist(network_df, 'preferredName_A', 'preferredName_B', ['score'])
    degree_dict = dict(G.degree())
    betweenness_dict = nx.betweenness_centrality(G)
    
    df = pd.DataFrame({
        'node': list(G.nodes()),
        'degree': [degree_dict.get(n, 0) for n in G.nodes()],
        'betweenness': [betweenness_dict.get(n, 0) for n in G.nodes()]
    })
    
    df = df.sort_values('degree', ascending=False)
    df.to_csv(os.path.join(RESULTS_DIR, f"topology_{timepoint}.csv"), index=False)
    return df

def export_to_cytoscape(network_df, degs_df, tp):
    """
    Calculates the dynamic interaction probability and exports .sif and attributes.
    """
    if network_df.empty:
        return
        
    logger.info(f"  Calculating dynamic interaction probabilities for {tp} (SIF Export)...")
    
    # Pre-calculate RNA weights (Log2FC normalized to 0-1 confidence loosely)
    # We use min-max scaling of absolute log2FC as a simple weight
    abs_logfc = degs_df.set_index('gene')['log2fc'].abs()
    if abs_logfc.max() > 0:
        rna_weights = abs_logfc / abs_logfc.max()
    else:
        rna_weights = abs_logfc
        
    edges_sif = []
    edge_attrs = ["Edge\tInteractionProbability\tInteractionType"]
    
    for _, row in network_df.iterrows():
        nodeA = row['preferredName_A']
        nodeB = row['preferredName_B']
        string_score = row['score'] # 0.0 to 1.0 (STRING API already returns 0-1 scale)
        
        # Calculate RNA Weight
        rnaA = rna_weights.get(nodeA, 0.5)
        rnaB = rna_weights.get(nodeB, 0.5)
        
        # Combined Node Weights (RNA only)
        weightA = rnaA
        weightB = rnaB
        
        # Interaction Probability
        prob = string_score * ((weightA + weightB) / 2.0)
        
        # Format for SIF
        edges_sif.append(f"{nodeA}\tpp\t{nodeB}")
        
        # Format for Attributes
        edge_name = f"{nodeA} (pp) {nodeB}"
        edge_attrs.append(f"{edge_name}\t{prob:.4f}\tDynamic_Proteogenomic")
        
    # Write SIF
    sif_path = os.path.join(CYTOSCAPE_DIR, f"network_{tp}.sif")
    with open(sif_path, 'w') as f:
        f.write("\n".join(edges_sif) + "\n")
        
    # Write Attributes
    attr_path = os.path.join(CYTOSCAPE_DIR, f"edge_attributes_{tp}.txt")
    with open(attr_path, 'w') as f:
        f.write("\n".join(edge_attrs) + "\n")
        
    logger.info(f"  -> Exported {len(edges_sif)} interactions to {sif_path}")

def main():
    df = load_data()
    degs_by_tp = calculate_differential_expression(df)
    summary_data = []
    
    for tp in TIME_POINTS:
        if tp not in degs_by_tp:
            continue
            
        logger.info(f"Building interactome for {tp}...")
        degs_df = degs_by_tp[tp]
        degs_list = degs_df['gene'].tolist()
        
        network_df = get_string_network(degs_list)
        logger.info(f"  Network for {tp} (Human-Human): {len(network_df)} interactions")
        
        # Inject Viral-Host Interactions based on timepoint
        current_interactions = EARLY_RSV_INTERACTIONS if tp in ["12hrs", "24hrs"] else LATE_RSV_INTERACTIONS
        
        viral_edges = []
        for v_prot, h_targets in current_interactions.items():
            for h_prot in h_targets:
                # Inyección forzada: muchos blancos virales (como quinasas) no son DEGs, pero son atacados.
                viral_edges.append({
                    'preferredName_A': v_prot,
                    'preferredName_B': h_prot,
                    'score': 1.0  # Máxima confianza para interacciones validadas
                })
        
        if viral_edges:
            viral_df = pd.DataFrame(viral_edges)
            network_df = pd.concat([network_df, viral_df], ignore_index=True)
            logger.info(f"  Inyectadas {len(viral_edges)} interacciones VSR-Huésped")
            
        if not network_df.empty:
            network_df.to_csv(os.path.join(RESULTS_DIR, f"network_{tp}.csv"), index=False)
            topo_df = analyze_topology(network_df, tp)
            
            # Export to Cytoscape (SIF) with probabilities
            export_to_cytoscape(network_df, degs_df, tp)
            
            top_hubs = topo_df.head(5)['node'].tolist() if not topo_df.empty else []
            summary_data.append({
                'Timepoint': tp,
                'DEGs': len(degs_list),
                'Interactions': len(network_df),
                'Top_Hubs': ', '.join(top_hubs)
            })
            
    if summary_data:
        pd.DataFrame(summary_data).to_csv(os.path.join(RESULTS_DIR, "interactome_summary.csv"), index=False)
        logger.info("Interactome analysis complete! Check the results/interactomes directory.")

if __name__ == "__main__":
    main()
