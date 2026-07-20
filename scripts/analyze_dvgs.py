import os
import urllib.request
import gzip
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

DATA_DIR = "data/dvgs/GSE281185"
RESULTS_DIR = "results/dvgs"
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

FTP_BASE = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281185/suppl/"
FILE_TO_DOWNLOAD = "GSE281185_ViReMa_R2-8U.csv.gz"

def download_data():
    local_path = os.path.join(DATA_DIR, FILE_TO_DOWNLOAD)
    if not os.path.exists(local_path):
        logger.info(f"Descargando {FILE_TO_DOWNLOAD} desde GEO FTP...")
        urllib.request.urlretrieve(FTP_BASE + FILE_TO_DOWNLOAD, local_path)
    return local_path

def analyze_dvgs(filepath):
    logger.info("Analizando cinética de Genomas Defectivos Virales (DVGs)...")
    
    # Read compressed CSV
    with gzip.open(filepath, 'rt') as f:
        df = pd.read_csv(f)
        
    logger.info(f"Columnas disponibles en ViReMa output: {list(df.columns)}")
    logger.info(f"Total de registros de recombinación viral: {len(df)}")
    
    # Simple summary of the data
    # ViReMa usually outputs 'Start', 'Stop', 'Defect_Type', 'Read_Count' etc.
    # We will generate a basic summary if 'Read_Count' or similar exists.
    
    # Check if we can group by sample or type
    if 'Sample' in df.columns:
        summary = df.groupby('Sample').size().reset_index(name='Total_Recombination_Events')
        summary.to_csv(os.path.join(RESULTS_DIR, "dvg_summary_by_sample.csv"), index=False)
        logger.info("Resumen guardado en results/dvgs/dvg_summary_by_sample.csv")
    else:
        # Just save the first 100 rows to see what it looks like
        df.head(100).to_csv(os.path.join(RESULTS_DIR, "dvg_preview.csv"), index=False)
        logger.info("Vista previa guardada en results/dvgs/dvg_preview.csv (No se encontró columna 'Sample')")
        
    return df

def main():
    filepath = download_data()
    analyze_dvgs(filepath)
    
if __name__ == "__main__":
    main()
