import time
import logging
from typing import List, Dict, Any
from urllib.error import HTTPError
import pandas as pd
from Bio import Entrez

# 1. Configuración y Constantes
ENTREZ_EMAIL = "robertnquiroz@gmail.com"  # Cambia esto por tu correo
SEARCH_TERM = (
    '"respiratory syncytial virus"[All Fields] AND '
    '"Expression profiling by high throughput sequencing"[DataSet Type] AND '
    'GSE[Entry Type]'
)
OUTPUT_FILE = "todos_estudios_VSR_GEO.csv"
BATCH_SIZE = 100
MAX_RETRIES = 3
RETRY_DELAY = 5  # Segundos

# Configurar logging para mejor registro de la ejecución
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

def fetch_with_retry(func, *args, **kwargs) -> Any:
    """Ejecuta una función de Entrez con reintentos en caso de fallo de conexión."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return func(*args, **kwargs)
        except (HTTPError, Exception) as e:
            if attempt == MAX_RETRIES:
                logger.error(f"Fallo definitivo tras {MAX_RETRIES} intentos: {e}")
                raise e
            logger.warning(f"Error de conexión (Intento {attempt}/{MAX_RETRIES}): {e}. Reintentando en {RETRY_DELAY}s...")
            time.sleep(RETRY_DELAY)

def main():
    Entrez.email = ENTREZ_EMAIL
    
    logger.info("Iniciando búsqueda de estudios para RSV en GEO...")

    # 2. Primera búsqueda: Solo para saber el TOTAL de resultados (retmax=1)
    try:
        handle = fetch_with_retry(Entrez.esearch, db="gds", term=SEARCH_TERM, retmax=1)
        search_results = Entrez.read(handle)
        handle.close()
    except Exception as e:
        logger.error("No se pudo obtener el total de estudios. Terminando ejecución.")
        return

    total_estudios = int(search_results.get("Count", 0))
    if total_estudios == 0:
        logger.info("No se encontraron estudios con esos criterios.")
        return
        
    logger.info(f"¡Se encontraron {total_estudios} estudios en total!")

    # 3. Segunda búsqueda: Pedimos TODOS los IDs
    logger.info("Obteniendo lista de IDs de todos los estudios...")
    try:
        handle = fetch_with_retry(Entrez.esearch, db="gds", term=SEARCH_TERM, retmax=total_estudios)
        search_results = Entrez.read(handle)
        handle.close()
    except Exception as e:
        logger.error("No se pudieron obtener los IDs de los estudios.")
        return

    id_list = search_results.get("IdList", [])

    # 4. Obtener detalles en LOTES para no saturar el servidor
    estudios: List[Dict[str, Any]] = []
    
    logger.info(f"Descargando metadatos en lotes de {BATCH_SIZE}...")

    for start in range(0, total_estudios, BATCH_SIZE):
        end = min(total_estudios, start + BATCH_SIZE)
        logger.info(f"Procesando estudios del {start + 1} al {end} de {total_estudios}...")
        
        batch_ids = id_list[start:end]
        
        try:
            handle_summary = fetch_with_retry(
                Entrez.esummary, 
                db="gds", 
                id=",".join(batch_ids)
            )
            summaries = Entrez.read(handle_summary)
            handle_summary.close()
            
            # Extraemos la información con manejo seguro usando .get()
            for summary in summaries:
                accession = "GSE" + summary["GSE"] if "GSE" in summary else summary.get("Accession", "N/A")
                estudios.append({
                    "GEO_ID": accession,
                    "Titulo": summary.get("title", "Sin título"),
                    "Muestras": summary.get("n_samples", "N/A"),
                    "Fecha_Publicacion": summary.get("PDAT", "N/A"),
                    "Taxonomia": summary.get("taxon", "N/A")
                })
        except Exception as e:
            logger.error(f"Fallo al procesar el lote del {start + 1} al {end}: {e}. Saltando al siguiente lote.")
            
        # Pausa de 1 segundo para evitar bloqueos del servidor (rate limiting)
        time.sleep(1)

    # 5. Convertir a tabla y guardar en CSV
    if not estudios:
        logger.error("No se extrajo información de ningún estudio.")
        return

    df_todos_los_estudios = pd.DataFrame(estudios)
    
    # Exportamos todo a un archivo CSV
    df_todos_los_estudios.to_csv(OUTPUT_FILE, index=False)
    
    logger.info(f"¡Descarga completada! Datos guardados en: {OUTPUT_FILE}")
    logger.info("--- Primeros 5 estudios encontrados ---")
    
    # Imprimir los primeros 5 de forma tabulada
    print(df_todos_los_estudios.head().to_string(index=False))

if __name__ == "__main__":
    main()
