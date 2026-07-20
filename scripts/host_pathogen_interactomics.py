import os
import io
import requests
import urllib.request
import zipfile
import pandas as pd
import numpy as np
import networkx as nx
# pyrefly: ignore [missing-import]
import decoupler as dc
# pyrefly: ignore [missing-import]
import py4cytoscape as p4c
from Bio import Entrez, SeqIO
import warnings

warnings.filterwarnings('ignore')

# =====================================================================
# FASE 1: BÚSQUEDA DE DATOS PROTEÓMICOS REALES (PRIDE Archive API)
# =====================================================================
def buscar_datasets_pride(termino_busqueda="Respiratory Syncytial Virus"):
    """
    Consulta la API REST de PRIDE para encontrar proyectos públicos de proteómica 
    relacionados con el VSR.
    """
    print(f"--- FASE 1: Buscando datasets proteómicos en PRIDE para '{termino_busqueda}' ---")
    url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword={termino_busqueda}"
    
    try:
        respuesta = requests.get(url)
        respuesta.raise_for_status()
        datos = respuesta.json()
        
        # Extraer los proyectos encontrados
        proyectos = []
        for doc in datos[:5]: # Top 5
            proyectos.append({
                'Accession': doc.get('accession'),
                'Titulo': doc.get('title'),
                'Especies': ", ".join(doc.get('organisms', [])),
                'Instrumento': ", ".join(doc.get('instruments', []))
            })
            
        df_pride = pd.DataFrame(proyectos)
        print("\nTop proyectos encontrados listos para descarga:")
        print(df_pride.to_string(index=False))
        return df_pride
    except Exception as e:
        print(f"Error conectando a PRIDE: {e}")
        return None

def descargar_datos_vsr_pride():
    """
    Busca proyectos de proteómica del VSR en PRIDE y descarga los resultados 
    del proyecto más relevante (proteinGroups.txt) a la carpeta local 'datos/'.
    """
    termino = "Respiratory Syncytial Virus"
    print(f"\n--- FASE 1.5: Descargando datos experimentales reales de PRIDE ---")
    url_busqueda = f"https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword={termino}"
    
    try:
        respuesta = requests.get(url_busqueda)
        respuesta.raise_for_status()
        datos_busqueda = respuesta.json()
        
        proyectos = datos_busqueda.get('_embedded', {}).get('compactprojects', [])
        if not proyectos:
            proyectos = datos_busqueda if isinstance(datos_busqueda, list) else []
            
        if not proyectos:
            print("❌ No se encontraron proyectos.")
            return False
            
        proyecto_top = proyectos[0]
        accession_pxd = proyecto_top.get('accession')
        print(f"✅ Proyecto Seleccionado: {accession_pxd} - {proyecto_top.get('title')}")
        
        url_archivos = f"https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={accession_pxd}"
        resp_archivos = requests.get(url_archivos)
        resp_archivos.raise_for_status()
        archivos = resp_archivos.json()
        
        ftp_link, nombre_archivo = None, None
        
        for archivo in archivos:
            categoria = archivo.get('fileCategory', {}).get('value', '')
            nombre = archivo.get('fileName', '')
            if 'RESULT' in categoria or nombre.endswith('.zip') or 'txt' in nombre:
                locaciones = archivo.get('publicFileLocations', [])
                if locaciones:
                    ftp_link = locaciones[0].get('value')
                    nombre_archivo = nombre
                    break
                    
        if ftp_link:
            directorio_destino = "datos"
            os.makedirs(directorio_destino, exist_ok=True)
            ruta_descarga = os.path.join(directorio_destino, nombre_archivo)
            
            print(f"Descargando matriz de resultados: {nombre_archivo}...")
            urllib.request.urlretrieve(ftp_link, ruta_descarga)
            
            if ruta_descarga.endswith('.zip'):
                print("Descomprimiendo archivo...")
                with zipfile.ZipFile(ruta_descarga, 'r') as zip_ref:
                    for file_info in zip_ref.infolist():
                        if 'proteinGroups.txt' in file_info.filename:
                            file_info.filename = "proteinGroups.txt"
                            zip_ref.extract(file_info, directorio_destino)
                            print("✅ Matriz proteinGroups.txt extraída con éxito en /datos.")
                            return True
            else:
                print(f"✅ Archivo descargado en: {ruta_descarga}")
                return True
        else:
            print(f"❌ El proyecto {accession_pxd} no tiene archivos txt públicos accesibles.")
            return False
            
    except Exception as e:
        print(f"Error en la descarga: {e}")
        return False

# =====================================================================
# FASE 2: SECUENCIAS Y ENSAMBLAJE FASTA (Biopython)
# =====================================================================
def descargar_genoma_vsr(email):
    """
    Descarga el proteoma de referencia del VSR desde NCBI.
    Indispensable para generar el archivo FASTA de búsqueda para motores como MaxQuant.
    """
    print("\n--- FASE 2: Descarga de secuencias de referencia (NCBI) ---")
    Entrez.email = email
    id_referencia = "NC_038235"
    
    try:
        handle = Entrez.efetch(db="nucleotide", id=id_referencia, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        
        proteinas = []
        for feature in record.features:
            if feature.type == "CDS":
                gen = feature.qualifiers.get('gene', ['Unknown'])[0]
                secuencia = feature.qualifiers.get('translation', [''])[0]
                proteinas.append({'Viral_Gene': gen, 'Sequence': secuencia})
                
        df_virus = pd.DataFrame(proteinas)
        print(f"Extraídas {len(df_virus)} proteínas virales del genoma {id_referencia}.")
        return df_virus
    except Exception as e:
        print(f"Error en NCBI: {e}")
        return pd.DataFrame()

# =====================================================================
# FASE 3: PROCESAMIENTO DE MATRIZ CUANTITATIVA (Pandas)
# =====================================================================
def procesar_matriz_maxquant_real(ruta_archivo, columnas_mock, columnas_vsr):
    """
    Carga un archivo proteinGroups.txt real de MaxQuant, elimina contaminantes, 
    resuelve nombres de genes y calcula el Log2 Fold Change.
    
    Parámetros:
    - ruta_archivo: str, ruta local al archivo (ej. 'datos/proteinGroups.txt')
    - columnas_mock: list, nombres exactos de las columnas de control (células sanas)
    - columnas_vsr: list, nombres exactos de las columnas de infección (VSR)
    """
    print(f"\n--- FASE 3: Procesamiento Real de {ruta_archivo} ---")
    
    try:
        # 1. CARGA DE DATOS
        # Se usa low_memory=False por la gran cantidad de columnas de metadatos en MaxQuant
        df = pd.read_csv(ruta_archivo, sep='\t', low_memory=False)
        filas_originales = len(df)
        
        # 2. CONTROL DE CALIDAD (Filtro de contaminantes y decoys)
        # MaxQuant marca estas anomalías con un '+'
        filtros_qc = ['Reverse', 'Potential contaminant', 'Only identified by site']
        for col in filtros_qc:
            if col in df.columns:
                df = df[df[col] != '+']
                
        print(f"Control de Calidad: Retenidas {len(df)} proteínas (de {filas_originales} iniciales).")
        
        # 3. RESOLUCIÓN DE NOMBRES DE GENES
        # Eliminar filas que no mapean a ningún gen
        df = df.dropna(subset=['Genes'])
        # Conservamos el gen primario.
        df['Gene'] = df['Genes'].astype(str).str.split(';').str[0]
        
        # 4. TRANSFORMACIÓN DE INTENSIDADES (LFQ o iBAQ)
        # Reemplazar los ceros absolutos por NaN. En espectrometría, un '0' significa 
        # que el péptido no superó el límite de detección, no que la concentración biológica sea cero.
        todas_las_columnas = columnas_mock + columnas_vsr
        for col in todas_las_columnas:
            if col not in df.columns:
                raise ValueError(f"Error: La columna '{col}' no existe en el archivo proteinGroups.txt")
                
        df[todas_las_columnas] = df[todas_las_columnas].replace(0, np.nan)
        
        # 5. CÁLCULO DE LOG2 FOLD CHANGE
        # Aplicamos logaritmo base 2 para normalizar la distribución de intensidades
        df_log = np.log2(df[todas_las_columnas])
        
        # Calculamos la media de las réplicas biológicas/técnicas para cada condición
        df['Media_Mock'] = df_log[columnas_mock].mean(axis=1)
        df['Media_VSR'] = df_log[columnas_vsr].mean(axis=1)
        
        # Fold Change = Condición (VSR) - Control (Mock)
        df['Log2FC'] = df['Media_VSR'] - df['Media_Mock']
        
        # 6. AGREGACIÓN FINAL
        # Eliminar proteínas que no tienen un FC válido (ej. solo se detectaron en una condición)
        df = df.dropna(subset=['Log2FC'])
        
        # Si múltiples proteínas (isoformas) mapean al mismo genotipo, promediamos su comportamiento
        df_final = df.groupby('Gene')['Log2FC'].mean().reset_index()
        
        # Establecer el gen como índice para que Decoupler y NetworkX lo interpreten correctamente
        df_final = df_final.set_index('Gene')
        
        print(f"Matriz lista: {len(df_final)} genes únicos listos para inferencia de redes.")
        return df_final
        
    except FileNotFoundError:
        print(f"❌ Error: No se encontró el archivo en {ruta_archivo}")
        return pd.DataFrame()
    except BaseException as e:
        print(f"❌ Error procesando la matriz: {e}")
        return pd.DataFrame()

# =====================================================================
# FASE 4: INFERENCIA DE VÍAS DE SEÑALIZACIÓN (Decoupler Nativo)
# =====================================================================
def inferir_actividad_tfs(df_proteomica):
    """
    Utiliza el Univariate Linear Model (ULM) para inferir factores de transcripción.
    Extrae las redes CollecTRI y DoRothEA de forma nativa a través de decoupleR
    (evadiendo la API REST externa), creando una red Transcriptional unificada.
    """
    print("\n--- FASE 4: Inferencia de Factores de Transcripción (ULM) ---")
    ruta_local = "red_transcriptional_local.csv"
    
    # 1. OBTENCIÓN DE LA RED REGULATORIA (CollecTRI + DoRothEA)
    if os.path.exists(ruta_local):
        print("Cargando red Transcriptional (CollecTRI + DoRothEA) desde el disco local...")
        net_conocimiento = pd.read_csv(ruta_local)
    else:
        print("Ensamblando red Transcriptional nativa vía decoupleR (evitando API externa)...")
        try:
            print("  - Cargando CollecTRI...")
            net_col = dc.op.collectri(organism='human')
            print("  - Cargando DoRothEA...")
            net_dor = dc.op.dorothea(organism='human')
            
            # Unir ambas redes para tener el set 'transcriptional' completo
            net_conocimiento = pd.concat([
                net_col[['source', 'target', 'weight']], 
                net_dor[['source', 'target', 'weight']]
            ])
            # Eliminar redundancias priorizando la primera aparición
            net_conocimiento = net_conocimiento.drop_duplicates(subset=['source', 'target'])
            
            # Guardar en caché local
            net_conocimiento.to_csv(ruta_local, index=False)
            print(f"✅ Red Transcriptional ensamblada con éxito ({len(net_conocimiento)} interacciones).")
            
        except Exception as e:
            print(f"\n❌ Fallo crítico al ensamblar la red localmente con decoupleR: {e}")
            return None, None
            
    # 2. INFERENCIA ESTADÍSTICA
    print("Iniciando cálculo del modelo lineal univariado (ULM)...")
    try:
        # Para versiones antiguas de decoupleR (1.x)
        estimacion, p_valores = dc.run_ulm(
            mat=df_proteomica.T, 
            net=net_conocimiento,
            source='source', target='target', weight='weight',
            verbose=False
        )
    except AttributeError:
        # Para versiones modernas de decoupleR (>= 2.x)
        estimacion, p_valores = dc.mt.ulm(
            df_proteomica.T, 
            net_conocimiento,
            verbose=False
        )
    
    # 3. EXTRACCIÓN DE RESULTADOS
    tfs_activos = estimacion.T
    tfs_activos.columns = ['Score_ULM']
    tfs_activos = tfs_activos.sort_values(by='Score_ULM', ascending=False)
    
    # Extraer los Top 10 TFs para tener una red enorme y representativa
    top_n = 10
    tfs_top_lista = tfs_activos.head(top_n).index.tolist()
    
    print(f"\nTop {top_n} Factores de Transcripción Inferidos:")
    print(tfs_activos.head(top_n))
    print(f"\nFactores de transcripción dominantes que pasarán a la Fase 5: {tfs_top_lista}")
    
    return tfs_top_lista, net_conocimiento

# =====================================================================
# FASE 5: INTEGRACIÓN TOPOLÓGICA (NetworkX)
# =====================================================================
def construir_red_interactoma(tfs_top_lista, net_collectri, df_proteomica):
    """
    Construye un grafo dirigido combinando la red regulatoria de los TFs principales 
    con las interacciones virales conocidas.
    """
    print(f"\n--- FASE 5: Ensamblando interactoma para los Top {len(tfs_top_lista)} TFs ---")
    
    # 1. Extraer subred biológica de los TFs desde la red de conocimiento
    subred_df = net_collectri[net_collectri['source'].isin(tfs_top_lista)].copy()
    G = nx.from_pandas_edgelist(
        subred_df, source='source', target='target', edge_attr='weight', create_using=nx.DiGraph()
    )
    
    # 2. Agregar el Log2FC experimental a los nodos
    dict_log2fc = df_proteomica['Log2FC'].to_dict()
    nx.set_node_attributes(G, dict_log2fc, 'Log2FC')
    
    # Asegurar que los TFs centrales y otros nodos sin cuantificación no tengan NaNs (vital para Cytoscape)
    for nodo in G.nodes():
        if 'Log2FC' not in G.nodes[nodo] or pd.isna(G.nodes[nodo]['Log2FC']):
            G.nodes[nodo]['Log2FC'] = 0.0
            
    # 3. Agregar interacciones virales directas (VSR_F y VSR_G)
    # Por ejemplo, la proteína G del VSR interactúa con CX3CR1, que activa cascadas inmunes
    interacciones_virales = [
        ("VSR_G", "CX3CR1", {"weight": 1.0, "tipo": "adhesion_viral"}),
        ("VSR_F", "TLR4", {"weight": 1.0, "tipo": "entrada_viral"})
    ]
    # Conectar el virus a los TFs inferidos
    for tf in tfs_top_lista:
        interacciones_virales.append(("TLR4", tf, {"weight": 1.0, "tipo": "cascada_inmune"}))
        
    G.add_edges_from([(u, v, attr) for u, v, attr in interacciones_virales])
    
    # Asignar Log2FC de 0.0 a las proteínas virales para la visualización
    G.nodes["VSR_G"]['Log2FC'] = 0.0
    G.nodes["VSR_F"]['Log2FC'] = 0.0
    G.nodes["TLR4"]['Log2FC'] = 0.0
    G.nodes["CX3CR1"]['Log2FC'] = 0.0
    
    print("Calculando métricas topológicas (Betweenness Centrality)...")
    dict_betweenness = nx.betweenness_centrality(G)
    nx.set_node_attributes(G, dict_betweenness, 'Betweenness')
    
    print(f"Red ensamblada: {G.number_of_nodes()} nodos y {G.number_of_edges()} aristas.")
    return G
# =====================================================================
# FASE 6: RENDERIZADO Y EXPORTACIÓN PARA PUBLICACIÓN (py4cytoscape)
# =====================================================================
def publicar_en_cytoscape(G, tfs_top_lista):
    """
    Transmite la red a Cytoscape, aplica gradientes continuos basados en Log2FC,
    mapea la Centralidad (Betweenness) al tamaño del nodo, y exporta la figura.
    """
    print("\n--- FASE 6: Renderizado en Cytoscape y Exportación ---")
    try:
        p4c.cytoscape_ping()
    except Exception:
        print("❌ Cytoscape no está ejecutándose. Abre Cytoscape e intenta de nuevo.")
        return
        
    estilo = "Estilo_VSR_Q1"
    
    # Transmitir red
    titulo_red = f"Interactoma VSR - Top {len(tfs_top_lista)} TFs"
    p4c.create_network_from_networkx(G, title=titulo_red, collection="VSR_Omics")
    
    # Configurar Estilo Base
    if estilo in p4c.get_visual_style_names():
        p4c.delete_visual_style(estilo)
    p4c.create_visual_style(estilo)
    
    p4c.set_node_shape_default('ELLIPSE', style_name=estilo)
    p4c.set_node_label_mapping('name', style_name=estilo)
    p4c.set_visual_property_default({'visualProperty': 'NETWORK_BACKGROUND_PAINT', 'value': '#FFFFFF'}, style_name=estilo)
    
    # Mapeo Continuo de Proteómica (Azul = Inhibido, Rojo = Sobreexpresado)
    p4c.set_node_color_mapping(
        table_column='Log2FC', 
        table_column_values=[-2.0, 0.0, 2.0], 
        colors=['#3182BD', '#DDDDDD', '#DE2D26'], 
        mapping_type='c', style_name=estilo
    )
    
    # Mapeo Continuo de Tamaño basado en la Centralidad (Cuellos de Botella)
    betweenness_vals = [d.get('Betweenness', 0) for n, d in G.nodes(data=True)]
    max_b = max(betweenness_vals) if betweenness_vals else 0.1
    p4c.set_node_size_mapping(
        table_column='Betweenness',
        table_column_values=[0.0, max_b],
        sizes=[25, 120],
        mapping_type='c', style_name=estilo
    )
    
    # Mapeo Continuo de Aristas (Naranja = Inhibición, Verde = Activación)
    p4c.set_edge_color_mapping(
        table_column='weight', table_column_values=[-1.0, 0.0, 1.0], 
        colors=['#E6550D', '#DDDDDD', '#31A354'], mapping_type='c', style_name=estilo
    )
    
    p4c.set_visual_style(estilo)
    p4c.layout_network('force-directed')
    
    # Exportación
    ruta_salida = os.path.join(os.getcwd(), "Figura_Interactoma_VSR_300DPI.png")
    p4c.export_image(filename=ruta_salida, type='PNG', resolution=300, zoom=100.0)
    print(f"✅ Figura de publicación exportada a: {ruta_salida}")

# =====================================================================
# EJECUCIÓN DEL PIPELINE
# =====================================================================
if __name__ == "__main__":
    # 1. Explorar PRIDE (Localizar datos)
    df_proyectos = buscar_datasets_pride()
    
    # 1.5. Descargar los datos experimentales reales
    # datos_descargados = descargar_datos_vsr_pride()
    datos_descargados = True # Usando dataset DIA-NN manual (PXD039212)
    
    # 2. Obtener secuencias 
    df_secuencias = descargar_genoma_vsr("investigacion@tdea.edu.co")
    
    if datos_descargados:
        # 3. Procesar datos reales del MS (Log2FC)
        # ATENCIÓN: Estos nombres deben coincidir con las columnas del archivo descargado
        columnas_control_reales = [
            r'E:\Data\Morgan_W_Mann\hSAEC_RSV_ZL0454_GBUP\DIA_hSAEC_RSV_ZL54_NoDrugCon_S1_200ng_3_Slot2-34_1_3954.d',
            r'E:\Data\Morgan_W_Mann\hSAEC_RSV_ZL0454_GBUP\DIA_hSAEC_RSV_ZL54_NoDrugCon_S2_200ng_3_Slot2-35_1_3956.d',
            r'E:\Data\Morgan_W_Mann\hSAEC_RSV_ZL0454_GBUP\DIA_hSAEC_RSV_ZL54_NoDrugCon_S3_200ng_2_Slot2-36_1_3957.d'
        ]
        columnas_infeccion_reales = [
            r'E:\Data\Morgan_W_Mann\hSAEC_RSV_ZL0454_GBUP\DIA_hSAEC_RSV_ZL54_NoDrugVirus_S1_200ng_1_Slot2-40_1_3946.d',
            r'E:\Data\Morgan_W_Mann\hSAEC_RSV_ZL0454_GBUP\DIA_hSAEC_RSV_ZL54_NoDrugVirus_S2_200ng_1_Slot2-41_1_3947.d',
            r'E:\Data\Morgan_W_Mann\hSAEC_RSV_ZL0454_GBUP\DIA_hSAEC_RSV_ZL54_NoDrugVirus_S3_200ng_1_Slot2-42_1_3948.d'
        ]
        
        df_cuantificacion = procesar_matriz_maxquant_real(
            ruta_archivo="data/proteomics/DIA-NN Output/report.pg_matrix.tsv", 
            columnas_mock=columnas_control_reales, 
            columnas_vsr=columnas_infeccion_reales
        )
        
        # 4. Inferir TFs
        if not df_cuantificacion.empty:
            tfs_top_lista, red_conocimiento = inferir_actividad_tfs(df_cuantificacion)
            
            if tfs_top_lista is not None and red_conocimiento is not None:
                # 5. Modelado Matemático
                grafo = construir_red_interactoma(tfs_top_lista, red_conocimiento, df_cuantificacion)
                
                # 6. Visualización (Requiere Cytoscape abierto)
                publicar_en_cytoscape(grafo, tfs_top_lista)
            else:
                print("❌ No se construyó el interactoma porque no se pudo inferir la actividad de los TFs.")
        else:
            print("❌ No se construyó el interactoma porque la matriz quedó vacía tras el filtrado.")
    else:
        print("❌ Pipeline detenido: No se pudo descargar la matriz experimental de PRIDE.")
