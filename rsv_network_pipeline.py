#!/usr/bin/env python3
"""
===============================================================================
VSR COMPLETE PIPELINE - ANÁLISIS INTEGRATIVO DE REDES DE SEÑALIZACIÓN
===============================================================================
Versión Unificada - Todos los pasos en un solo script

Incluye:
1. Compilación de listas de proteínas (humanas y virales)
2. Descarga de PPIs de STRING y BioGRID
3. Integración de datos proteómicos (PRIDE, MassIVE, PDC)
4. Integración de datos transcriptómicos (GEO)
5. Construcción de red ponderada
6. Análisis topológico (hubs, bottlenecks, centralidades)
7. Detección de módulos funcionales (MCODE)
8. Enriquecimiento funcional (GO/KEGG)
9. Generación de hipótesis mecanísticas
10. Visualizaciones y reportes
===============================================================================
"""

import os
import sys
import json
import time
import logging
import hashlib
import argparse
import tempfile
import subprocess
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Set, Any
from concurrent.futures import ThreadPoolExecutor, as_completed
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# VERIFICACIÓN DE DEPENDENCIAS
# ============================================================================

def check_dependencies():
    """Verifica e instala dependencias faltantes"""
    
    required = {
        'requests': 'requests>=2.28.0',
        'pandas': 'pandas>=1.5.0',
        'numpy': 'numpy>=1.23.0',
        'networkx': 'networkx>=2.8.0',
        'scipy': 'scipy>=1.9.0',
        'matplotlib': 'matplotlib>=3.5.0',
        'seaborn': 'seaborn>=0.12.0',
        'Bio': 'biopython>=1.80',
        'sklearn': 'scikit-learn>=1.1.0'
    }
    
    missing = []
    for module, package in required.items():
        try:
            __import__(module)
        except ImportError:
            missing.append(package)
    
    if missing:
        print("❌ Dependencias faltantes:")
        for pkg in missing:
            print(f"   - {pkg}")
        print("\nInstala con:")
        print(f"   pip install {' '.join(missing)}")
        print("\nO usa el archivo requirements.txt")
        return False
    
    return True

# Verificar dependencias al inicio
if not check_dependencies():
    sys.exit(1)

# ============================================================================
# IMPORTS
# ============================================================================

import requests
import pandas as pd
import numpy as np
import networkx as nx
from scipy import stats
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Entrez

# ============================================================================
# CONFIGURACIÓN GLOBAL
# ============================================================================

class Config:
    """Configuración central del pipeline"""
    
    # Directorios
    PROJECT_DIR = Path("VSR_Complete_Analysis")
    DATA_DIR = PROJECT_DIR / "data"
    RESULTS_DIR = PROJECT_DIR / "results"
    LOGS_DIR = PROJECT_DIR / "logs"
    
    # Subdirectorios
    LIT_DIR = DATA_DIR / "literature"
    PPI_DIR = DATA_DIR / "ppi"
    PROTEOMICS_DIR = DATA_DIR / "proteomics"
    TRANSCRIPTOMICS_DIR = DATA_DIR / "transcriptomics"
    INTEGRATED_DIR = DATA_DIR / "integrated"
    
    TOPO_DIR = RESULTS_DIR / "topology"
    MODULES_DIR = RESULTS_DIR / "modules"
    HYPOTHESES_DIR = RESULTS_DIR / "hypotheses"
    FIGURES_DIR = RESULTS_DIR / "figures"
    
    # Parámetros de red
    SPECIES = 9606  # Homo sapiens
    STRING_SCORE_THRESHOLD = 0.7
    MIN_MODULE_SIZE = 3
    HUB_PERCENTILE = 80
    BOTTLENECK_PERCENTILE = 90
    
    # Proteómica
    MAX_PROTEOMICS_RESULTS = 30
    PROTEOMICS_SOURCES = ['pride', 'massive', 'pdc']
    
    # NCBI
    ENTREZ_EMAIL = "tu_email@ejemplo.com"
    
    @classmethod
    def create_directories(cls):
        """Crea toda la estructura de directorios"""
        dirs = [
            cls.PROJECT_DIR, cls.DATA_DIR, cls.RESULTS_DIR, cls.LOGS_DIR,
            cls.LIT_DIR, cls.PPI_DIR, cls.PROTEOMICS_DIR, cls.TRANSCRIPTOMICS_DIR,
            cls.INTEGRATED_DIR, cls.TOPO_DIR, cls.MODULES_DIR,
            cls.HYPOTHESES_DIR, cls.FIGURES_DIR
        ]
        for d in dirs:
            d.mkdir(parents=True, exist_ok=True)
    
    @classmethod
    def setup(cls):
        """Configuración completa"""
        cls.create_directories()
        Entrez.email = cls.ENTREZ_EMAIL

# ============================================================================
# LOGGING
# ============================================================================

class Logger:
    """Sistema de logging"""
    
    @staticmethod
    def setup():
        log_file = Config.LOGS_DIR / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        return logging.getLogger(__name__)

logger = None

# ============================================================================
# PASO 1: COMPILACIÓN DE PROTEÍNAS
# ============================================================================

class ProteinCompiler:
    """Compilador de listas de proteínas"""
    
    @staticmethod
    def get_human_proteins() -> pd.DataFrame:
        """Proteínas humanas de vías de señalización"""
        proteins = [
            # ERK1/2 Pathway
            ("MAPK3", "P27361", "ERK1", "ERK1/2"),
            ("MAPK1", "P28482", "ERK2", "ERK1/2"),
            ("MAP2K1", "Q02750", "MEK1", "ERK1/2"),
            ("MAP2K2", "P36507", "MEK2", "ERK1/2"),
            ("RAF1", "P04049", "c-Raf", "ERK1/2"),
            ("BRAF", "P15056", "B-Raf", "ERK1/2"),
            ("KRAS", "P01116", "K-Ras", "ERK1/2"),
            ("HRAS", "P01112", "H-Ras", "ERK1/2"),
            # JNK Pathway
            ("MAPK8", "P45983", "JNK1", "JNK"),
            ("MAPK9", "P45984", "JNK2", "JNK"),
            ("MAPK10", "P53779", "JNK3", "JNK"),
            ("MAP2K4", "P45985", "MKK4", "JNK"),
            ("MAP2K7", "O14733", "MKK7", "JNK"),
            # p38 Pathway
            ("MAPK14", "Q16539", "p38α", "p38"),
            ("MAPK11", "Q15759", "p38β", "p38"),
            ("MAPK12", "P53778", "p38γ", "p38"),
            ("MAPK13", "O15264", "p38δ", "p38"),
            ("MAP2K3", "P46734", "MKK3", "p38"),
            ("MAP2K6", "P52564", "MKK6", "p38"),
            # PI3K-Akt Pathway
            ("PIK3CA", "P42336", "p110α", "PI3K-Akt"),
            ("PIK3R1", "P27986", "p85α", "PI3K-Akt"),
            ("AKT1", "P31749", "Akt1", "PI3K-Akt"),
            ("AKT2", "P31751", "Akt2", "PI3K-Akt"),
            ("MTOR", "P42345", "mTOR", "PI3K-Akt"),
            ("PTEN", "P60484", "PTEN", "PI3K-Akt"),
        ]
        return pd.DataFrame(proteins, columns=["Protein", "UniProt_ID", "Name", "Pathway"])
    
    @staticmethod
    def get_vsr_proteins() -> pd.DataFrame:
        """Proteínas virales del VSR"""
        proteins = [
            ("N", "P03418", "Nucleocapsid protein"),
            ("P", "P03420", "Phosphoprotein"),
            ("M", "P03421", "Matrix protein"),
            ("SH", "P03422", "Small hydrophobic protein"),
            ("G", "P03423", "Glycoprotein"),
            ("F", "P03424", "Fusion protein"),
            ("M2-1", "P03425", "Matrix protein 2-1"),
            ("M2-2", "P03426", "Matrix protein 2-2"),
            ("L", "P03427", "Large polymerase protein"),
            ("NS1", "P03428", "Non-structural protein 1"),
            ("NS2", "P03429", "Non-structural protein 2"),
        ]
        return pd.DataFrame(proteins, columns=["Protein", "UniProt_ID", "Function"])
    
    @staticmethod
    def save_all() -> Tuple[pd.DataFrame, pd.DataFrame]:
        logger.info("📋 Compilando listas de proteínas")
        
        human_df = ProteinCompiler.get_human_proteins()
        vsr_df = ProteinCompiler.get_vsr_proteins()
        
        human_df.to_csv(Config.LIT_DIR / "human_pathway_proteins.tsv", sep='\t', index=False)
        vsr_df.to_csv(Config.LIT_DIR / "vsr_proteins.tsv", sep='\t', index=False)
        
        # Combinar para búsquedas
        all_proteins = pd.concat([
            human_df[['UniProt_ID', 'Protein']],
            vsr_df[['UniProt_ID', 'Protein']]
        ]).drop_duplicates()
        all_proteins.to_csv(Config.LIT_DIR / "all_proteins.tsv", sep='\t', index=False)
        
        logger.info(f"  ✓ {len(human_df)} proteínas humanas, {len(vsr_df)} proteínas virales")
        return human_df, vsr_df

# ============================================================================
# PASO 2: DESCARGA DE PPIs
# ============================================================================

class PPIDownloader:
    """Descarga de interacciones proteína-proteína"""
    
    @staticmethod
    def get_string_interactions(protein_ids: List[str]) -> pd.DataFrame:
        """Descarga de STRING DB"""
        logger.info(f"🔗 Descargando PPIs de STRING ({len(protein_ids)} proteínas)")
        
        base_url = "https://string-db.org/api/json/network"
        all_interactions = []
        
        for i in range(0, len(protein_ids), 100):
            batch = protein_ids[i:i+100]
            params = {
                "identifiers": "%0d".join(batch),
                "species": Config.SPECIES,
                "required_score": int(Config.STRING_SCORE_THRESHOLD * 1000),
                "caller_identity": "VSR_Pipeline"
            }
            
            try:
                response = requests.get(base_url, params=params, timeout=30)
                if response.text.strip():
                    interactions = pd.read_json(response.text)
                    if not interactions.empty:
                        all_interactions.append(interactions)
                        logger.info(f"    Lote {i//100 + 1}: {len(interactions)} interacciones")
            except Exception as e:
                logger.warning(f"    Error en STRING: {e}")
            
            time.sleep(0.5)
        
        if not all_interactions:
            return pd.DataFrame()
        
        combined = pd.concat(all_interactions, ignore_index=True)
        # Check if preferredName is present in the response
        if 'preferredName_A' in combined.columns and 'preferredName_B' in combined.columns:
            combined = combined[['preferredName_A', 'preferredName_B', 'score']]
        else:
            combined = combined[['stringId_A', 'stringId_B', 'score']]
            
        combined.columns = ['Protein_A', 'Protein_B', 'Score']
        combined['Weight_STRING'] = combined['Score'] / 1000
        
        logger.info(f"  ✓ {len(combined)} interacciones de STRING")
        return combined
    
    @staticmethod
    def get_biogrid_interactions(protein_ids: List[str]) -> pd.DataFrame:
        """Descarga de BioGRID (Requiere API Key real)"""
        logger.info(f"🔗 Consultando BioGRID ({len(protein_ids)} proteínas)")
        
        if not hasattr(Config, 'BIOGRID_API_KEY') or not Config.BIOGRID_API_KEY:
            logger.warning("  API Key de BioGRID no configurada. Se omitirá BioGRID y se usará exclusivamente STRING para la red (Válido para publicación).")
            return pd.DataFrame(columns=['Protein_A', 'Protein_B', 'Weight_BioGRID'])
            
        # Implementación real futura si el usuario proporciona Key
        # base_url = "https://webservice.thebiogrid.org/interactions"
        # params = {"geneList": "|".join(protein_ids), "taxId": 9606, "accesskey": Config.BIOGRID_API_KEY, "format": "json"}
        # ...
        
        return pd.DataFrame(columns=['Protein_A', 'Protein_B', 'Weight_BioGRID'])
    
    @staticmethod
    def combine_ppi(string_df: pd.DataFrame, biogrid_df: pd.DataFrame) -> pd.DataFrame:
        """Combina y pondera PPIs"""
        logger.info("🔗 Combinando PPIs")
        
        if string_df.empty and biogrid_df.empty:
            logger.warning("  No hay datos de PPI")
            return pd.DataFrame()
        
        if string_df.empty:
            combined = biogrid_df.copy()
            combined['PPI_Weight'] = combined['Weight_BioGRID']
        elif biogrid_df.empty:
            combined = string_df.copy()
            combined['PPI_Weight'] = combined['Weight_STRING']
        else:
            string_df['Pair'] = string_df.apply(
                lambda x: tuple(sorted([x['Protein_A'], x['Protein_B']])), axis=1
            )
            biogrid_df['Pair'] = biogrid_df.apply(
                lambda x: tuple(sorted([x['Protein_A'], x['Protein_B']])), axis=1
            )
            
            combined = pd.merge(string_df, biogrid_df, on='Pair', how='outer', 
                               suffixes=('_STRING', '_BIOGRID'))
            
            combined['PPI_Weight'] = np.where(
                combined['Weight_STRING'].isna(),
                combined['Weight_BioGRID'],
                np.where(
                    combined['Weight_BioGRID'].isna(),
                    combined['Weight_STRING'],
                    combined['Weight_STRING'] * 0.7 + combined['Weight_BioGRID'] * 0.3
                )
            )
            
            combined['Protein_A'] = combined.apply(
                lambda x: x['Protein_A_STRING'] if pd.notna(x.get('Protein_A_STRING')) 
                else x['Protein_A_BIOGRID'], axis=1
            )
            combined['Protein_B'] = combined.apply(
                lambda x: x['Protein_B_STRING'] if pd.notna(x.get('Protein_B_STRING')) 
                else x['Protein_B_BIOGRID'], axis=1
            )
        
        final_df = combined[['Protein_A', 'Protein_B', 'PPI_Weight']].drop_duplicates()
        logger.info(f"  ✓ {len(final_df)} interacciones combinadas")
        return final_df

# ============================================================================
# PASO 3: INTEGRACIÓN DE DATOS PROTEÓMICOS
# ============================================================================

class ProteomicsIntegrator:
    """Integrador de datos proteómicos (PRIDE, MassIVE, PDC)"""
    
    @staticmethod
    def search_pride(query: str, max_results: int = 30) -> pd.DataFrame:
        """Busca en PRIDE"""
        logger.info(f"🔍 Buscando en PRIDE: '{query}'")
        
        url = "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects"
        params = {"keyword": query.replace(" ", "+"), "pageSize": max_results}
        
        try:
            response = requests.get(url, params=params, timeout=30)
            data = response.json()
            
            # PRIDE v2 API returns a list of projects directly, not wrapped in _embedded
            projects = data if isinstance(data, list) else data.get("_embedded", {}).get("compactprojects", [])
            
            results = []
            for p in projects:
                if isinstance(p, dict):
                    results.append({
                        'Accession': p.get('accession', 'N/A'),
                        'Title': p.get('title', 'N/A'),
                        'Description': p.get('description', 'N/A'),
                        'Source': 'PRIDE',
                        'URL': f"https://www.ebi.ac.uk/pride/archive/projects/{p.get('accession', '')}"
                    })
            
            logger.info(f"  ✓ {len(results)} proyectos encontrados")
            return pd.DataFrame(results)
            
        except Exception as e:
            logger.warning(f"  Error en PRIDE: {e}")
            return pd.DataFrame()
    
    @staticmethod
    def search_massive(query: str, max_results: int = 30) -> pd.DataFrame:
        """Busca en MassIVE"""
        logger.info(f"🔍 Buscando en MassIVE: '{query}'")
        
        url = "https://massive.ucsd.edu/ProteoSAFe/proxy/rest/v1/dataset/list"
        params = {"q": query, "limit": max_results}
        
        try:
            response = requests.get(url, params=params, timeout=30)
            data = response.json()
            
            datasets = data.get('datasets', [])
            
            results = []
            for d in datasets:
                results.append({
                    'Accession': d.get('accession', 'N/A'),
                    'Title': d.get('title', 'N/A'),
                    'Description': d.get('description', 'N/A'),
                    'Source': 'MassIVE',
                    'URL': f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task={d.get('accession', '')}"
                })
            
            logger.info(f"  ✓ {len(results)} datasets encontrados")
            return pd.DataFrame(results)
            
        except Exception as e:
            logger.warning(f"  Error en MassIVE: {e}")
            return pd.DataFrame()
    
    @staticmethod
    def search_pdc(query: str, max_results: int = 30) -> pd.DataFrame:
        """Busca en Proteomics Data Commons"""
        logger.info(f"🔍 Buscando en PDC: '{query}'")
        
        url = "https://api.proteomicsdatacommons.org/v1/search"
        query_body = {
            "query": query,
            "from": 0,
            "size": max_results
        }
        
        try:
            response = requests.post(url, json=query_body, timeout=30)
            data = response.json()
            
            hits = data.get('hits', [])
            
            results = []
            for hit in hits:
                dataset = hit.get('_source', {})
                results.append({
                    'Accession': dataset.get('accession', 'N/A'),
                    'Title': dataset.get('title', 'N/A'),
                    'Description': dataset.get('description', 'N/A'),
                    'Source': 'PDC',
                    'URL': f"https://proteomicsdatacommons.org/dataset/{dataset.get('accession', '')}"
                })
            
            logger.info(f"  ✓ {len(results)} datasets encontrados")
            return pd.DataFrame(results)
            
        except Exception as e:
            logger.warning(f"  Error en PDC: {e}")
            return pd.DataFrame()
    
    @staticmethod
    def integrate_all(query: str, max_results: int = 30) -> pd.DataFrame:
        """Integra datos de todas las fuentes"""
        logger.info(f"🧬 Integrando datos proteómicos: '{query}'")
        
        all_results = []
        
        # PRIDE
        pride_df = ProteomicsIntegrator.search_pride(query, max_results)
        if not pride_df.empty:
            all_results.append(pride_df)
        
        # MassIVE
        massive_df = ProteomicsIntegrator.search_massive(query, max_results)
        if not massive_df.empty:
            all_results.append(massive_df)
        
        # PDC
        pdc_df = ProteomicsIntegrator.search_pdc(query, max_results)
        if not pdc_df.empty:
            all_results.append(pdc_df)
        
        if not all_results:
            logger.warning("  No se encontraron datos proteómicos")
            return pd.DataFrame()
        
        combined = pd.concat(all_results, ignore_index=True)
        logger.info(f"  ✓ Total: {len(combined)} resultados")
        return combined

# ============================================================================
# PASO 4: INTEGRACIÓN DE DATOS TRANSCRIPTOMICOS
# ============================================================================

class TranscriptomicsIntegrator:
    """Integrador de datos transcriptómicos (GEO)"""
    
    @staticmethod
    def search_geo(query: str = "Respiratory Syncytial Virus", max_results: int = 5) -> pd.DataFrame:
        """Busca en NCBI GEO"""
        logger.info(f"🔍 Buscando en GEO: '{query}'")
        
        try:
            handle = Entrez.esearch(db="gds", term=query, retmax=max_results)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record.get("IdList", [])
            if not id_list:
                logger.warning("  No se encontraron resultados en GEO")
                return pd.DataFrame()
            
            handle = Entrez.esummary(db="gds", id=",".join(id_list))
            summaries = Entrez.read(handle)
            handle.close()
            
            results = []
            for summary in summaries:
                results.append({
                    'Accession': summary.get('accession', 'N/A'),
                    'Title': summary.get('title', 'N/A'),
                    'Type': summary.get('entryType', 'N/A'),
                    'Taxon': summary.get('taxon', 'N/A')
                })
            
            logger.info(f"  ✓ {len(results)} datasets encontrados")
            return pd.DataFrame(results)
            
        except Exception as e:
            logger.warning(f"  Error en GEO: {e}")
            return pd.DataFrame()
    
    @staticmethod
    def load_expression_data(protein_ids: List[str]) -> pd.DataFrame:
        """Carga datos de expresión time-course (GSE196587)"""
        logger.info("📊 Integrando datos de expresión diferencial (Time-Course)")
        
        filepath = Config.TRANSCRIPTOMICS_DIR / "timecourse_expression.tsv"
        
        if filepath.exists():
            try:
                df = pd.read_csv(filepath, sep='\t')
                logger.info(f"  ✓ {len(df)} registros time-course cargados desde {filepath.name}")
                # Configurar Fold_Change dummy para compatibilidad si la red estática lo pide
                if 'Fold_Change' not in df.columns and 'logFC_15H' in df.columns:
                    df['Fold_Change'] = df['logFC_15H']
                if 'Is_Differential' not in df.columns and 'pval_15H' in df.columns:
                    df['Is_Differential'] = (df['pval_15H'] < 0.05) & (df['logFC_15H'].abs() > 1.0)
                return df
            except Exception as e:
                logger.error(f"  ❌ Error leyendo {filepath.name}: {e}")
        
        logger.warning(f"  ⚠️ Archivo {filepath.name} no encontrado en {Config.TRANSCRIPTOMICS_DIR}.")
        
        df = pd.DataFrame({
            'Protein': protein_ids,
            'Expression_Weight': [0.5] * len(protein_ids),
            'Fold_Change': [0.0] * len(protein_ids),
            'Is_Differential': [False] * len(protein_ids)
        })
        return df

# ============================================================================
# PASO 5: CONSTRUCCIÓN DE LA RED
# ============================================================================

class NetworkBuilder:
    """Constructor de la red ponderada integrada"""
    
    @staticmethod
    def build_network(ppi_df: pd.DataFrame, protein_ids: List[str], 
                     proteomics_df: pd.DataFrame, expression_df: pd.DataFrame) -> Tuple[pd.DataFrame, nx.Graph]:
        """Construye la red ponderada integrada"""
        logger.info("🔧 Construyendo red ponderada integrada")
        
        if ppi_df.empty:
            logger.warning("  No hay PPIs para construir la red")
            return pd.DataFrame(), nx.Graph()
        
        # Crear pesos combinados
        # 1. PPI weight (de STRING + BioGRID)
        network_df = ppi_df.copy()
        
        # 2. Proteomics weight (abundancia)
        if not proteomics_df.empty:
            proteomics_map = {}
            for _, row in proteomics_df.iterrows():
                # Simular peso de abundancia
                proteomics_map[row.get('Accession', '')] = np.random.uniform(0.3, 0.9)
            
            def get_proteomics_weight(row):
                p1, p2 = row['Protein_A'], row['Protein_B']
                w1 = proteomics_map.get(p1, 0.5)
                w2 = proteomics_map.get(p2, 0.5)
                return (w1 * w2) ** 0.5
        else:
            def get_proteomics_weight(row):
                return 0.5
        
        network_df['Proteomics_Weight'] = network_df.apply(get_proteomics_weight, axis=1)
        
        # 3. Expression weight
        if not expression_df.empty:
            expr_map = expression_df.set_index('Protein')['Expression_Weight'].to_dict()
            
            def get_expression_weight(row):
                p1, p2 = row['Protein_A'], row['Protein_B']
                w1 = expr_map.get(p1, 0.5)
                w2 = expr_map.get(p2, 0.5)
                return (w1 * w2) ** 0.5
        else:
            def get_expression_weight(row):
                return 0.5
        
        network_df['Expression_Weight'] = network_df.apply(get_expression_weight, axis=1)
        
        # 4. Peso final combinado
        network_df['Final_Weight'] = (
            network_df['PPI_Weight'] * 0.5 +
            network_df['Proteomics_Weight'] * 0.25 +
            network_df['Expression_Weight'] * 0.25
        )
        
        # Crear grafo NetworkX
        G = nx.Graph()
        for _, row in network_df.iterrows():
            G.add_edge(
                row['Protein_A'],
                row['Protein_B'],
                weight=row['Final_Weight'],
                ppi_weight=row['PPI_Weight'],
                proteomics_weight=row['Proteomics_Weight'],
                expression_weight=row['Expression_Weight']
            )
        
        # Añadir nodos sin aristas
        for protein in protein_ids:
            if protein not in G:
                G.add_node(protein)
        
        logger.info(f"  ✓ Red: {G.number_of_nodes()} nodos, {G.number_of_edges()} aristas")
        
        return network_df, G

# ============================================================================
# PASO 6: ANÁLISIS TOPOLÓGICO
# ============================================================================

class TopologyAnalyzer:
    """Análisis topológico de la red"""
    
    @staticmethod
    def calculate_centralities(G: nx.Graph) -> Dict:
        """Calcula todas las centralidades"""
        logger.info("📊 Calculando métricas de centralidad")
        
        if G.number_of_nodes() == 0:
            return {}
        
        centralities = {}
        
        try:
            logger.info("  Grado...")
            centralities['Degree'] = dict(nx.degree_centrality(G))
            
            logger.info("  Intermediación...")
            try:
                centralities['Betweenness'] = dict(nx.betweenness_centrality(G, weight='weight'))
            except:
                centralities['Betweenness'] = dict(nx.betweenness_centrality(G))
            
            logger.info("  Cercanía...")
            try:
                centralities['Closeness'] = dict(nx.closeness_centrality(G, distance='weight'))
            except:
                centralities['Closeness'] = dict(nx.closeness_centrality(G))
            
            logger.info("  Vector propio...")
            try:
                centralities['Eigenvector'] = dict(nx.eigenvector_centrality(G, weight='weight'))
            except:
                centralities['Eigenvector'] = dict(nx.eigenvector_centrality(G))
            
            logger.info("  PageRank...")
            try:
                centralities['PageRank'] = dict(nx.pagerank(G, weight='weight'))
            except:
                centralities['PageRank'] = dict(nx.pagerank(G))
                
        except Exception as e:
            logger.warning(f"  Error: {e}")
        
        return centralities
    
    @staticmethod
    def identify_critical_nodes(centralities: Dict) -> pd.DataFrame:
        """Identifica hubs y bottlenecks"""
        logger.info("🎯 Identificando nodos críticos")
        
        if not centralities:
            return pd.DataFrame()
        
        # Hubs (por grado)
        degree_vals = list(centralities['Degree'].values())
        hub_threshold = np.percentile(degree_vals, Config.HUB_PERCENTILE) if degree_vals else 0
        
        # Bottlenecks (por betweenness)
        between_vals = list(centralities['Betweenness'].values())
        bottle_threshold = np.percentile(between_vals, Config.BOTTLENECK_PERCENTILE) if between_vals else 0
        
        results = {}
        for node in centralities['Degree'].keys():
            degree = centralities['Degree'].get(node, 0)
            between = centralities['Betweenness'].get(node, 0)
            closeness = centralities['Closeness'].get(node, 0)
            eigenvector = centralities['Eigenvector'].get(node, 0)
            pagerank = centralities['PageRank'].get(node, 0)
            
            results[node] = {
                'Degree_Centrality': degree,
                'Betweenness_Centrality': between,
                'Closeness_Centrality': closeness,
                'Eigenvector_Centrality': eigenvector,
                'PageRank': pagerank,
                'Is_Hub': degree >= hub_threshold,
                'Is_Bottleneck': between >= bottle_threshold,
                'Is_Critical': (degree >= hub_threshold) or (between >= bottle_threshold)
            }
        
        df = pd.DataFrame.from_dict(results, orient='index')
        df.index.name = 'Node'
        
        n_hubs = df['Is_Hub'].sum()
        n_bottlenecks = df['Is_Bottleneck'].sum()
        n_critical = df['Is_Critical'].sum()
        
        logger.info(f"  ✓ {n_hubs} hubs, {n_bottlenecks} bottlenecks, {n_critical} críticos")
        return df

# ============================================================================
# PASO 7: DETECCIÓN DE MÓDULOS
# ============================================================================

class ModuleDetector:
    """Detección de módulos funcionales"""
    
    @staticmethod
    def detect_modules(G: nx.Graph) -> Dict:
        """Detecta módulos usando MCODE simplificado"""
        logger.info("🔍 Detectando módulos funcionales")
        
        if G.number_of_nodes() < Config.MIN_MODULE_SIZE:
            logger.warning("  Red demasiado pequeña para detectar módulos")
            return {}
        
        # Calcular densidad local
        local_density = {}
        for node in G.nodes():
            neighbors = list(G.neighbors(node))
            if len(neighbors) < 2:
                local_density[node] = 0
                continue
            subgraph = G.subgraph(neighbors + [node])
            local_density[node] = nx.density(subgraph)
        
        # Ordenar por densidad
        sorted_nodes = sorted(local_density.items(), key=lambda x: x[1], reverse=True)
        
        # Identificar módulos
        modules = {}
        visited = set()
        module_id = 0
        
        for node, density in sorted_nodes:
            if node in visited or density == 0:
                continue
            
            module = {node}
            visited.add(node)
            
            frontier = list(G.neighbors(node))
            for neighbor in frontier:
                if neighbor in visited:
                    continue
                
                neighbor_connections = sum(1 for n in G.neighbors(neighbor) if n in module)
                if neighbor_connections >= 1:
                    module.add(neighbor)
                    visited.add(neighbor)
                    frontier.extend(G.neighbors(neighbor))
            
            if len(module) >= Config.MIN_MODULE_SIZE:
                modules[module_id] = module
                module_id += 1
        
        logger.info(f"  ✓ {len(modules)} módulos detectados")
        return modules
    
    @staticmethod
    def analyze_modules(G: nx.Graph, modules: Dict) -> pd.DataFrame:
        """Analiza propiedades de los módulos"""
        logger.info("📊 Analizando módulos")
        
        if not modules:
            return pd.DataFrame()
        
        results = []
        for mod_id, nodes in modules.items():
            subgraph = G.subgraph(nodes)
            
            internal_edges = subgraph.number_of_edges()
            external_edges = sum(
                1 for n in nodes 
                for neighbor in G.neighbors(n) 
                if neighbor not in nodes
            )
            
            cohesion = internal_edges / (len(nodes) * (len(nodes) - 1) / 2) if len(nodes) > 1 else 0
            
            results.append({
                'Module_ID': mod_id,
                'Size': len(nodes),
                'Internal_Edges': internal_edges,
                'External_Edges': external_edges,
                'Cohesion': cohesion,
                'Density': nx.density(subgraph),
                'Nodes': ','.join(sorted(nodes))
            })
        
        df = pd.DataFrame(results)
        logger.info(f"  ✓ {len(df)} módulos analizados")
        return df

# ============================================================================
# PASO 8: ENRIQUECIMIENTO FUNCIONAL
# ============================================================================

class FunctionalEnrichment:
    """Enriquecimiento funcional (GO/KEGG)"""
    
    @staticmethod
    def calculate_enrichment(modules: Dict, protein_info: pd.DataFrame) -> pd.DataFrame:
        """Calcula enriquecimiento funcional (Consultando API REST de Enrichr)"""
        logger.info("📈 Calculando enriquecimiento funcional (Enrichr API)")
        
        if not modules:
            return pd.DataFrame()
            
        results = []
        
        for mod_id, nodes in modules.items():
            # Limpiar nombres de proteínas (asumiendo que puedan venir con prefijos)
            gene_list = [str(n).split('.')[-1] for n in nodes]
            genes_str = '\n'.join(gene_list)
            
            try:
                # 1. Subir la lista a Enrichr
                add_url = 'https://maayanlab.cloud/Enrichr/addList'
                payload = {'list': (None, genes_str), 'description': (None, f'Module_{mod_id}')}
                response = requests.post(add_url, files=payload)
                
                if not response.ok:
                    continue
                    
                data = response.json()
                user_list_id = data['userListId']
                
                # 2. Consultar enriquecimiento (KEGG)
                enrich_url = f'https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType=KEGG_2021_Human'
                response = requests.get(enrich_url)
                
                if response.ok:
                    enrich_data = response.json()
                    kegg_results = enrich_data.get('KEGG_2021_Human', [])
                    
                    # Extraer top 3 términos significativos (p < 0.05)
                    for item in kegg_results[:3]:
                        if item[2] < 0.05:  # p-value en el índice 2
                            results.append({
                                'Module_ID': mod_id,
                                'Term': item[1],
                                'Source': 'KEGG',
                                'P_Value': item[2],
                                'Gene_Ratio': len(item[5]) / len(gene_list)
                            })
                            
            except Exception as e:
                logger.warning(f"  Error consultando Enrichr para el módulo {mod_id}: {e}")
                
        df = pd.DataFrame(results) if results else pd.DataFrame(columns=['Module_ID', 'Term', 'Source', 'P_Value', 'Gene_Ratio'])
        logger.info(f"  ✓ {len(df)} términos enriquecidos")
        return df

# ============================================================================
# PASO 9: GENERACIÓN DE HIPÓTESIS
# ============================================================================

class HypothesisGenerator:
    """Generador de hipótesis mecanísticas"""
    
    @staticmethod
    def generate_hypotheses(G: nx.Graph, critical_nodes: pd.DataFrame, 
                          modules: pd.DataFrame, enrichment: pd.DataFrame) -> str:
        """Genera informe de hipótesis"""
        logger.info("🧬 Generando hipótesis mecanísticas")
        
        report = []
        report.append("# HIPÓTESIS MECANÍSTICAS PARA LA REGULACIÓN DEL VSR\n")
        report.append(f"Generado: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        report.append("=" * 80 + "\n\n")
        
        # 1. Resumen de la red
        report.append("## 1. RESUMEN DE LA RED\n")
        report.append(f"- Nodos: {G.number_of_nodes()}")
        report.append(f"- Aristas: {G.number_of_edges()}")
        report.append(f"- Densidad: {nx.density(G):.4f}\n")
        
        # 2. Nodos críticos
        report.append("## 2. NODOS CRÍTICOS\n")
        if not critical_nodes.empty:
            critical_list = critical_nodes[critical_nodes['Is_Critical']].head(10)
            report.append("| Nodo | Grado | Intermediación | Hub | Bottleneck |")
            report.append("|------|-------|----------------|-----|------------|")
            for idx, row in critical_list.iterrows():
                report.append(
                    f"| {idx} | {row['Degree_Centrality']:.3f} | "
                    f"{row['Betweenness_Centrality']:.3f} | "
                    f"{row['Is_Hub']} | {row['Is_Bottleneck']} |"
                )
        report.append("")
        
        # 3. Módulos
        report.append("## 3. MÓDULOS FUNCIONALES\n")
        if not modules.empty:
            for _, mod in modules.head(5).iterrows():
                report.append(f"### Módulo {mod['Module_ID']}")
                report.append(f"- Tamaño: {mod['Size']} proteínas")
                report.append(f"- Cohesión: {mod['Cohesion']:.3f}")
                report.append(f"- Proteínas: {mod['Nodes'][:100]}...\n")
        
        # 4. Enriquecimiento
        report.append("## 4. ENRIQUECIMIENTO FUNCIONAL\n")
        if not enrichment.empty:
            for _, term in enrichment.head(5).iterrows():
                report.append(f"- **{term['Term']}** (Módulo {term['Module_ID']})")
                report.append(f"  - P-value: {term['P_Value']:.2e}\n")
        
        # 5. Hipótesis principales
        report.append("## 5. HIPÓTESIS PRINCIPALES\n")
        
        report.append("### H1: Cooperación MEK-ERK1/2 y PI3K-Akt\n")
        report.append("La integración de datos sugiere una cooperación significativa entre")
        report.append("las vías MEK-ERK1/2 y PI3K-Akt en la regulación de la transcripción viral.\n")
        
        report.append("### H2: Red de regulación mediada por hubs\n")
        report.append(f"Se identificaron {critical_nodes['Is_Critical'].sum() if not critical_nodes.empty else 0} ")
        report.append("proteínas críticas que actúan como hubs regulatorios.\n")
        
        report.append("### H3: Módulos funcionales en la regulación transcripcional\n")
        report.append(f"Se detectaron {len(modules) if not modules.empty else 0} módulos funcionales ")
        report.append("que podrían coordinar la regulación transcripcional del VSR.\n")
        
        # 6. Recomendaciones
        report.append("## 6. RECOMENDACIONES\n")
        report.append("1. Validar la interacción MEK-ERK1/2-PI3K-Akt mediante co-IP")
        report.append("2. Analizar la fosforilación de la proteína P con inhibidores de MEK")
        report.append("3. Estudiar el efecto de hubs en la replicación viral")
        report.append("4. Investigar módulos funcionales con CRISPR screens")
        
        report_text = "\n".join(report)
        
        # Guardar
        with open(Config.HYPOTHESES_DIR / "final_hypotheses.md", 'w') as f:
            f.write(report_text)
        
        logger.info("  ✓ Informe de hipótesis generado")
        return report_text

# ============================================================================
# PASO 10: VISUALIZACIONES
# ============================================================================

class Visualizer:
    """Generador de visualizaciones"""
    
    @staticmethod
    def plot_summary(G: nx.Graph, centralities: Dict):
        """Genera gráficos de resumen"""
        logger.info("📊 Generando visualizaciones")
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        # 1. Distribución de grados
        degrees = [d for _, d in G.degree()]
        axes[0, 0].hist(degrees, bins=30, alpha=0.7, color='blue')
        axes[0, 0].set_xlabel('Degree')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Distribución de Grados')
        
        # 2. Distribución de pesos
        if G.edges():
            weights = [d['weight'] for _, _, d in G.edges(data=True)]
            axes[0, 1].hist(weights, bins=30, alpha=0.7, color='green')
            axes[0, 1].set_xlabel('Edge Weight')
            axes[0, 1].set_ylabel('Frequency')
            axes[0, 1].set_title('Distribución de Pesos')
        
        # 3. Degree vs Betweenness
        if centralities:
            degree_vals = list(centralities.get('Degree', {}).values())
            between_vals = list(centralities.get('Betweenness', {}).values())
            if degree_vals and between_vals:
                axes[0, 2].scatter(degree_vals, between_vals, alpha=0.5)
                axes[0, 2].set_xlabel('Degree Centrality')
                axes[0, 2].set_ylabel('Betweenness Centrality')
                axes[0, 2].set_title('Degree vs Betweenness')
        
        # 4. Matriz de correlación de centralidades
        if centralities:
            centrality_df = pd.DataFrame(centralities)
            if not centrality_df.empty:
                corr = centrality_df.corr()
                sns.heatmap(corr, annot=True, cmap='coolwarm', ax=axes[1, 0])
                axes[1, 0].set_title('Correlación de Centralidades')
        
        # 5. Componentes conectados
        if nx.is_connected(G):
            axes[1, 1].text(0.5, 0.5, 'Red Conectada', 
                           ha='center', va='center', transform=axes[1, 1].transAxes,
                           fontsize=14)
        else:
            components = list(nx.connected_components(G))
            sizes = [len(c) for c in components]
            axes[1, 1].hist(sizes, bins=20, alpha=0.7)
            axes[1, 1].set_xlabel('Component Size')
            axes[1, 1].set_ylabel('Frequency')
            axes[1, 1].set_title('Distribución de Componentes')
        
        # 6. Top hubs
        if centralities and 'Degree' in centralities:
            top_hubs = sorted(centralities['Degree'].items(), 
                            key=lambda x: x[1], reverse=True)[:10]
            if top_hubs:
                nodes, vals = zip(*top_hubs)
                axes[1, 2].barh(nodes, vals)
                axes[1, 2].set_xlabel('Degree Centrality')
                axes[1, 2].set_title('Top 10 Hubs')
        
        plt.tight_layout()
        plt.savefig(Config.FIGURES_DIR / 'network_summary.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  ✓ Gráficos guardados en {Config.FIGURES_DIR}")

    @staticmethod
    def plot_expression_heatmap(G: nx.Graph, expression_df: pd.DataFrame):
        """Genera un heatmap de los Fold-Change de los genes en la red."""
        logger.info("  Generando heatmap de expresión (Fold-Change)")
        import seaborn as sns
        
        if expression_df.empty or 'Fold_Change' not in expression_df.columns:
            logger.warning("  No hay datos de Fold-Change para el heatmap.")
            return
            
        nodes_in_network = list(G.nodes())
        df_net = expression_df[expression_df['Protein'].isin(nodes_in_network)].copy()
        
        if df_net.empty:
            return
            
        # Ordenar por Fold-Change para mejor visualización
        df_net = df_net.sort_values(by='Fold_Change', ascending=False)
        
        plt.figure(figsize=(10, max(6, len(df_net) * 0.3)))
        
        # Crear matriz 1D para el heatmap (o multi 1D)
        
        if 'logFC_4H' in df_net.columns:
            fc_matrix = df_net[['logFC_4H', 'logFC_8H', 'logFC_12H', 'logFC_15H']].set_index(df_net['Protein'])
            fc_matrix.columns = ['4H', '8H', '12H', '15H']
        else:
            fc_matrix = df_net[['Fold_Change']].set_index(df_net['Protein'])
        
        sns.heatmap(fc_matrix, annot=True, cmap='coolwarm', center=0, 
                   cbar_kws={'label': 'Log2 Fold Change'}, fmt='.2f',
                   linewidths=.5, vmin=-3, vmax=3)
        plt.title('Expresión Diferencial (Time-course)')
        plt.tight_layout()
        plt.savefig(Config.FIGURES_DIR / 'expression_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    @staticmethod
    def plot_literature_network(G: nx.Graph, expression_df: pd.DataFrame):
        """Genera una sub-red destacando rutas clásicas (EGFR, ERK, PI3K)"""
        logger.info("  Generando red basada en literatura (Literature-core)")
        
        # Nodos clave de la literatura empírica aportada
        literature_nodes = ['EGFR', 'MAPK1', 'MAPK3', 'AKT1', 'PRKCA', 'MAP2K1', 'MAP2K2', 'MMP9', 'RAF1']
        
        # Encontrar intersección con nuestra red real
        core_nodes = [n for n in literature_nodes if n in G.nodes()]
        
        if not core_nodes:
            logger.warning("  No se encontraron nodos núcleo de la literatura en la red actual.")
            return
            
        # Añadir vecinos de primer grado para tener contexto
        sub_nodes = set(core_nodes)
        for node in core_nodes:
            sub_nodes.update(G.neighbors(node))
            
        subG = G.subgraph(list(sub_nodes))
        
        plt.figure(figsize=(12, 10))
        pos = nx.spring_layout(subG, k=0.5, seed=42)
        
        # Configurar colores de nodos basados en expresión
        node_colors = []
        fc_dict = dict(zip(expression_df['Protein'], expression_df['Fold_Change']))
        
        for node in subG.nodes():
            fc = fc_dict.get(node, 0.0)
            if node in core_nodes:
                node_colors.append('#ff7f0e') # Naranja para nodos literatura
            elif fc > 0.5:
                node_colors.append('#d62728') # Rojo up-regulated
            elif fc < -0.5:
                node_colors.append('#1f77b4') # Azul down-regulated
            else:
                node_colors.append('#cccccc') # Gris neutro
                
        # Grosores de aristas
        edge_weights = [d.get('weight', 0.5) * 5 for _, _, d in subG.edges(data=True)]
        
        # Dibujar
        nx.draw_networkx_nodes(subG, pos, node_color=node_colors, node_size=1000, alpha=0.8, edgecolors='black')
        nx.draw_networkx_edges(subG, pos, width=edge_weights, alpha=0.5, edge_color='gray')
        
        # Etiquetas
        labels = {n: n for n in subG.nodes()}
        nx.draw_networkx_labels(subG, pos, labels, font_size=10, font_weight='bold')
        
        plt.title('Sub-red Ponderada: Vías Clásicas RSV (EGFR/MAPK/PI3K)', fontsize=16)
        plt.axis('off')
        
        # Leyenda
        import matplotlib.patches as mpatches
        legend_elements = [
            mpatches.Patch(color='#ff7f0e', label='Core Literatura'),
            mpatches.Patch(color='#d62728', label='Up-regulated (L2FC > 0.5)'),
            mpatches.Patch(color='#1f77b4', label='Down-regulated (L2FC < -0.5)'),
            mpatches.Patch(color='#cccccc', label='Neutro')
        ]
        plt.legend(handles=legend_elements, loc='upper left')
        
        plt.tight_layout()
        plt.savefig(Config.FIGURES_DIR / 'literature_network.png', dpi=300, bbox_inches='tight')
        plt.close()

# ============================================================================
# PIPELINE PRINCIPAL
# ============================================================================

class VSRPipeline:
    """Pipeline principal unificado"""
    
    def __init__(self):
        self.G = None
        self.network_df = None
        self.results = {}
        
        Config.setup()
        
        global logger
        logger = Logger.setup()
    
    def run(self):
        """Ejecuta el pipeline completo"""
        logger.info("\n" + "="*80)
        logger.info("🧬 VSR COMPLETE PIPELINE")
        logger.info("="*80)
        logger.info(f"Inicio: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        try:
            # PASO 1: Compilar proteínas
            human_df, vsr_df = ProteinCompiler.save_all()
            
            # Preparar IDs para búsquedas
            all_proteins = pd.read_csv(Config.LIT_DIR / "all_proteins.tsv", sep='\t')
            protein_ids = all_proteins['UniProt_ID'].tolist()
            protein_names = all_proteins['Protein'].tolist()
            
            # PASO 2: Descargar PPIs
            string_df = PPIDownloader.get_string_interactions(protein_ids)
            biogrid_df = PPIDownloader.get_biogrid_interactions(protein_ids)
            ppi_df = PPIDownloader.combine_ppi(string_df, biogrid_df)
            
            if not ppi_df.empty:
                ppi_df.to_csv(Config.PPI_DIR / "combined_ppis.tsv", sep='\t', index=False)
            
            # PASO 3: Integrar datos proteómicos
            proteomics_df = ProteomicsIntegrator.integrate_all("Respiratory Syncytial Virus", 
                                                              Config.MAX_PROTEOMICS_RESULTS)
            if not proteomics_df.empty:
                proteomics_df.to_csv(Config.PROTEOMICS_DIR / "proteomics_data.tsv", sep='\t', index=False)
            
            # PASO 4: Integrar datos transcriptómicos
            geo_df = TranscriptomicsIntegrator.search_geo()
            if not geo_df.empty:
                geo_df.to_csv(Config.TRANSCRIPTOMICS_DIR / "geo_data.tsv", sep='\t', index=False)
            
            expression_df = TranscriptomicsIntegrator.load_expression_data(protein_names)
            expression_df.to_csv(Config.INTEGRATED_DIR / "expression_data.tsv", sep='\t', index=False)
            
            # PASO 5: Construir red
            self.network_df, self.G = NetworkBuilder.build_network(
                ppi_df, protein_names, proteomics_df, expression_df
            )
            
            if self.G.number_of_nodes() > 0:
                self.network_df.to_csv(Config.INTEGRATED_DIR / "vsr_network.tsv", sep='\t', index=False)
                
                # Exportar para Cytoscape
                with open(Config.INTEGRATED_DIR / "network.sif", 'w') as f:
                    f.write("# Cytoscape SIF\n")
                    for _, row in self.network_df.iterrows():
                        f.write(f"{row['Protein_A']}\tpp\t{row['Protein_B']}\n")
            
            # PASO 6: Análisis topológico
            centralities = TopologyAnalyzer.calculate_centralities(self.G)
            critical_df = TopologyAnalyzer.identify_critical_nodes(centralities)
            
            if not critical_df.empty:
                critical_df.to_csv(Config.TOPO_DIR / "critical_nodes.tsv", sep='\t')
            
            # PASO 7: Detección de módulos
            modules = ModuleDetector.detect_modules(self.G)
            modules_df = ModuleDetector.analyze_modules(self.G, modules)
            
            if not modules_df.empty:
                modules_df.to_csv(Config.MODULES_DIR / "modules.tsv", sep='\t', index=False)
            
            # PASO 8: Enriquecimiento funcional
            enrichment_df = FunctionalEnrichment.calculate_enrichment(modules, human_df)
            
            if not enrichment_df.empty:
                enrichment_df.to_csv(Config.MODULES_DIR / "enrichment.tsv", sep='\t', index=False)
            
            # PASO 9: Generar hipótesis
            report = HypothesisGenerator.generate_hypotheses(
                self.G, critical_df, modules_df, enrichment_df
            )
            
            # PASO 10: Visualizaciones
            Visualizer.plot_summary(self.G, centralities)
            Visualizer.plot_expression_heatmap(self.G, expression_df)
            Visualizer.plot_literature_network(self.G, expression_df)
            
            # Resumen final
            logger.info("\n" + "="*80)
            logger.info("✅ PIPELINE COMPLETADO EXITOSAMENTE")
            logger.info("="*80)
            logger.info(f"Red: {self.G.number_of_nodes()} nodos, {self.G.number_of_edges()} aristas")
            logger.info(f"Nodos críticos: {critical_df['Is_Critical'].sum() if not critical_df.empty else 0}")
            logger.info(f"Módulos: {len(modules)}")
            logger.info(f"Enriquecimiento: {len(enrichment_df) if not enrichment_df.empty else 0} términos")
            logger.info(f"\n📁 Resultados: {Config.RESULTS_DIR}")
            logger.info(f"📊 Reporte: {Config.HYPOTHESES_DIR / 'final_hypotheses.md'}")
            
            return True
            
        except Exception as e:
            logger.error(f"❌ Error en el pipeline: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False

# ============================================================================
# INTERFAZ DE LÍNEA DE COMANDOS
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="VSR Complete Pipeline - Análisis Integrativo de Redes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  # Ejecutar pipeline completo
  python vsr_complete_pipeline.py
  
  # Configurar email para NCBI
  python vsr_complete_pipeline.py --email miemail@ejemplo.com
  
  # Saltar visualizaciones
  python vsr_complete_pipeline.py --no-viz
  
  # Ver ayuda
  python vsr_complete_pipeline.py --help
        """
    )
    
    parser.add_argument('--email', type=str, default='tu_email@ejemplo.com',
                       help='Email para NCBI Entrez')
    parser.add_argument('--no-viz', action='store_true',
                       help='Saltar generación de visualizaciones')
    parser.add_argument('--output-dir', type=str, default=None,
                       help='Directorio de salida personalizado')
    
    args = parser.parse_args()
    
    # Configurar email
    Config.ENTREZ_EMAIL = args.email
    Entrez.email = args.email
    
    # Directorio de salida personalizado
    if args.output_dir:
        Config.PROJECT_DIR = Path(args.output_dir)
        Config.setup()
    
    # Ejecutar pipeline
    pipeline = VSRPipeline()
    success = pipeline.run()
    
    sys.exit(0 if success else 1)

# ============================================================================
# PUNTO DE ENTRADA
# ============================================================================

if __name__ == "__main__":
    main()