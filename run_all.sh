#!/bin/bash
# ===============================================================================
# PIPELINE MAESTRO: ANÁLISIS INTEGRATIVO VSR MULTI-ÓMICA
# ===============================================================================
# Este script orquesta la ejecución completa de todos los módulos bioinformáticos,
# desde la adquisición de datos crudos hasta la integración multi-ómica.
# ===============================================================================

set -e # Detener el script si hay algún error
set -u # Detener si se usa una variable no definida

echo "==============================================================================="
echo "🧬 INICIANDO PIPELINE MULTI-ÓMICO VSR"
echo "==============================================================================="

# Asegurar directorios base
mkdir -p data results logs scripts

# -------------------------------------------------------------------------------
# FASE 1: ADQUISICIÓN DE DATOS Y METADATOS
# -------------------------------------------------------------------------------
echo ""
echo "[FASE 1] ADQUISICIÓN DE DATOS (Data Fetching)"
echo "-------------------------------------------------------------------------------"

echo "[1/3] Descargando catálogo de estudios GEO..."
python3 scripts/fetch_rsv_geo.py

echo "[2/3] Descargando genomas virales defectivos (ViReMa) de GSE281185..."
python3 scripts/analyze_dvgs.py

echo "[3/3] Descargando y preparando dataset scRNA-seq (GSE281623)..."
python3 scripts/run_scrnaseq.py

# -------------------------------------------------------------------------------
# FASE 2: PROCESAMIENTO CORE (Bulk & Redes Estáticas)
# -------------------------------------------------------------------------------
echo ""
echo "[FASE 2] PROCESAMIENTO CORE Y REDES (Processing)"
echo "-------------------------------------------------------------------------------"
echo "[1/3] Infiriendo interactomas temporales (12h a 120h) de GSE247298..."
# Asume que los archivos .tsv.gz ya fueron descomprimidos de GSE247298_RAW.tar
python3 scripts/infer_timecourse_interactome.py

echo "[2/3] Calculando módulos de coexpresión (WGCNA) sobre GSE247298..."
python3 scripts/run_wgcna.py
echo "[3/3] Analizando Interacciones Huésped-Patógeno (Interactómica)..."
python3 scripts/host_pathogen_interactomics.py


echo "[Bonus] Formulando Hipótesis Mecanísticas (Objetivos 1 y 4)..."
python3 scripts/formulate_hypotheses.py

echo "[Bonus] Generando visualizaciones previas de las Redes..."
python3 scripts/plot_simplified_networks.py

# -------------------------------------------------------------------------------
# FASE 3: INTEGRACIÓN Y VISUALIZACIÓN FINAL
# -------------------------------------------------------------------------------
echo ""
echo "[FASE 3] INTEGRACIÓN Y VISUALIZACIÓN"
echo "-------------------------------------------------------------------------------"

echo "[1/2] Generando panel dual integrado (Cinética Huésped vs. Virus)..."
python3 scripts/plot_integrated_kinetics.py

echo "[2/2] Construyendo Dashboard HTML Final..."
python3 scripts/generate_master_report.py

echo ""
echo "==============================================================================="
echo "✅ PIPELINE COMPLETADO EXITOSAMENTE"
echo "Resultados principales ubicados en:"
echo " - Redes e interacciones: results/interactomes/"
echo " - Módulos de Coexpresión: results/wgcna/"
echo " - Resumen de Mutaciones Virales: results/dvgs/"
echo " - Figuras listas para publicación: results/figures/"
echo "==============================================================================="
