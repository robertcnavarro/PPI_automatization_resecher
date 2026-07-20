# RSV Multi-Omics Interactomics Pipeline

Este repositorio contiene un pipeline bioinformático completamente automatizado para el análisis de redes de señalización celular y la interacción huésped-patógeno durante la infección por el **Virus Sincitial Respiratorio (VSR)**.

El objetivo principal es integrar datos experimentales masivos (Proteómica y Transcriptómica) para inferir la actividad de Factores de Transcripción (TFs) dominantes e identificar cuellos de botella topológicos que puedan servir como blancos terapéuticos.

## 🧬 Características Principales

1. **Ingesta Multi-Ómica:** Procesamiento directo de matrices de cuantificación generadas por **DIA-NN** / MaxQuant, así como automatización de búsquedas en bases de datos públicas (PRIDE Archive y NCBI GEO).
2. **Inferencia de Factores de Transcripción (TFs):** Uso de `decoupleR` (Univariate Linear Model - ULM) cruzando la expresión experimental contra redes de conocimiento transcripcional validadas biológicamente (**CollecTRI** y **DoRothEA**).
3. **Análisis Topológico de Redes:** Construcción de interactomas masivos usando `NetworkX`, con cálculo automatizado de métricas topológicas (ej. *Betweenness Centrality* y *Degree*) para identificar los "cuellos de botella" (Hubs) de la infección.
4. **Exportación y Renderizado (Cytoscape):** Conexión directa mediante la API REST (`py4cytoscape`) para dibujar redes visualmente impactantes, con mapeo continuo de log2FC en el color de los nodos, grosor de aristas según confianza estadística, y tamaño de nodo escalado por centralidad matemática.
5. **Orquestación Global:** Scripts adicionales para análisis dinámico a lo largo del tiempo (Time-Course) y co-expresión génica con **WGCNA**.

## 🛠️ Requisitos Previos

Para ejecutar el pipeline necesitas tener instalados los siguientes entornos:

* **Python 3.9+** con las librerías: `pandas`, `numpy`, `networkx`, `decoupler`, `py4cytoscape`, `biopython`, `requests`.
* **Cytoscape Desktop Application** (debe estar abierto y ejecutándose en segundo plano para la Fase 6).

## 🚀 Uso del Pipeline

El núcleo del análisis interactómico se encuentra en el script principal de Python.

1. Abre **Cytoscape** en tu ordenador.
2. Ejecuta el script principal de proteogenómica:

```bash
python scripts/host_pathogen_interactomics.py
```

### Orquestación Completa
Si deseas correr todo el flujo de trabajo (incluyendo inferencias temporales y módulos de WGCNA), puedes ejecutar el bash script maestro:
```bash
./run_all.sh
```

## 📊 Arquitectura del Pipeline (`host_pathogen_interactomics.py`)

* **Fase 1 & 2:** Búsqueda en bases de datos (PRIDE) y extracción de genomas virales (NCBI).
* **Fase 3:** Control de Calidad y limpieza de la matriz de cuantificación proteómica (ej. DIA-NN output).
* **Fase 4:** Inferencia estadística con `decoupleR`. Se carga la red transcripcional local (para evitar caídas de red) y se extraen los **Top 10 TFs** más activos en la respuesta inmune (ej. *STAT2, STAT1, IRF9, RELA*).
* **Fase 5:** Modelado Matemático con `NetworkX`. Se calcula el *Betweenness Centrality* para medir la importancia biológica de las proteínas cuantificadas.
* **Fase 6:** Mapeo Visual. Se envía el grafo a Cytoscape por medio de CyREST, exportando la figura de publicación en alta calidad (`PNG` a 300 DPI).

## 📝 Licencia y Autores
Este proyecto forma parte de la investigación automatizada de redes de Interacción Proteína-Proteína (PPI) aplicadas a virología clínica.
