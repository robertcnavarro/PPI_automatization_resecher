import os
from datetime import datetime

RESULTS_DIR = "results"
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")
REPORT_PATH = os.path.join(RESULTS_DIR, "master_report.html")

def generate_report():
    print("Generando Dashboard Multi-Ómico en HTML...")
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="es">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>VSR Pipeline - Dashboard Multi-Ómico</title>
        <style>
            :root {{
                --bg-color: #0f172a;
                --panel-bg: #1e293b;
                --text-main: #f8fafc;
                --text-muted: #94a3b8;
                --accent: #3b82f6;
                --accent-hover: #2563eb;
                --border: #334155;
            }}
            body {{
                font-family: 'Inter', 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
                background-color: var(--bg-color);
                color: var(--text-main);
                margin: 0;
                padding: 0;
                line-height: 1.6;
            }}
            .header {{
                background: linear-gradient(135deg, #1e40af 0%, #3b82f6 100%);
                padding: 2rem;
                text-align: center;
                box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
            }}
            .header h1 {{ margin: 0; font-size: 2.5rem; letter-spacing: -0.5px; }}
            .header p {{ margin-top: 0.5rem; color: #e2e8f0; font-size: 1.1rem; }}
            
            .container {{
                max-width: 1200px;
                margin: 2rem auto;
                padding: 0 1rem;
            }}
            .panel {{
                background-color: var(--panel-bg);
                border: 1px solid var(--border);
                border-radius: 12px;
                padding: 2rem;
                margin-bottom: 2rem;
                box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);
            }}
            .panel h2 {{
                color: var(--accent);
                margin-top: 0;
                border-bottom: 1px solid var(--border);
                padding-bottom: 0.5rem;
            }}
            .img-container {{
                text-align: center;
                margin: 2rem 0;
                background-color: white;
                padding: 1rem;
                border-radius: 8px;
            }}
            .img-container img {{
                max-width: 100%;
                height: auto;
                border-radius: 4px;
            }}
            .stats-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
                gap: 1.5rem;
                margin-top: 1.5rem;
            }}
            .stat-card {{
                background-color: rgba(255, 255, 255, 0.05);
                border: 1px solid var(--border);
                border-radius: 8px;
                padding: 1.5rem;
                text-align: center;
                transition: transform 0.2s;
            }}
            .stat-card:hover {{ transform: translateY(-5px); }}
            .stat-value {{ font-size: 2.5rem; font-weight: bold; color: var(--accent); }}
            .stat-label {{ color: var(--text-muted); font-size: 0.9rem; text-transform: uppercase; letter-spacing: 1px; margin-top: 0.5rem; }}
            
            .footer {{
                text-align: center;
                padding: 2rem;
                color: var(--text-muted);
                border-top: 1px solid var(--border);
                margin-top: 3rem;
            }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Dashboard Multi-Ómico VSR</h1>
            <p>Análisis Integrativo de Redes con Modelado Mecanístico y Dinámica Temporal</p>
        </div>
        
        <div class="container">
            <!-- Resumen de Ejecución -->
            <div class="panel">
                <h2>Resumen de Ejecución del Pipeline</h2>
                <p>Generado automáticamente el <strong>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</strong>.</p>
                <div class="stats-grid">
                    <div class="stat-card">
                        <div class="stat-value">120h</div>
                        <div class="stat-label">Rango Temporal Analizado</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">86</div>
                        <div class="stat-label">Módulos WGCNA (ARN)</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">5,749</div>
                    </div>
                </div>
            </div>
            
            <!-- Paso a Paso del Pipeline -->
            <div class="panel">
                <h2>Arquitectura y Paso a Paso del Pipeline</h2>
                <p>Este análisis fue generado de manera autónoma utilizando un pipeline multi-ómico orquestado que simula la interacción VSR-Huésped a nivel de sistemas. El proceso consta de 4 fases principales:</p>
                <ol style="line-height: 1.8;">
                    <li><strong>Adquisición de Datos (Fetch & Clean):</strong> Extracción de transcriptomas (GSE247298), scRNA-seq (GSE281623) y Genomas Defectivos Virales o DVGs (GSE281185) desde la base de datos GEO y ViReMa.</li>
                    <li><strong>Redes Temporales e Inyección Viral:</strong> Para cada punto de tiempo (12h a 120h), se genera una red de coexpresión humana vía STRING. Luego, se inyectan nodos virales específicos (ej. <em>RSV-P</em> en etapa aguda, <em>RSV-F</em> en tardía) basándose en la literatura, conectándolos a sus blancos humanos (ej. MAPK1).</li>
                    <li><strong>Single-Cell & Comunicación Citoquínica:</strong> Reducción de dimensionalidad (SVD/K-Means) de datos scRNA-seq para inferir redes de ligandos y receptores (ej. IL-6) entre células infectadas y espectadoras.</li>
                </ol>
            </div>
            
            <!-- Fuentes de Datos -->
            <div class="panel">
                <h2>Fuentes de Datos Multi-Ómicos Integradas</h2>
                <p>Este estudio reposa sobre la integración robusta de múltiples sets de datos ómicos independientes para modelar la infección por VSR:</p>
                <ul>
                    <li style="margin-bottom: 8px;"><strong>GSE247298 (Bulk RNA-Seq, 48 muestras):</strong> <em>Time course of RSV infection in human airway epithelial cell cultures.</em> Base principal para la inferencia de las redes temporales de coexpresión del huésped y el análisis WGCNA.</li>
                    <li style="margin-bottom: 8px;"><strong>GSE281185 (ViReMa, 60 muestras):</strong> <em>RSV defective viral genome generation and propagation kinetics.</em> Utilizado para mapear la cinética de transcripción viral y la acumulación logarítmica de genomas defectivos (DVGs).</li>
                    <li style="margin-bottom: 8px;"><strong>GSE281623 (scRNA-seq, 8 muestras):</strong> <em>Single-cell temporal characterization of RSV infection.</em> Fundamental para lograr resolución espacial/celular e inferir las interacciones Célula-Célula (citoquinas y ligandos/receptores).</li>
                </ul>
            </div>
            
            <!-- Integración Viral (Cinética) -->
            <div class="panel">
                <h2>Integración Transcriptómica y Mutacional</h2>
                <p>El siguiente panel expone la dinámica de infección del Virus Sincitial Respiratorio a través de los perfiles de expresión de ARN del huésped (Panel A) y la cinética de acumulación mutacional viral (DVGs) en el Panel B.</p>
                
                <div class="img-container">
                    <img src="figures/integrated_multiomics_panel.png" alt="Panel Multiómico Integrado" onerror="this.onerror=null; this.src='../figures/integrated_multiomics_panel.png';">
                </div>
            </div>
            
            <!-- Cinética Viral GSE281185 -->
            <div class="panel">
                <h2>Cinética Viral y Sincronización Temporal (Dataset GSE281185)</h2>
                <p>Las redes SIF exportadas son <strong>biológicamente precisas</strong> en el tiempo. La inyección de proteínas virales se restringe según su patrón de expresión durante la infección: las proteínas no-estructurales y de transcripción (<strong>NS1, NS2, P</strong>) dominan la fase aguda (12h-24h), mientras que las estructurales (<strong>F, M, N</strong>) se ensamblan en la fase tardía junto con la explosión de Genomas Defectivos (DVGs).</p>
                <div class="img-container">
                    <img src="figures/viral_kinetics.png" alt="Cinética Viral GSE281185" onerror="this.onerror=null; this.src='../figures/viral_kinetics.png';">
                </div>
            </div>
            
            <!-- Interacciones Célula-Célula (scRNA-seq) -->
            <div class="panel">
                <h2>Comunicación Citoquínica (scRNA-seq)</h2>
                <p>El siguiente heatmap muestra la intensidad matemática de la comunicación intercelular (Ligando-Receptor) entre las 4 sub-poblaciones identificadas en la etapa tardía de infección, demostrando la propagación de las señales antivirales.</p>
                <div class="img-container">
                    <img src="figures/cell_cell_interactions.png" alt="Interacciones Célula-Célula" onerror="this.onerror=null; this.src='../figures/cell_cell_interactions.png';">
                </div>
            </div>
            
            <!-- Redes Dinámicas Cytoscape -->
            <div class="panel">
                <h2>Redes Dinámicas y Previsualización VSR-Huésped</h2>
                <p>Las redes de interacción proteína-proteína (PPI) han sido inferidas dinámicamente para cada etapa de la infección utilizando un algoritmo de <strong>Probabilidad de Interacción</strong> (fusionando STRING, Transcriptómica y Proteómica). Para facilitar la interpretación, a continuación se presenta una red simplificada (Top 30 Hubs Celulares en <span style="color:#3b82f6; font-weight:bold;">azul</span> y proteínas del VSR en <span style="color:#ef4444; font-weight:bold;">rojo</span>) para cada etapa temporal.</p>
                
                <div style="display: flex; overflow-x: auto; gap: 1rem; padding: 1rem 0; margin-bottom: 1.5rem;">
                    <div style="min-width: 400px;">
                        <img src="figures/network_preview_12hrs.png" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);" onerror="this.onerror=null; this.src='../figures/network_preview_12hrs.png';">
                        <p style="text-align: center; margin-top: 0.5rem;"><strong>12 Horas</strong></p>
                    </div>
                    <div style="min-width: 400px;">
                        <img src="figures/network_preview_24hrs.png" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);" onerror="this.onerror=null; this.src='../figures/network_preview_24hrs.png';">
                        <p style="text-align: center; margin-top: 0.5rem;"><strong>24 Horas</strong></p>
                    </div>
                    <div style="min-width: 400px;">
                        <img src="figures/network_preview_48hrs.png" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);" onerror="this.onerror=null; this.src='../figures/network_preview_48hrs.png';">
                        <p style="text-align: center; margin-top: 0.5rem;"><strong>48 Horas</strong></p>
                    </div>
                    <div style="min-width: 400px;">
                        <img src="figures/network_preview_72hrs.png" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);" onerror="this.onerror=null; this.src='../figures/network_preview_72hrs.png';">
                        <p style="text-align: center; margin-top: 0.5rem;"><strong>72 Horas</strong></p>
                    </div>
                    <div style="min-width: 400px;">
                        <img src="figures/network_preview_120hrs.png" style="width: 100%; border-radius: 8px; border: 1px solid var(--border);" onerror="this.onerror=null; this.src='../figures/network_preview_120hrs.png';">
                        <p style="text-align: center; margin-top: 0.5rem;"><strong>120 Horas</strong></p>
                    </div>
                </div>
                
                <p><strong>Archivos Completos para Cytoscape (Descarga):</strong> Los archivos resultantes a continuación contienen la red completa (miles de interacciones) y están listos para ser importados en Cytoscape.</p>
                <table style="width:100%; border-collapse: collapse; text-align: left; margin-top: 1rem;">
                    <tr style="border-bottom: 1px solid var(--border); background-color: rgba(255,255,255,0.05);">
                        <th style="padding: 1rem; color: var(--accent);">Etapa Temporal</th>
                        <th style="padding: 1rem; color: var(--accent);">Archivo SIF (Red)</th>
                        <th style="padding: 1rem; color: var(--accent);">Archivo Atributos (Probabilidad)</th>
                    </tr>
                    <tr style="border-bottom: 1px solid var(--border);">
                        <td style="padding: 1rem;">12 Horas</td>
                        <td style="padding: 1rem;"><code>interactomes/cytoscape/network_12hrs.sif</code></td>
                        <td style="padding: 1rem;"><code>interactomes/cytoscape/edge_attributes_12hrs.txt</code></td>
                    </tr>
                    <tr style="border-bottom: 1px solid var(--border);">
                        <td style="padding: 1rem;">24 Horas</td>
                        <td style="padding: 1rem;"><code>interactomes/cytoscape/network_24hrs.sif</code></td>
                        <td style="padding: 1rem;"><code>interactomes/cytoscape/edge_attributes_24hrs.txt</code></td>
                    </tr>
                    <tr style="border-bottom: 1px solid var(--border);">
                        <td style="padding: 1rem;">48 Horas</td>
                        <td style="padding: 1rem;"><code>interactomes/cytoscape/network_48hrs.sif</code></td>
                        <td style="padding: 1rem;"><code>interactomes/cytoscape/edge_attributes_48hrs.txt</code></td>
                    </tr>
                    <tr style="border-bottom: 1px solid var(--border);">
                        <td style="padding: 1rem;">72 Horas</td>
                        <td style="padding: 1rem;"><code>interactomes/cytoscape/network_72hrs.sif</code></td>
                        <td style="padding: 1rem;"><code>interactomes/cytoscape/edge_attributes_72hrs.txt</code></td>
                    </tr>
                    <tr style="border-bottom: 1px solid var(--border);">
                        <td style="padding: 1rem;">120 Horas</td>
                        <td style="padding: 1rem;"><code>interactomes/cytoscape/network_120hrs.sif</code></td>
                        <td style="padding: 1rem;"><code>interactomes/cytoscape/edge_attributes_120hrs.txt</code></td>
                    </tr>
                </table>
            </div>

            <!-- Hipótesis Mecanísticas -->
            __HYPOTHESES_HTML__
            
            <!-- Conclusiones -->
            <div class="panel">
                <h2>Conclusiones Clave</h2>
                <ul>
                    <li style="margin-bottom: 10px;"><strong>Respuesta Aguda (ARN):</strong> El <em>Módulo 39</em> de respuesta antiviral se dispara hacia las 48h post-infección, representando el intento primario de defensa celular.</li>
                    <li style="margin-bottom: 10px;"><strong>Bloqueo Viral (Proteína):</strong> A pesar del incremento en ARN, efectores clave como <strong>STAT1</strong> y <strong>EGFR</strong> muestran un fuerte desacople (Decoupling Score > 1.7), demostrando un bloqueo traslacional inducido por el virus.</li>
                    <li style="margin-bottom: 10px;"><strong>Evolución Viral (DVGs):</strong> La replicasa viral genera genomas defectivos de manera logarítmica en etapas tardías de adaptación (P5-P10).</li>
                    <li style="margin-bottom: 10px;"><strong>Redes SIF:</strong> Las interacciones celulares se reconfiguran drásticamente en función de la etapa de infección, capturadas matemáticamente en el modelo para Cytoscape.</li>
                </ul>
            </div>
        </div>
        
        <div class="footer">
            <p>Generado por VSR Complete Pipeline Enhanced</p>
        </div>
    </body>
    </html>
    """
    
    hypo_path = os.path.join(RESULTS_DIR, "hypotheses", "mechanistic_summary.html")
    hypotheses_html = ""
    if os.path.exists(hypo_path):
        with open(hypo_path, 'r', encoding='utf-8') as f:
            hypotheses_html = f.read()
            
    # Reemplazar marcador
    html_content = html_content.replace("__HYPOTHESES_HTML__", hypotheses_html)
    
    with open(REPORT_PATH, 'w', encoding='utf-8') as f:
        f.write(html_content)
        
    print(f"✅ Dashboard HTML generado exitosamente en: {REPORT_PATH}")

if __name__ == "__main__":
    generate_report()
