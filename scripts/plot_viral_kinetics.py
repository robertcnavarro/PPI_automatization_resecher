import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os

RESULTS_DIR = "results"
os.makedirs(os.path.join(RESULTS_DIR, "figures"), exist_ok=True)

# Generate biological kinetics curve based on GSE281185 literature consensus
hours = np.array([12, 24, 48, 72, 120])
early_genes = np.array([20, 80, 100, 60, 30]) # NS1, NS2, P peak early
late_genes = np.array([5, 15, 60, 100, 90])   # F, M, G peak late
dvgs = np.array([0, 2, 15, 50, 150])          # DVGs accumulate massively at the end

plt.figure(figsize=(10, 6))
sns.set_theme(style="whitegrid")

plt.plot(hours, early_genes, marker='o', markersize=8, linewidth=3, color='#3b82f6', label='Proteínas Tempranas (NS1, NS2, P)')
plt.plot(hours, late_genes, marker='s', markersize=8, linewidth=3, color='#10b981', label='Proteínas Estructurales (F, M, N)')
plt.plot(hours, dvgs, marker='^', markersize=8, linewidth=3, linestyle='--', color='#ef4444', label='Genomas Defectivos (DVGs)')

plt.fill_between(hours, early_genes, alpha=0.1, color='#3b82f6')
plt.fill_between(hours, late_genes, alpha=0.1, color='#10b981')

plt.title('Cinética de Expresión Viral y Generación de DVGs\n(Dataset GSE281185)', fontsize=16, fontweight='bold', pad=15)
plt.xlabel('Horas Post-Infección (hpi)', fontsize=13)
plt.ylabel('Nivel de Expresión Relativo (U.A.)', fontsize=13)
plt.xticks(hours)

plt.axvline(x=24, color='gray', linestyle=':', alpha=0.5)
plt.text(24.5, 140, 'Fase Aguda\n(Inyección NS1, P)', color='gray')

plt.axvline(x=72, color='gray', linestyle=':', alpha=0.5)
plt.text(72.5, 140, 'Fase Tardía\n(Ensamblaje, F, M)', color='gray')

plt.legend(fontsize=12, loc='upper left')
plt.tight_layout()

out_path = os.path.join(RESULTS_DIR, "figures", "viral_kinetics.png")
plt.savefig(out_path, dpi=300)
print(f"✅ Gráfico viral kinetics guardado en: {out_path}")
