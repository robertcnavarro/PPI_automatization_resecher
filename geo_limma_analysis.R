# Instalar dependencias si no existen
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="http://cran.us.r-project.org")
if (!requireNamespace("GEOquery", quietly = TRUE))
    BiocManager::install("GEOquery", update=FALSE, ask=FALSE)
if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma", update=FALSE, ask=FALSE)
if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr", repos="http://cran.us.r-project.org")

library(GEOquery)
library(limma)
library(dplyr)

# 1. Descargar dataset GSE42026
cat("Descargando dataset GSE42026 desde NCBI GEO...\n")
gse <- getGEO("GSE42026", GSEMatrix = TRUE)
if (length(gse) > 1) {
  idx <- grep("GPL6244", attr(gse, "names")) # Affymetrix Human Gene 1.0 ST Array
  if(length(idx) == 0) idx <- 1
  gse <- gse[[idx]]
} else {
  gse <- gse[[1]]
}

# 2. Extraer expresión y metadata
exprs_mat <- exprs(gse)
pdata <- pData(gse)

# 3. Definir grupos (Mock vs RSV)
cat("Procesando matriz de diseño experimental...\n")
# Extraer la descripción fenotípica (ajustar de ser necesario para otros datasets)
conditions <- as.character(pdata$title)
groups <- ifelse(grepl("mock|control|uninfected", conditions, ignore.case=TRUE), "Mock", "RSV")
groups <- factor(groups, levels=c("Mock", "RSV"))

design <- model.matrix(~ 0 + groups)
colnames(design) <- c("Mock", "RSV")

# 4. Análisis de expresión diferencial (limma)
cat("Ejecutando análisis estadístico con limma...\n")
fit <- lmFit(exprs_mat, design)
contrast.matrix <- makeContrasts(RSV_vs_Mock = RSV - Mock, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Obtener tabla de resultados top
results <- topTable(fit2, adjust.method="fdr", number=Inf)

# 5. Mapear sondas a genes
cat("Mapeando IDs a Símbolos Genéticos...\n")
fdata <- fData(gse)

if("gene_assignment" %in% colnames(fdata)) {
  # Extraer Gene Symbol del campo gene_assignment de Affymetrix
  results$Protein <- sapply(strsplit(as.character(fdata[rownames(results), "gene_assignment"]), " // "), function(x) x[2])
} else if("Gene Symbol" %in% colnames(fdata)) {
  results$Protein <- fdata[rownames(results), "Gene Symbol"]
} else if("Symbol" %in% colnames(fdata)) {
  results$Protein <- fdata[rownames(results), "Symbol"]
} else {
  # Fallback: usar nombres de las sondas
  results$Protein <- rownames(results)
}

# 6. Formatear y exportar para Python (rsv_network_pipeline.py)
cat("Formateando resultados para pipeline de Python...\n")

# Filtrar genes válidos
final_results <- results %>%
  filter(!is.na(Protein) & Protein != "" & Protein != "---") %>%
  group_by(Protein) %>%
  # Si hay múltiples sondas para un gen, tomar la del menor P-Value
  arrange(adj.P.Val) %>%
  slice(1) %>%
  ungroup()

# Renombrar y calcular columnas requeridas
final_results <- final_results %>%
  mutate(
    Fold_Change = logFC,
    Is_Differential = (adj.P.Val < 0.05) & (abs(Fold_Change) > 1.0),
    # Weight entre 0 y 1 derivado del p-value (menor p-value -> mayor peso)
    Expression_Weight = 1 - (adj.P.Val)
  ) %>%
  select(Protein, Expression_Weight, Fold_Change, Is_Differential)

# Guardar en el directorio que espera el script Python
out_dir <- "VSR_Complete_Analysis/data/transcriptomics"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
out_file <- file.path(out_dir, "differential_expression.tsv")

write.table(final_results, file=out_file, sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("✅ Análisis finalizado con éxito. Archivo guardado en: %s\n", out_file))
