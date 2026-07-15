# geo_timecourse_analysis.R
# Análisis de series de tiempo para VSR (GSE196587)
# Contrastes: 4H vs 0H, 8H vs 0H, 12H vs 0H, 15H vs 0H

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="http://cran.us.r-project.org")
if (!require("GEOquery")) BiocManager::install("GEOquery", update=FALSE)
if (!require("limma")) BiocManager::install("limma", update=FALSE)
if (!require("dplyr")) install.packages("dplyr", repos="http://cran.us.r-project.org")

library(GEOquery)
library(limma)
library(dplyr)

# 1. Descargar dataset GSE196587
cat("Descargando dataset de time-course GSE196587 desde NCBI GEO...\n")
gse <- getGEO("GSE196587", GSEMatrix=TRUE, AnnotGPL=TRUE)
if (length(gse) > 1) {
  idx <- grep("GPL570", attr(gse, "names"))
  if (length(idx) == 0) idx <- 1
  gse <- gse[[idx]]
} else {
  gse <- gse[[1]]
}

# 2. Extraer datos
pdata <- pData(gse)
eset <- exprs(gse)
fdata <- fData(gse)

# 3. Identificar tiempos
# Los títulos son del tipo "HRSV at 0H, biological rep2"
conditions <- as.character(pdata$title)
timepoints <- rep("Unknown", length(conditions))
timepoints[grepl("0H", conditions)] <- "H0"
timepoints[grepl("4H", conditions)] <- "H4"
timepoints[grepl("8H", conditions)] <- "H8"
timepoints[grepl("12H", conditions)] <- "H12"
timepoints[grepl("15H", conditions)] <- "H15"

# Asegurarse de que tenemos un factor ordenado
time_factor <- factor(timepoints, levels=c("H0", "H4", "H8", "H12", "H15"))

# 4. Diseño experimental
cat("Procesando matriz de diseño experimental (Time-Course)...\n")
design <- model.matrix(~ 0 + time_factor)
colnames(design) <- levels(time_factor)

fit <- lmFit(eset, design)

# 5. Contrastes (Cada tiempo contra H0)
contrast.matrix <- makeContrasts(
  T4_vs_0 = H4 - H0,
  T8_vs_0 = H8 - H0,
  T12_vs_0 = H12 - H0,
  T15_vs_0 = H15 - H0,
  levels=design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 6. Extraer resultados para todos los contrastes
cat("Extrayendo resultados por cada punto de tiempo...\n")

get_results_for_contrast <- function(fit2, coef_name, pdata_symbol_col="Gene symbol") {
  res <- topTable(fit2, coef=coef_name, number=Inf, sort.by="none")
  if (pdata_symbol_col %in% colnames(fdata)) {
    res$Protein <- fdata[rownames(res), pdata_symbol_col]
  } else {
    res$Protein <- rownames(res)
  }
  
  # Limpiar símbolos múltiples (ej. "AKT1///AKT2" -> "AKT1")
  res$Protein <- sapply(strsplit(as.character(res$Protein), "///"), `[`, 1)
  res$Protein <- trimws(res$Protein)
  
  # Seleccionar lo necesario
  df <- res %>% 
    filter(!is.na(Protein) & Protein != "" & Protein != "---") %>%
    group_by(Protein) %>%
    arrange(adj.P.Val) %>%
    slice(1) %>%
    ungroup()
  
  return(df)
}

res_4h <- get_results_for_contrast(fit2, "T4_vs_0")
res_8h <- get_results_for_contrast(fit2, "T8_vs_0")
res_12h <- get_results_for_contrast(fit2, "T12_vs_0")
res_15h <- get_results_for_contrast(fit2, "T15_vs_0")

# 7. Unir todo en una sola matriz de tiempo
final_results <- res_4h %>% select(Protein, logFC_4H = logFC, pval_4H = adj.P.Val) %>%
  inner_join(res_8h %>% select(Protein, logFC_8H = logFC, pval_8H = adj.P.Val), by="Protein") %>%
  inner_join(res_12h %>% select(Protein, logFC_12H = logFC, pval_12H = adj.P.Val), by="Protein") %>%
  inner_join(res_15h %>% select(Protein, logFC_15H = logFC, pval_15H = adj.P.Val), by="Protein")

# Promediar el peso para compatibilidad con red estática si se requiere
final_results <- final_results %>%
  mutate(Expression_Weight = 1 - ((pval_4H + pval_8H + pval_12H + pval_15H) / 4))

out_dir <- "VSR_Complete_Analysis/data/transcriptomics"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
out_file <- file.path(out_dir, "timecourse_expression.tsv")

write.table(final_results, file=out_file, sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("✅ Time-course finalizado con éxito. Archivo guardado en: %s\n", out_file))
