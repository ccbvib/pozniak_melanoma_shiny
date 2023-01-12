## This script converts the Seurat objects from the Pozniak, 2022 publication to
## pseudobulk datasets that are used to display in the Shiny apps found in this
## github repository

## load libaries
library("dplyr")
library("tidyr")
library("tibble")
library("Seurat")

## set server path to data files
fp <- paste0("/media/seq-srv-05/vrc/Project/Project_Theo/pan_cancer_shiny/")

## setlocal path to files
# fp <- paste0("~/Documents/tmp/2022/202210_joanna/counts/")

############################# Malignant_cells ################################## 
# read Seurat object
seu <- readRDS(file = paste0(fp, "Malignant_cells_v2.rds"))
# seu <- readRDS(file = paste0(fp, "data/Malignant_cells.rds"))
# gl <- readRDS(file = paste0(fp, "common_genes.rds"))

## extract Seurat metadata
seu@meta.data %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(cell_id = names(.)[1]) %>%
  as.data.frame() -> seu_meta

## pre-filtering step to reduce the size of the assay
cts <- GetAssayData(seu, assay = "RNA")
# cts <- cts[(which(rownames(cts) %in% gl)), ]
cts <- cts[(which(rownames(cts))), ]
seu <- subset(seu, features = rownames(cts))

## extract Seurat raw cts
raw_cts <- as.data.frame(seu@assays$RNA@counts)

## extract cell metadata for each unique cell label within a list
samp_list <- list()
for (i in 1:length(sort(unique(seu_meta$Malignant_clusters)))) {
seu_meta %>%
  dplyr::filter(Malignant_clusters == sort(unique(seu_meta$Malignant_clusters))[i]) %>%
  as.data.frame() -> samp_list[[i]]
}

## aggregate list of unique cells in order to rearrange and rename columns
dplyr::bind_rows(samp_list) %>%
  dplyr::mutate(stem_name = paste0(
    Malignant_clusters, "_",
    Response, "_", ## new category
    `BT/OT`, "_", `orig.ident`, "_", `GC number`, "_",
    row_number())) %>%
  as.data.frame() -> cell_list

## rearrange and rename columns of raw cts
raw_cts %>%
  dplyr::select(cell_list$cell_id) %>%
  tibble::rownames_to_column() %>%
  as.data.frame() -> raw_cts_sorted

## give new column names
names(raw_cts_sorted) <- c("gene", cell_list$stem_name)

## aggregate cts in long format
raw_cts_sorted %>%
  tidyr::gather(key = "cell_name", value = "expression", -gene) %>%
  dplyr::mutate(cell_name_st = gsub(" ", "", gsub("_[^_]+$", "", cell_name))) %>%
  dplyr::arrange(gene) %>%
  dplyr::group_by(cell_name_st, gene) %>%
  dplyr::summarise(mean_exp = mean(expression),
                   cell_count = n()) %>%
  dplyr::arrange(desc(cell_count)) %>%
  dplyr::mutate(cell_sample = paste0(cell_name_st, "_", cell_count)) %>%
  ungroup() %>%
  dplyr::select(-c(cell_count, cell_name_st)) %>%
  as.data.frame() -> agg_cts

## create a vector of unique cell type samples
uniq_samples <- sort(unique(agg_cts$cell_sample))

## create count columns for each cell_sample
cell_samp_list <- list()
for (i in 1:length(uniq_samples)) {
agg_cts %>%
  ungroup() %>%
  dplyr::filter(cell_sample == uniq_samples[i]) %>%
  tidyr::spread(key = cell_sample, value = mean_exp, fill = NA) %>%
  tibble::column_to_rownames(var = "gene") %>%
  as.data.frame() -> cell_samp_list[[i]]
}

## create wide count matrix
dplyr::bind_cols(sort(unique(agg_cts$gene)), cell_samp_list) %>%
  dplyr::rename(gene = names(.)[1]) %>%
  as.data.frame() -> wide_agg_cts

## create wide count matrix without pseudogenes or viral transcripts
wide_agg_cts %>%
  dplyr::filter(stringr::str_count(gene, "[.]") < 1,
                !grepl("-", gene)) -> wide_agg_cts_clean

## save aggregated count matrix
# saveRDS(object = wide_agg_cts, file = paste0(fp, "mal_agg_counts.rds"))
saveRDS(object = wide_agg_cts, file = paste0(fp, "mal_agg_counts_v2.rds"))

# fp <- "~/Documents/tmp/2022/202210_joanna2/pozniak_melanoma_shiny/"
# mc <- readRDS(file = paste0(fp, "counts/mal_agg_counts_v2.rds"))
# names(mc) <- gsub("Antigen_Presentation", "AntigenPresentation", names(mc))
# names(mc) <- gsub("Interferon_Alpha_Beta_Response", "InterferonAlphaBetaResponse", names(mc))
# names(mc) <- gsub("Mesenchymal_like", "Mesenchymallike", names(mc))
# names(mc) <- gsub("Mitochondrial\\(low_quality\\)", "Mitochondriallowquality", names(mc))
# names(mc) <- gsub("Neural_Crest_like", "NeuralCrestlike", names(mc))
# names(mc) <- gsub("Patient_specific_", "Patientspecific", names(mc))
# names(mc) <- gsub("Stress\\(HypoxiaResponse\\)", "StressHypoxiaResponse", names(mc))
# names(mc) <- gsub("Stress\\(p53Response\\)", "Stressp53Response", names(mc))
# saveRDS(object = mc, file = paste0(fp, "counts/mal_agg_counts_fixed_v2.rds"))

############################# NRAS13_all_43k ################################### 
seu <- readRDS(file = paste0(fp, "data/NRAS13_all_43k.rds"))

## extract Seurat metadata
seu@meta.data %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(cell_id = names(.)[1]) %>%
  dplyr::select(sampleID, cell_type, cell_id) %>%
  as.data.frame() -> seu_meta

## extract Seurat raw cts
raw_cts <- as.data.frame(seu@assays$SCT@counts)

## extract cell metadata for each unique cell label within a list
samp_list <- list()
for (i in 1:length(sort(unique(seu_meta$cell_type)))) {
seu_meta %>%
  dplyr::filter(cell_type == sort(unique(seu_meta$cell_type))[i]) %>%
  as.data.frame() -> samp_list[[i]]
}

## aggregate list of unique cells in order to rearrange and rename columns
dplyr::bind_rows(samp_list) %>%
  dplyr::mutate(stem_name = paste0(
    cell_type, "_", sampleID, "_", row_number())) %>%
  as.data.frame() -> cell_list

## rearrange and rename columns of raw cts
raw_cts %>%
  dplyr::select(cell_list$cell_id) %>%
  tibble::rownames_to_column() %>%
  as.data.frame() -> raw_cts_sorted

## give new column names
names(raw_cts_sorted) <- c("gene", cell_list$stem_name)

## aggregate cts in long format
raw_cts_sorted %>%
  tidyr::gather(key = "cell_name", value = "expression", -gene) %>%
  dplyr::mutate(cell_name_st = gsub(" ", "", gsub("_[^_]+$", "", cell_name))) %>%
  dplyr::arrange(gene) %>%
  dplyr::group_by(cell_name_st, gene) %>%
  dplyr::summarise(mean_exp = mean(expression),
                   cell_count = n()) %>%
  dplyr::arrange(desc(cell_count)) %>%
  dplyr::mutate(cell_sample = paste0(cell_name_st, "_", cell_count)) %>%
  ungroup() %>%
  dplyr::select(-c(cell_count, cell_name_st)) %>%
  as.data.frame() -> agg_cts

## create a vector of unique cell type samples
uniq_samples <- sort(unique(agg_cts$cell_sample))

## create count columns for each cell_sample
cell_samp_list <- list()
for (i in 1:length(uniq_samples)) {
agg_cts %>%
  ungroup() %>%
  dplyr::filter(cell_sample == uniq_samples[i]) %>%
  tidyr::spread(key = cell_sample, value = mean_exp, fill = NA) %>%
  tibble::column_to_rownames(var = "gene") %>%
  as.data.frame() -> cell_samp_list[[i]]
}

## create wide count matrix
dplyr::bind_cols(sort(unique(agg_cts$gene)), cell_samp_list) %>%
  dplyr::rename(gene = names(.)[1]) %>%
  as.data.frame() -> wide_agg_cts

## create wide count matrix without pseudogenes or viral transcripts
wide_agg_cts %>%
  dplyr::filter(stringr::str_count(gene, "[.]") < 1,
                !grepl("-", gene)) -> wide_agg_cts_clean

## save aggregated count matrix
saveRDS(object = wide_agg_cts, file = paste0(fp, "NRAS13_all_43k_agg_counts.rds"))

############################# NRAS13_malign_preprint ###########################
seu <- readRDS(file = paste0(fp, "data/NRAS13_malign_preprint.rds"))

## extract Seurat metadata
seu@meta.data %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(cell_id = names(.)[1]) %>%
  dplyr::select(orig.ident, seurat_clusters, cell_id) %>%
  dplyr::mutate(cluster = paste0("cluster_", seurat_clusters)) %>%
  as.data.frame() -> seu_meta

## extract Seurat raw cts
raw_cts <- as.data.frame(seu@assays$SCT@counts)

## extract cell metadata for each unique cell label within a list
samp_list <- list()
for (i in 1:length(sort(unique(seu_meta$cluster)))) {
seu_meta %>%
  dplyr::filter(cluster == sort(unique(seu_meta$cluster))[i]) %>%
  as.data.frame() -> samp_list[[i]]
}

## aggregate list of unique cells in order to rearrange and rename columns
dplyr::bind_rows(samp_list) %>%
  dplyr::mutate(stem_name = paste0(
    orig.ident, "_", cluster, "_", row_number())) %>%
  as.data.frame() -> cell_list

## rearrange and rename columns of raw cts
raw_cts %>%
  dplyr::select(cell_list$cell_id) %>%
  tibble::rownames_to_column() %>%
  as.data.frame() -> raw_cts_sorted

## give new column names
names(raw_cts_sorted) <- c("gene", cell_list$stem_name)

## aggregate cts in long format
raw_cts_sorted %>%
  tidyr::gather(key = "cell_name", value = "expression", -gene) %>%
  dplyr::mutate(cell_name_st = gsub(" ", "", gsub("_[^_]+$", "", cell_name))) %>%
  dplyr::arrange(gene) %>%
  dplyr::group_by(cell_name_st, gene) %>%
  dplyr::summarise(mean_exp = mean(expression),
                   cell_count = n()) %>%
  dplyr::arrange(desc(cell_count)) %>%
  dplyr::mutate(cell_sample = paste0(cell_name_st, "_", cell_count)) %>%
  ungroup() %>%
  dplyr::select(-c(cell_count, cell_name_st)) %>%
  as.data.frame() -> agg_cts

## create a vector of unique cell type samples
uniq_samples <- sort(unique(agg_cts$cell_sample))

## create count columns for each cell_sample
cell_samp_list <- list()
for (i in 1:length(uniq_samples)) {
agg_cts %>%
  ungroup() %>%
  dplyr::filter(cell_sample == uniq_samples[i]) %>%
  tidyr::spread(key = cell_sample, value = mean_exp, fill = NA) %>%
  tibble::column_to_rownames(var = "gene") %>%
  as.data.frame() -> cell_samp_list[[i]]
}

## create wide count matrix
dplyr::bind_cols(sort(unique(agg_cts$gene)), cell_samp_list) %>%
  dplyr::rename(gene = names(.)[1]) %>%
  as.data.frame() -> wide_agg_cts

## create wide count matrix without pseudogenes or viral transcripts
wide_agg_cts %>%
  dplyr::filter(stringr::str_count(gene, "[.]") < 1,
                !grepl("-", gene)) -> wide_agg_cts_clean

## save aggregated count matrix
saveRDS(object = wide_agg_cts, file = paste0(fp, "NRAS13_malign_preprint_agg_counts.rds"))