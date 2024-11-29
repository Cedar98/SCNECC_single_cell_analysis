source("http://bioconductor.org/biocLite.R")
biocLite()

biocLite("monocle")

library(monocle)

install.packages("devtools")
devtools::install_github("cole-trapnell-lab/monocle-release@develop")

biocLite(c("DDRTree", "pheatmap"))

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd)

cth <- newCellTypeHierarchy()

MYF5_id <- row.names(subset(fData(cds), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(cds), gene_short_name == "ANPEP"))

cth <- addCellType(cth, "Myoblast", classify_func =
                     function(x) { x[MYF5_id,] >= 1 })
cth <- addCellType(cth, "Fibroblast", classify_func =
                     function(x) { x[MYF5_id,] < 1 & x[ANPEP_id,] > 1 } )

cds <- classifyCells(cds, cth, 0.1)

cds <- clusterCells(cds)

disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds)
cds <- orderCells(cds)

diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Media")
sig_genes <- subset(diff_test_res, qval < 0.1)

#do not run
HSMM_expr_matrix <- read.table("fpkm_matrix.txt")
HSMM_sample_sheet <- read.delim("cell_sample_sheet.txt")
HSMM_gene_annotation <- read.delim("gene_annotations.txt")

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                       phenoData = pd, featureData = fd)

# Where 'data_to_be_imported' can either be a Seurat object
# or an SCESet.

importCDS(data_to_be_imported)

# We can set the parameter 'import_all' to TRUE if we'd like to
# import all the slots from our Seurat object or SCESet.
# (Default is FALSE or only keep minimal dataset)

importCDS(data_to_be_imported, import_all = TRUE)

lung <- load_lung()

# To convert to Seurat object
lung_seurat <- exportCDS(lung, 'Seurat')

# To convert to SCESet
lung_SCESet <- exportCDS(lung, 'Scater')

#Do not run
HSMM <- newCellDataSet(count_matrix,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily=negbinomial.size())

HSMM <- newCellDataSet(as(umi_matrix, "sparseMatrix"),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

cellranger_pipestance_path <- "/path/to/your/pipeline/output/directory"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)

fd <- fData(gbm)

# The number 2 is picked arbitrarily in the line below.
# Where "2" is placed you should place the column number that corresponds to your
# featureData's gene short names.

colnames(fd)[2] <- "gene_short_name"

gbm_cds <- newCellDataSet(exprs(gbm),
                          phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                          featureData = new("AnnotatedDataFrame", data = fd),
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)

# First create a CellDataSet from the relative expression levels
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.1,
                       expressionFamily = tobit(Lower = 0.1))

# Next, use it to estimate RNA counts
rpc_matrix <- relative2abs(HSMM, method = "num_genes")

# Now, make a new CellDataSet using the RNA counts
HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))

HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

print(head(pData(HSMM)))

valid_cells <- row.names(subset(pData(HSMM),
                                Cells.in.Well == 1 &
                                  Control == FALSE &
                                  Clump == FALSE &
                                  Debris == FALSE &
                                  Mapped.Fragments > 1000000))
HSMM <- HSMM[,valid_cells]

pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(HSMM), color = Hours, geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
               pData(HSMM)$Total_mRNAs < upper_bound]
HSMM <- detectGenes(HSMM, min_expr = 0.1)

# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(HSMM),
                             gene_short_name == "ANPEP"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Myoblast", classify_func =
                     function(x) { x[MYF5_id,] >= 1 })
cth <- addCellType(cth, "Fibroblast", classify_func = function(x)
{ x[MYF5_id,] < 1 & x[ANPEP_id,] > 1 })

HSMM <- classifyCells(HSMM, cth, 0.1)

table(pData(HSMM)$CellType)

pie <- ggplot(pData(HSMM),
              aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)

# HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM, return_all = F) # norm_method='log'

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "CellType",
                   markers = c("MYF5", "ANPEP"))

plot_cell_clusters(HSMM, 1, 2, color = "Media")

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~Media + num_genes_expressed",
                        verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "CellType")

HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "Cluster") +
  facet_wrap(~CellType)

marker_diff <- markerDiffTable(HSMM[expressed_genes,],
                               cth,
                               residualModelFormulaStr = "~Media + num_genes_expressed",
                               cores = 1)

candidate_clustering_genes <-
  row.names(subset(marker_diff, qval < 0.01))
marker_spec <-
  calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3))

semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
plot_ordering_genes(HSMM)

plot_pc_variance_explained(HSMM, return_all = F)

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 3,
                        norm_method = 'log',
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~Media + num_genes_expressed",
                        verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "CellType")

HSMM <- clusterCells(HSMM,
                     num_clusters = 2,
                     frequency_thresh = 0.1,
                     cell_type_hierarchy = cth)
plot_cell_clusters(HSMM, 1, 2, color = "CellType",
                   markers = c("MYF5", "ANPEP"))

pie <- ggplot(pData(HSMM), aes(x = factor(1), fill =
                                 factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                      fullModelFormulaStr = "~Media")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
plot_ordering_genes(HSMM_myo)

HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,
                            method = 'DDRTree')

HSMM_myo <- orderCells(HSMM_myo)

plot_cell_trajectory(HSMM_myo, color_by = "Hours")

plot_cell_trajectory(HSMM_myo, color_by = "State")

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
HSMM_myo <- orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")

plot_cell_trajectory(HSMM_myo, color_by = "State") +
  facet_wrap(~State, nrow = 1)

blast_genes <- row.names(subset(fData(HSMM_myo),
                                gene_short_name %in% c("CCNB2", "MYOD1", "MYOG")))
plot_genes_jitter(HSMM_myo[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)

HSMM_expressed_genes <-  row.names(subset(fData(HSMM_myo),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "Hours")

HSMM_myo <- detectGenes(HSMM_myo, min_expr = 0.1)
fData(HSMM_myo)$use_for_ordering <-
  fData(HSMM_myo)$num_cells_expressed > 0.05 * ncol(HSMM_myo)

plot_pc_variance_explained(HSMM_myo, return_all = F)

HSMM_myo <- reduceDimension(HSMM_myo,
                            max_components = 2,
                            norm_method = 'log',
                            num_dim = 3,
                            reduction_method = 'tSNE',
                            verbose = T)

HSMM_myo <- clusterCells(HSMM_myo, verbose = F)

plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)')

plot_rho_delta(HSMM_myo, rho_threshold = 2, delta_threshold = 4 )

HSMM_myo <- clusterCells(HSMM_myo,
                         rho_threshold = 2,
                         delta_threshold = 4,
                         skip_rho_sigma = T,
                         verbose = F)

plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)')

clustering_DEG_genes <-
  differentialGeneTest(HSMM_myo[HSMM_expressed_genes,],
                       fullModelFormulaStr = '~Cluster',
                       cores = 1)

HSMM_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

HSMM_myo <-
  setOrderingFilter(HSMM_myo,
                    ordering_genes = HSMM_ordering_genes)

HSMM_myo <-
  reduceDimension(HSMM_myo, method = 'DDRTree')

HSMM_myo <-
  orderCells(HSMM_myo)

HSMM_myo <-
  orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))

plot_cell_trajectory(HSMM_myo, color_by = "Hours")

disp_table <- dispersionTable(HSMM_myo)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.5 &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id

CCNB2_id <-
  row.names(subset(fData(HSMM_myo), gene_short_name == "CCNB2"))
MYH3_id <-
  row.names(subset(fData(HSMM_myo), gene_short_name == "MYH3"))

cth <- newCellTypeHierarchy()

cth <- addCellType(cth,
                   "Cycling myoblast",
                   classify_func = function(x) { x[CCNB2_id,] >= 1 })

cth <- addCellType(cth,
                   "Myotube",
                   classify_func = function(x) { x[MYH3_id,] >= 1 })

cth <- addCellType(cth,
                   "Reserve cell",
                   classify_func =
                     function(x) { x[MYH3_id,] == 0 & x[CCNB2_id,] == 0 })

HSMM_myo <- classifyCells(HSMM_myo, cth)

marker_diff <- markerDiffTable(HSMM_myo[HSMM_expressed_genes,],
                               cth,
                               cores = 1)
#semisup_clustering_genes <-
#row.names(subset(marker_diff, qval < 0.05))
semisup_clustering_genes <-
  row.names(marker_diff)[order(marker_diff$qval)][1:1000]

HSMM_myo <- setOrderingFilter(HSMM_myo, semisup_clustering_genes)
#plot_ordering_genes(HSMM_myo)
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,
                            method = 'DDRTree', norm_method = 'log')
HSMM_myo <- orderCells(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by = "CellType") +
  theme(legend.position = "right")

HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]

my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))

cds_subset <- HSMM_filtered[my_genes,]
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "Hours",
                               ncol = 1)

marker_genes <- row.names(subset(fData(HSMM_myo),
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5",
                                                        "ANPEP", "PDGFRA","MYOG",
                                                        "TPM1",  "TPM2",  "MYH2",
                                                        "MYH3",  "NCAM1", "TNNT1",
                                                        "TNNT2", "TNNC1", "CDK1",
                                                        "CDK2",  "CCNB1", "CCNB2",
                                                        "CCND1", "CCNA1", "ID1")))

diff_test_res <- differentialGeneTest(HSMM_myo[marker_genes,],
                                      fullModelFormulaStr = "~Media")

# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)

sig_genes[,c("gene_short_name", "pval", "qval")]

MYOG_ID1 <- HSMM_myo[row.names(subset(fData(HSMM_myo),
                                      gene_short_name %in% c("MYOG", "CCNB2"))),]
plot_genes_jitter(MYOG_ID1, grouping = "Media", ncol= 2)

to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("UBC", "NCAM1", "ANPEP")))
cds_subset <- HSMM[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~CellType")
diff_test_res[,c("gene_short_name", "pval", "qval")]

plot_genes_jitter(cds_subset,
                  grouping = "CellType",
                  color_by = "CellType",
                  nrow= 1,
                  ncol = NULL,
                  plot_trend = TRUE)

full_model_fits <-
  fitModel(cds_subset,  modelFormulaStr = "~CellType")
reduced_model_fits <- fitModel(cds_subset, modelFormulaStr = "~1")
diff_test_res <- compareModels(full_model_fits, reduced_model_fits)
diff_test_res

to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("MYH3", "MEF2C", "CCNB2", "TNNT1")))
cds_subset <- HSMM_myo[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res[,c("gene_short_name", "pval", "qval")]

plot_genes_in_pseudotime(cds_subset, color_by = "Hours")

diff_test_res <- differentialGeneTest(HSMM_myo[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)

to_be_tested <-
  row.names(subset(fData(HSMM),
                   gene_short_name %in% c("TPM1", "MYH3", "CCNB2", "GAPDH")))

cds_subset <- HSMM[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~CellType + Hours",
                                      reducedModelFormulaStr = "~Hours")
diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_jitter(cds_subset,
                  grouping = "Hours", color_by = "CellType", plot_trend = TRUE) +
  facet_wrap( ~ feature_label, scales= "free_y")

lung <- load_lung()
plot_cell_trajectory(lung, color_by = "Time")

BEAM_res <- BEAM(lung, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

plot_genes_branched_heatmap(lung[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

lung_genes <- row.names(subset(fData(lung),
                               gene_short_name %in% c("Ccnd2", "Sftpb", "Pdpn")))
plot_genes_branched_pseudotime(lung[lung_genes,],
                               branch_point = 1,
                               color_by = "Time",
                               ncol = 1)

