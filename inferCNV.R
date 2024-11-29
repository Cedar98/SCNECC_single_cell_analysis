# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="singleCell.counts.matrix",
                                    annotations_file="cellAnnotations.txt",
                                    delim="\t",
                                    gene_order_file="gene_ordering_file.txt",
                                    ref_group_names=c("normal"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T
)