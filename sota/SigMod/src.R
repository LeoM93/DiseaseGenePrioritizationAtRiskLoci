require( igraph )
if(packageVersion("igraph") < "1.0.0") {
    stop("Need to install igraph version 1.0.0 or above")
}

R_dir = '/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/sota/SigMod/R' ## the path to the SigMod_v2\R folder
## load all R functions in this folder
file.sources = list.files(R_dir, pattern='*.R$',full.names=TRUE, ignore.case=TRUE)
tmp=sapply(file.sources, source, .GlobalEnv)

args <- commandArgs(trailingOnly = TRUE)

network_file = args[1]
network_data = read.table(network_file,header = TRUE, stringsAsFactors = FALSE,sep='\t', quote='"')
############## Read gene p-values data ##############
gene_ps_file = args[2]

gene_ps = read.table(gene_ps_file, stringsAsFactors = FALSE, header =TRUE, sep='\t',)
colnames(gene_ps) = c("gene", "p")
############## Construct a scored network ##############
## assign correct value to interaction_indices
interaction_indices = c(1,2)
## assign correct value to weight_index
weight_index = NULL
## construct a scored-network
scored_net = construct_scored_net(network_data,
interaction_indices= interaction_indices, gene_ps=gene_ps)
res_info = SigMod_bisection( net=scored_net)
data_frame <- rbind(V(res_info$opt_module$selected)$name, V(res_info$opt_module$selected_next)$name)
write.table(V(res_info$opt_module$selected_next)$name, file=args[3], quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
