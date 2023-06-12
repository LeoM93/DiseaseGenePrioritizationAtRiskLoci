source("/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/sota/dmGWAS/R/dms.R")
source("/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/sota/dmGWAS/R/dms_2.4.R")
args <- commandArgs(trailingOnly = TRUE)

network <- read.csv(args[1], sep='\t', header=FALSE)
gene.weight <- read.csv(args[2], sep='\t', header=FALSE)
res.list = dms_2.4(network, gene.weight)

print(res.list)
write.table(res.list[["zi.ordered"]], file=args[3], quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
