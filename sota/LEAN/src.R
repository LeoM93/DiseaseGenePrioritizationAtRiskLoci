library(LEANR)

args <- commandArgs(trailingOnly = TRUE)
res<-run.lean(ranking=args[1], network=args[2],
add.scored.genes=T, verbose=T, n_reps=1000, ncores=3)
write.csv(res$restab,args[3], row.names = TRUE)
