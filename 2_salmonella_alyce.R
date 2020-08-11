

setwd("/Users/gh11/alyce_salmonella/")

## summarise the results from calculating the pairwise distance between every two genomes


## do it for one then move on to the other
dists = read.table("summary_specific_paraB.tab", sep = ",", comment.char = "", stringsAsFactors = F,
                        header = T)

final_distance_per_gene = rep("",dim(dists)[1])
for (i in 1:dim(dists)[1]) {
  curr = table(t(dists[i,-c(1,2)]))
}




