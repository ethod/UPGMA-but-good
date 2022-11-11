source('upgma_source.R')

library(ape)
library(Biostrings)

require("GEOquery")

#initialization

gds <- getGEO('GDS3850')

GDS_dist_mat <- init_gds(gds)

image(x=1:nrow(GDS_dist_mat), y=1:nrow(GDS_dist_mat), GDS_dist_mat)

tmp <- read.table("mat_01.txt", stringsAsFactors=FALSE, header=TRUE)

tmp_UPGMA <- upgma(tmp)

GDS_UPGMA <- upgma(GDS_dist_mat)


#plot(nj(GDS_dist_mat)) for reference

GDS_family.mat <- matrix(c(unlist(GDS_UPGMA[1])), ncol = 3 )
GDS_desc <- matrix (NA, ncol = ncol(GDS_family.mat)-1, nrow = nrow(GDS_family.mat))

node_nb <- unlist(GDS_UPGMA[2])

for (i in c(1:3)){
  if (i <=2){
    GDS_desc[,i] <- (GDS_family.mat[,i])
  }else{
    GDS_dist <- c((GDS_family.mat[,i]))
  }
}

#normalize dist? 

for (i in c(1:length(GDS_dist))){
  GDS_dist[i] <- (GDS_dist[i]-min(GDS_dist))/(max(GDS_dist)-min(GDS_dist))
}

descendance <- list()

class(descendance) <- "phylo"

descendance$edge <-GDS_desc

descendance$Nnode <- node_nb

descendance$node.label <- c(1:node_nb)

descendance$tip.label <- c(letters, LETTERS, c('aa', 'bb'))

descendance$edge.length <- GDS_dist

#plot.phylo(descendance, show.node.label = TRUE)
#i would advise against running this line it makes my R session crash 

##########