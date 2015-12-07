# Created by Suthinan Rujirapipat, 7 December 2015
# As part of master project @ University of Sunderland, Information Technology Management

library(igraph, lib=".\\library")
library(ProNet, lib=".\\library")

# http://string-db.org/api/psi-mi-tab/interactionsList?identifiers=APP&limit=10000&required_score=900

APP <- read.table(file = ".\\data\\APPPSIMI.csv", header = FALSE)
APP <- APP[,c("V3","V4")]
APP[,1] = toupper(APP[,1])
APP[,2] = toupper(APP[,2])
APP <- unique(APP)

g <- graph.data.frame(APP, directed=FALSE) # For MCODE
g2 <- graph.data.frame(APP, directed=FALSE) # For MCL

setwd(".\\output")

#======================= MCODE generation and visualisation =======================

# MCODE network clustering
# graph	<- An igraph object.
# vwp	  <-Vertex weight percentage. Default value is 0.5.
# haircut <-	Boolean value, whether to remove singly-connected nodes from clusters (TRUE) or not (FALSE).
# fluff	<- Boolean value, whether to spand cluster cores by one neighbour shell outwards (TRUE) or not (FALSE).
# fdt	  <- Cluster density cutoff. Default value is 0.8.
# loops	<-  Boolean value, whether to include self-loops (TRUE) or not (FALSE).

mcg <- mcode(g, loops = TRUE, vwp=0.2, haircut= TRUE, fluff = TRUE, fdt= 0.2)

# Create and visualise of induced subgraphs (clusters), just in case
# 
# cluster1 <- induced.subgraph(g, mcg$COMPLEX[[1]])
# writeClipboard(V(cluster1)$name)
# 
# visualization(cluster1,
#               node.size=4,
#               node.label=V(cluster1)$name,
#               node.label.color="blue",
#               edge.color="gray", 
#               edge.width=0.1,
#               node.fill.color = "red")
# 
summary(mcg$COMPLEX)

index <- which(!is.na(mcg$score))
membership <- rep(0, vcount(g))

for (i in 1:length(index)) {
  membership[mcg$COMPLEX[[index[i]]]] <- i
}
# names(membership) <- V(g)$name

# Find top 5 largest clusters; This will ignore community overlapping nodes using first in first out assignment
sort(table(unlist(membership)), decreasing = TRUE)[1:5]

color <- "white" # Initialised
color[membership == 1] <- "red"
color[membership == 2] <- "pink"
color[membership == 4] <- "aquamarine"
color[membership == 5] <- "yellow"
color[membership == 9] <- "magenta"


# Create network graph of the result
png(file=paste(format(Sys.time(),"%H-%M-%S"),"MCODE", "png", sep = "."),height=1500, width=1500, bg="white")
visualization(graph = g, 
                layout="kamada.kawai",
                node.size=4,
                node.label.color="blue",
                edge.color = "gray",
                edge.width=0.1,
                node.fill.color = color)
              # node.fill.color = list(heat.colors(5)))
title(main = "Molecular Complex Detection Algorithm (MCODE), VWP = 0.2", cex.main = 4)
legend("bottomright", legend = c('color[membership == 1] <- "red"', 'color[membership == 2] <- "pink"',
                                 'color[membership == 4] <- "aquamarine"', 'color[membership == 5] <- "yellow"',
                                 'color[membership == 9] <- "magenta"'), pch = 1, title = "Cluster Size")
dev.off()

#======================= MCL generation and visualisation ==================================

# Convert igraph graph object to adjacency matrix

adj <- matrix(rep(0, length(V(g2))^2), nrow = length(V(g2)), ncol = length(V(g2)))
for (i in 1:length(V(g2))) {
  neighbors <- neighbors(g2, v = V(g2)$name[i], mode = "all")
  j <- match(neighbors$name, V(g2)$name, nomatch = 0)
  adj[i, j] = 1
}

# Clear temporary variables
rm(i)
rm(j)
rm(neighbors)

# Markov Cluster Algorithm
# x	<- an adjacency or (n x n) matrix
# addLoops	<-logical; if TRUE, self-loops with weight 1 are added to each vertex of x (see Details).
# expansion	<- numeric value > 1 for the expansion parameter
# inflation	<- numeric value > 0 for the inflation power coefficient
# allow1	<- logical; if TRUE, vertices are allowed to form their own cluster. 
# If FALSE, clusters of size 1 are interpreted as background noise and grouped in one cluster.
# max.iter	<- an interger, the maximum number of iterations for the MCL
# ESM	<- logical whether the equilibrium state matrix should be returned (default is FALSE)

lc <- mcl(adj, addLoops = TRUE, inflation = 1.8, allow1 = TRUE, ESM = FALSE)

# Check MCL complexes' size
lc

# Find top 5 largest clusters
sort(table(unlist(lc$Cluster)), decreasing = TRUE)[1:5]

# lc$Cluster <- lc$Cluster

# Assign color according to the top 5
color2 <- "white" # Initialised
color2[lc$Cluster == 3] <- "red"
color2[lc$Cluster == 1] <- "pink"
color2[lc$Cluster == 5] <- "aquamarine"
color2[lc$Cluster == 6] <- "yellow"
color2[lc$Cluster == 7] <- "magenta"

# Create network graph of the result
png(file=paste(format(Sys.time(),"%H-%M-%S"),"MCL", "png", sep = "."),height=1500, width=1500, bg="white")
visualization(graph = g2, 
                layout="kamada.kawai",
                node.size=4,
                node.label.color="blue",
                edge.color = "gray",
                edge.width=0.1,
                node.fill.color = color2)
title(main = "Markov Clustering (MCL), Top 5 largest clusters", cex.main = 4)
legend("bottomright", legend = c('color[membership == 3] <- "red"', 'color[membership == 1] <- "pink"',
                                 'color[membership == 5] <- "aquamarine"', 'color[membership == 6] <- "yellow"',
                                 'color[membership == 7] <- "magenta"'), pch = 1, title = "Cluster Size")
dev.off()

# generateMCL function, used To find the induced subgraph 
# generateMCL <- function (x, colorx) {
#   
#   mclCluster <- induced.subgraph(g2, lc$Cluster == x)
#   
#   visualization(graph = mclCluster, 
#                 layout="fruchterman.reingold",
#                 node.size=8,
#                 node.label=V(mclCluster)$name,
#                 node.label.color="blue",
#                 edge.color = "gray",
#                 edge.width=0.1,
#                 node.fill.color = colorx)
# }
# 
# generateMCL(2,"yellow")