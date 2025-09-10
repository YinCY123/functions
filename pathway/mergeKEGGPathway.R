library(tkWidgets)
library(KEGGgraph)
library(KEGG.db)

pathway_ids <- c("00010", "00061", "00020")

pathways = list.files(pattern="*.xml")
for (i in 1:length(pathways)) assign(pathways[i], parseKGML2Graph(pathways[i],expandGenes=TRUE))
rm(i);rm(pathways)
pathways<-objNameToList(objects(), parent.frame())
# merged <- ugraph(mergeKEGGgraphs(pathways, edgemode = "directed")) # undirected graph
merged <- mergeKEGGgraphs(pathways, edgemode = "directed") # directed graph