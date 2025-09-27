library(KEGGREST)
library(magrittr)
library(stringr)
library(stringr.plus)
library(dplyr)
library(tidyr)
library(qs)

pathways <- keggList("pathway", "hsa")
classification <- sub(" - Homo sapiens \\(human\\)", "", pathways)
brite <- keggGet("br:br08901")

strings <- str_split(brite, "\\nC {1,}")[[1]]

strings <- str_remove(strings, "\\+C\tMap number\n!") %>% 
    str_remove("\n!\n#\n#Last updated: June 3, 2025\n#&raquo; Japanese version")

pathway_info <- data.frame(
    order1 = str_extract(strings, "(\nA[:alnum:]{1,})( ?[:alnum:])*"),
    order2 = str_extract(strings, "\nB( +.{1,})+"),
    ids = str_extract(strings, "^[:digit:]{1,}"), 
    name = str_remove(strings, "^[:digit:]{1,} +") %>% str_remove_all("\nA.+|\nB.+")
) %>% 
    dplyr::mutate(order1 = str_remove(order1, "\nA"), 
        order2 = str_remove(order2, "\nB"))

pathway_info <- pathway_info %>% 
    tidyr::fill(order1, order2, .direction = "down")
pathway_info %>% head

pathway_info <- pathway_info[-1, ]
qsave(pathway_info, "../../datasets/kegg_pathways_information_hsa.qs")