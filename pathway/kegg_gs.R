gene_sets <- keggLink(target = "pathway", source = "hsa") 
gene_sets <- data.frame(pathway = gene_sets, 
                        entrez = names(gene_sets))
gene_sets <- gene_sets %>% 
  dplyr::mutate(pathway = str_replace(pathway, "path:", "") %>% str_trim(side = "both"), 
                entrez = str_replace(entrez, "hsa:", "") %>% str_trim(side = "both"))
gene_sets$symbol <- mapIds(org.Hs.eg.db, 
                           keys = gene_sets$entrez, keytype = "ENTREZID", column = "SYMBOL")

hsa_pathways <- keggList("pathway", "hsa")
gene_sets$pathway_name <- hsa_pathways[gene_sets$pathway]
gene_sets$pathway_name <- str_replace(gene_sets$pathway_name, " - Homo sapiens \\(human\\)", "") %>% str_trim(side = "both")
