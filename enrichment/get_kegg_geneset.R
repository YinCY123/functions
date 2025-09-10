get_kegg_geneset <- function(species, remove_suffix = TRUE, ...){
  # loading packages
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(KEGGREST))
  ifelse(species == "hsa", suppressPackageStartupMessages(require(org.Hs.eg.db)), 
    ifelse(species == "mmu", suppressPackageStartupMessages(require(org.Mm.eg.db)), species))

  gene_pathway <- keggLink(target = "pathway", source = species)
  pathway_name <- keggList(database = "pathway", organism = species)

  pathway_df <- data.frame(entrez = str_remove(names(gene_pathway), paste0("^", species, ":")), 
      pathway_id = str_remove(gene_pathway, pattern = "^path:")) %>% 
    dplyr::mutate(name = pathway_name[pathway_id])

  if(species == "hsa"){
    pathway_df <- pathway_df %>% 
      dplyr::mutate(symbol = mapIds(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID", column = "SYMBOL"))
  }else if(species == "mmu"){
    pathway_df <- pathway_df %>% 
      dplyr::mutate(symbol = mapIds(org.Mm.eg.db, keys = entrez, keytype = "ENTREZID", column = "SYMBOL"))
  }

  if(remove_suffix){
    pattern = ifelse(species == "hsa", " - Homo sapiens \\(human\\)", 
      ifelse(species == "mmu", " - Mus musculus \\(house mouse\\)", NA))
    pathway_df <- pathway_df %>% 
      dplyr::mutate(name = str_remove(name, pattern = pattern))
  }
  return(pathway_df)
}