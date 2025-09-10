# perform cytoking response analysis

GeneSetEnrichmentHyperTest = function(degs, input_celltype, species = "mouse") {
  `%notin%` = Negate(`%in%`)
  library(openxlsx)
  library(plyr)
  library(dplyr)
  set.seed(0)
  
  celltypes = c("ILC", readLines("sourceFiles/lig_seurat_data.txt"))
  if (input_celltype %notin% celltypes) stop(paste0("cell type must be one of the following: ", paste0(celltypes, collapse = ", ")))
  
  # Load the pre-computed significantly differentially expressed genes for each cytokine in each cell type
  ref_deg_sig_celltype = read.xlsx("dataFiles/SuppTable3_CytokineSignatures.xlsx", sheet = input_celltype)
  ref_expressed_genes = readRDS("dataFiles/ref_expressed_genes_per_celltype.Rda")
  ref_expressed_genes_celltype = subset(ref_expressed_genes, celltype == input_celltype)
  
  # subset into the celltype of interest
  samples = setdiff(read.xlsx("dataFiles/irea_cytokine_list.xlsx")$Cytokine_OriginalName,
                    c("IL2+IL15", "IL2+IFNg"))
  
  # Perform fisher's exact test
  res_pvals = c()
  res_ess = c()
  
  for (ss in samples) {
    markers_ss = subset(ref_deg_sig_celltype, Cytokine_Str == ss)
    
    num_overlap = length(intersect(degs, markers_ss$Gene))
    num_only_user = length(degs) - num_overlap
    num_only_ref = nrow(markers_ss) - num_overlap
    num_neither = length(setdiff(ref_expressed_genes_celltype$genes_expressed,
                                 c(markers_ss$Gene, degs)))
    
    mat_test = matrix(c(num_overlap, num_only_user, num_only_ref, num_neither), nrow = 2)
    
    test_pval = fisher.test(mat_test)$p.value
    test_es = num_overlap / length(degs)
    
    res_pvals = c(res_pvals, test_pval)
    res_ess = c(res_ess, test_es)
    
  }
  
  # Construct a result matrix
  df_irea = data.frame(Cytokine = samples, ES = res_ess, pval = res_pvals)
  
  # Perform multiple hypothesis testing correction
  df_irea$padj = p.adjust(df_irea$pval, method = "fdr")
  
  # Add pseudocount for log transform
  df_irea$nlog10_padj = pmax(0, -log10(df_irea$padj+10e-50))
  df_irea$Celltype = input_celltype
  
  # Order the results by p-value
  df_irea = df_irea[order(df_irea$pval), ]
  
  # Change to official cytokine name
  cytokine_spreadsheet = read.xlsx("dataFiles/irea_cytokine_list.xlsx")
  df_irea$Cytokine = mapvalues(df_irea$Cytokine,
                                from = cytokine_spreadsheet$Cytokine_OriginalName,
                                to = cytokine_spreadsheet$Cytokine_DisplayName, warn_missing = FALSE)
  
  return(df_irea)
}
