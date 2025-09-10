# perform polarization analysis

PolarizationGeneSetEnrichmentHyperTest = function(degs, 
    input_celltype, 
    species = "mouse", 
    ref_dir = NULL) {
  `%notin%` = Negate(`%in%`)
  library(openxlsx)
  library(plyr)
  library(dplyr)
  set.seed(101)
  
  celltypes = c("ILC", readLines("data/lig_seurat_listcells.txt"))
  if (input_celltype %notin% celltypes) stop(paste0("cell type must be one of the following: ", paste0(celltypes, collapse = ", ")))
  
  # celltype_spreadsheet = read.xlsx("dataFiles/celltype_list.xlsx")
  # input_celltype_display = mapvalues(input_celltype, from = celltype_spreadsheet$Celltype_OriginalName,
  # to = celltype_spreadsheet$Celltype_DisplayName, warn_missing = FALSE)
  
  print ("Hyper")
  print (input_celltype_display)
  
  # Load the pre-computed significantly differentially expressed genes for each cytokine in each cell type
  ref_deg_sig_celltype <- read.xlsx("data/polarization.xlsx", sheet = 1)
  ref_deg_sig_celltype <- ref_deg_sig_celltype %>% dplyr::filter(CellType == input_celltype)

  if(is.null(ref_dir)){
    stop("please give the referece dataset directory.")
  }
  
  ref_files <- dir_ls(ref_dir, type = "file")
  file <- ref_files[grepl(celltypes, ref_files)]
  if(length(file) < 1){
    stop("please enter a valide cell type.")
  }

  ref_expressed_genes = readRDS(file)
  # ref_expressed_genes_celltype = subset(ref_expressed_genes, celltype == input_celltype)
  
  # subset into the celltype of interest
  polarization_states = sort(unique(ref_deg_sig_celltype$Polarization))
  
  # Perform fisher's exact test
  res_pvals = c()
  res_ess = c()
  
  for (ss in polarization_states) {
    markers_ss = subset(ref_deg_sig_celltype, Polarization == ss)
    marker_genes <- markers_ss$Gene %>% str_split(",") %>% unlist
    
    num_overlap = length(intersect(degs, marker_genes))
    num_only_user = length(degs) - num_overlap
    num_only_ref = nrow(markers_ss) - num_overlap
    num_neither = length(setdiff(rownames(ref_expressed_genes),
                                 c(marker_genes, degs)))
    
    mat_test = matrix(c(num_overlap, num_only_user, num_only_ref, num_neither), nrow = 2)
    
    test_pval = fisher.test(mat_test)$p.value
    test_es = num_overlap / length(degs)
    
    res_pvals = c(res_pvals, test_pval)
    res_ess = c(res_ess, test_es)
    
  }
  
  # Construct a result matrix
  df_irea = data.frame(Polarization = polarization_states, ES = res_ess, pval = res_pvals)
  
  # Perform multiple hypothesis testing correction
  df_irea$padj = p.adjust(df_irea$pval, method = "fdr")
  
  # Add pseudocount for log transform
  df_irea$nlog10_padj = pmax(0, -log10(df_irea$padj+10e-50))
  df_irea$Celltype = input_celltype
  
  return(df_irea)
}