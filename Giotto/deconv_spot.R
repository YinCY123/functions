deconv_spot <- function(sces,
                        spes,
                        sc_cell_id_column,
                        sc_celltype_column,
                        sc_gene_id_column,
                        st_spot_ID_column = "cell_ID",
                        st_gene_ID_column = "feat_ID",
                        project = "ddls",
                        threads = 8,
                        ...){
  # loading required packages
  library(qs2)
  library(SpatialDDLS)
  library(SpatialExperiment)
  
  # create spatialDDLS object
  ddls <- createSpatialDDLSobject(
    sc.data = sces, 
    sc.cell.ID.column = sc_cell_id_column, 
    sc.cell.type.column = sc_celltype_column, 
    sc.gene.ID.column = sc_gene_id_column, 
    st.data = spes, 
    st.spot.ID.column = st_spot_ID_column, 
    st.gene.ID.column = st_gene_ID_column, 
    project = project
  )
  
  ddls <- genMixedCellProp(
    object = ddls, 
    cell.ID.column = sc_cell_id_column, 
    cell.type.column = sc_celltype_column, 
    num.sim.spots = 10000, 
    n.cells = 50
  )
  
  ddls <- simMixedProfiles(
    object = ddls, 
    type.data = "both", 
    threads = threads
  )
  
  ddls <- trainDeconvModel(
    object = ddls, 
    type.data.train = "mixed", 
    threads = threads
  )
  
  ddls <- deconvSpatialDDLS(
    object = ddls, 
    index.st = 1, 
    k.spot = 10, 
    fast.pca = TRUE
  )
  
  return(ddls)
}
