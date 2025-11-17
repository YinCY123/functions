import scvi
import cell2location
from cell2location.models import RegressionModel, Cell2location
from cell2location.utils.filtering import filter_genes

# Load your scRNA reference (AnnData)
adata_ref = scvi.data.read_h5ad("scRNA_reference.h5ad")

# Train reference model
RegressionModel.setup_anndata(adata=adata_ref, labels_key='cell_type')
mod = RegressionModel(adata_ref)
mod.train(max_epochs=250)

# Export reference signatures
adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000})
inf_aver = mod.get_reference_signatures(adata_ref)

# Load spatial data
adata_vis = scvi.data.read_h5ad("spatial_data.h5ad")

# Train cell2location model
Cell2location.setup_anndata(adata=adata_vis, batch_key="batch")
mod2 = Cell2location(adata_vis, cell_state_df=inf_aver)
mod2.train(max_epochs=30000)

# Save results
adata_vis = mod2.export_posterior(adata_vis, sample_kwargs={'num_samples': 1000})
adata_vis.write("cell2loc_results.h5ad")