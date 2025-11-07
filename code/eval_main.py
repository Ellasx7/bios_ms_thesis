#Thesis Main Script#

# required packages
import warnings
warnings.filterwarnings(action='ignore')
import numpy as np
import scanpy as sc
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.preprocessing import StandardScaler
import squidpy as sq




#-------------------- Input Data --------------------#
# 1. feature expression matrix, e.g.: gene expression matrix, image feature matrix
# 2. spatial coordinates of pairwise alignment/registration samples - coordinates of query sample is mapped to ref sample
# 3. 

# Pairwise, B1-B2, B2-B3, and so on
# gene expression cost for VALIS from scratch (no need to calculate distance_df ahead)

# # parameter settings
# ref_id = 'B1'
# ref_type = 'ref'
# qry_id = 'B2'
# qry_type = 'query'
# total_count = 1e4
# top_n_hvg = 100
# feature_type = 'gene expression'

# # directory settings
# spatial_coords_dir = './Desktop/HER2_data/test2_b_img/warped_spot_coords'
# feat_dir = './Desktop/HER2_data'+'/'+ref_type+'/'+ref_id

# spatial_coords_dir = './sample_data/spatial_coords'
# feat_dir = './sample_data/features'


### Inputs Handling ###
# 1. Spatial information: pixal coordinates of the reference section & ref-based warped pixal coordiantes of the adjacent(query) section.
	# The ref file should have two columns, namely 'x' and 'y', which store the original spatial coordinates, and indexed by spot id.
	# The qry file should have two columns, namely 'x' and 'y', which store the warped spatial coordinates, and indexed by spot id.
# 2. Feature expression matrices: Reference feature expression matrix & Query feature expression matrix. 
	# The expression matrix should be a nxm matrix, where it has n spots and m features.
# ----------------------------------------------------------------------------------------------------

# ref_spatial = pd.read_csv(f'{ref_spatial_dir}/{ref_id}_selection.tsv', delimiter = '\t')
# ref_spatial[['x','y']] = ref_spatial[['x','y']].astype(str)
# ref_spatial.index = [f'{x}x{y}' for x, y in zip(ref_spatial['x'], ref_spatial['y'])]
# ref_spatial = ref_spatial.loc[ref_idx]
# ref_coords = ref_spatial[['pixel_x', 'pixel_y']].rename(columns={'pixel_x': 'x', 'pixel_y': 'y'})
# qry_spatial = pd.read_csv(f'{qry_spatial_dir}/{qry_id}_selection.tsv', delimiter = '\t')
# qry_spatial[['x','y']] = qry_spatial[['x','y']].astype(str)
# qry_spatial.index = [f'{x}x{y}' for x, y in zip(qry_spatial['x'], qry_spatial['y'])]
# qry_spatial = qry_spatial.loc[qry_idx]
# qry_coords = qry_spatial[['registered_x', 'registered_y']].rename(columns={'registered_x': 'x', 'registered_y': 'y'})

# sc.pp.highly_variable_genes(ref_exp_adata, flavor="seurat_v3", n_top_genes=top_n_svg)
# ref_exp_adata = ref_exp_adata[:, ref_exp_adata.var['highly_variable']]
# sc.pp.highly_variable_genes(qry_exp_adata, flavor="seurat_v3", n_top_genes=top_n_svg)
# qry_exp_adata = qry_exp_adata[:, qry_exp_adata.var['highly_variable']]

#====================================================================================================

def rank_svg(adata):
	adata.obsm["spatial"] = adata.obs[["x", "y"]].values # array_x and array_y used to search neighbors
	sq.gr.spatial_neighbors(adata)
	sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names)
	morans_i = adata.uns["moranI"]
	return morans_i

def intersection_with_order(list1, list2):
    # Use a set for faster look-up in list2
    set2 = set(list2)
    return [item for item in list1 if item in set2]

def compute_sum_of_squares(idx, features_df):
    feature_vector = features_df.loc[idx].to_numpy()
    return np.sum(feature_vector ** 2)


### 03/02/2025 latest version to achieve symmetric cost###

def compute_weighted_cost(ref_id, qry_id, spatial_coords_dir, feat_dir,
						  feature_type = ['gene_expression','image_features'],
						  top_n_svg = 100, 
						  d = 500,
						  total_count = 1e4, 
						  verbose = True):
	'''
	This function takes the original spatial coordinates of a reference section and the mapped spatial coordinates of a query section and the expression matrices of both and returns the expression cost of the pairwise alignment/registration of these sections. 
	
	Parameters
	----------
	ref_id (str): Reference section id in pairwise alignment.
	qry_id (str): The id of the section to be warped in pairwise alignment.
	spatial_coords_dir (str): Directory to ref and qry spatial data.
	feat_dir (str): Directory to feature expression matrices.
	feature_type (str): The type of feature we use to evaluate the alignment performance.
	top_n_svg (int): The number of highly variable genes to keep. Default to 100.
	d (int): The spatial Euclidean distance limit. Default to 500.
	total_count (int): The total count per spot/cell to normalize to. Default to 10,000.
	verbose (bool): Show detailed messages. Default to True.
	
	Returns
	----------
	float: The weighted per-spot cost.
	'''
	
	# read in reference feature and spatial coordinates #
	if verbose == True:
		print(f'Reading in reference section {ref_id}\'s feature expression matrix and spatial coordinates...')
	if feature_type == 'gene_expression':
		ref_exp = pd.read_csv(f'{feat_dir}/{ref_id}_{feature_type}.tsv', delimiter = '\t', index_col=0)
	else:
		ref_exp = pd.read_csv(f'{feat_dir}/{ref_id}_{feature_type}.csv', index_col=0)
	ref_feat_names = ref_exp.columns.tolist() # feature name list
	ref_idx = ref_exp.index.tolist() # spot index/id list
	ref_spatial = pd.read_csv(f'{spatial_coords_dir}/{ref_id}_to_{qry_id}_{feature_type}_ref_coords.csv', index_col=0)
	ref_spatial = ref_spatial.loc[ref_idx]

	# read in query feature and warped spatial coordinates #
	if verbose == True:
		print(f'Reading in query section {qry_id}\'s feature expression matrix and spatial coordinates...')
	if feature_type == 'gene_expression':
		qry_exp = pd.read_csv(f'{feat_dir}/{qry_id}_{feature_type}.tsv', delimiter = '\t', index_col=0)
	else:
		qry_exp = pd.read_csv(f'{feat_dir}/{qry_id}_{feature_type}.csv', index_col=0)
	qry_feat_names = qry_exp.columns.tolist()
	qry_idx = qry_exp.index.tolist()
	qry_spatial = pd.read_csv(f'{spatial_coords_dir}/{ref_id}_to_{qry_id}_{feature_type}_qry_coords.csv', index_col=0)
	qry_spatial = qry_spatial.loc[qry_idx]
	
	if verbose == True:
		print('Processing feature expression data...')
	# feature preprocess: --> select SVGs --> Scale
	if feature_type == 'gene_expression':
		ref_exp_adata = sc.AnnData(ref_exp)
		sc.pp.normalize_total(ref_exp_adata, target_sum=total_count) # normalize to a total count of 10000 per spot
		sc.pp.log1p(ref_exp_adata)
		sc.pp.scale(ref_exp_adata)
		ref_exp_adata.obs = ref_spatial.copy()
		ref_morans_i = rank_svg(ref_exp_adata)
		ref_svg_list = ref_morans_i.index.tolist()
		
		qry_exp_adata = sc.AnnData(qry_exp)
		sc.pp.normalize_total(qry_exp_adata, target_sum=total_count)
		sc.pp.log1p(qry_exp_adata)
		sc.pp.scale(qry_exp_adata)
		qry_exp_adata.obs = qry_spatial.copy()
		qry_morans_i = rank_svg(qry_exp_adata)
		qry_svg_list = qry_morans_i.index.tolist()
		intersection = intersection_with_order(ref_svg_list,qry_svg_list)
		top_svg_list = intersection[0:top_n_svg]
		ref_exp_adata = ref_exp_adata[:, top_svg_list]
		# sc.pp.scale(ref_exp_adata)
		ref_exp = pd.DataFrame(ref_exp_adata.X.copy())
		ref_exp.index = ref_idx
		ref_exp.columns = ref_exp_adata.var.index.tolist()
		qry_exp_adata = qry_exp_adata[:, top_svg_list]
		# sc.pp.scale(qry_exp_adata)
		qry_exp = pd.DataFrame(qry_exp_adata.X.copy())
		qry_exp.index = qry_idx
		qry_exp.columns = qry_exp_adata.var.index.tolist()
		ref_feat_names = top_svg_list
		qry_feat_names = top_svg_list
	else:
		z_scaler = StandardScaler()
		ref_exp = z_scaler.fit_transform(ref_exp)
		ref_exp = pd.DataFrame(ref_exp)
		ref_exp.index = ref_idx
		ref_exp.columns = ref_feat_names
		qry_exp = z_scaler.fit_transform(qry_exp)
		qry_exp = pd.DataFrame(qry_exp)
		qry_exp.index = qry_idx
		qry_exp.columns = qry_feat_names
	if verbose == True:
		print(f'-------------------- Number of features used: {len(ref_feat_names)} --------------------')

	# Compute spatial Euclidean distances
	spatial_distance = cdist(ref_spatial.to_numpy(), qry_spatial.to_numpy(), metric='euclidean')
	spatial_distance = pd.DataFrame(spatial_distance, index=ref_idx, columns=qry_idx)

	# Initial check on the alignment
	all_qry_neighbor = {i: spatial_distance[i].idxmin() for i in qry_idx}
	qry_neighbor_within_d = sum([spatial_distance[i].min() <= d for i in all_qry_neighbor])
	percentage_within_d = qry_neighbor_within_d / len(all_qry_neighbor)
	if percentage_within_d >= 0.95:
		match_criteria = True
	else:
		match_criteria = False
	if verbose == True:
		print(f"Initial assessment of the alignment: {'Good' if match_criteria == True else 'Not good'}")

	# Find spatially nearest neighbor
	ref_neighbor = {i: spatial_distance.loc[i].idxmin() for i in ref_idx if spatial_distance.loc[i].min() <= d}
	qry_neighbor = {i: spatial_distance[i].idxmin() for i in qry_idx if spatial_distance[i].min() <= d}
	ref_set = {(key, value) for key, value in ref_neighbor.items()}
	qry_set = {(value, key) for key, value in qry_neighbor.items()}
	overlapping_set = ref_set & qry_set
	if verbose == True:
		print(f"Maximum number of possible pairs: {min(len(ref_idx),len(qry_idx))}, actual pairs that are matched: {len(overlapping_set)}.")
	ref_unmatch = set(ref_idx) - {item[0] for item in overlapping_set}
	qry_unmatch = set(qry_idx) - {item[1] for item in overlapping_set}

	# Compute feature Euclidean distances and expression cost
	feature_distance = cdist(ref_exp.to_numpy(), qry_exp.to_numpy(), metric='euclidean')
	feature_distance = pd.DataFrame(feature_distance, index=ref_idx,columns=qry_idx)
	overlapping_df = pd.DataFrame(list(overlapping_set), columns=["ref_index", "qry_index"]).sort_values(by = "ref_index")
	matched_cost = np.sum([feature_distance.loc[item[0], item[1]] for _, item in overlapping_df.iterrows()])
	unmatched_cost = 0
	if match_criteria:
		# If match criteria is met, use nearest neighbor for unmatched spots
		all_ref_neighbor = {i: spatial_distance.loc[i].idxmin() for i in ref_idx}
		for ref in ref_unmatch:
			unmatched_cost += feature_distance.loc[ref, all_ref_neighbor[ref]]
		for qry in qry_unmatch:
			unmatched_cost += feature_distance.loc[all_qry_neighbor[qry],qry]
	else:
		# Otherwise, use sum of squares for unmatched spots(add penalty)
		for ref in ref_unmatch:
			unmatched_cost += compute_sum_of_squares(ref, ref_exp)
		for qry in qry_unmatch:
			unmatched_cost += compute_sum_of_squares(qry, qry_exp)
	total_spots = len(ref_idx)+len(qry_idx)
	averaged_cost = (1/total_spots)*(2*matched_cost+unmatched_cost)
	if verbose == True:
		print(f'The calculated {feature_type} cost for {ref_id}-{qry_id} alignment is {averaged_cost:.2f}.')
		print()
	return round(averaged_cost,2)

#====================================================================================================
