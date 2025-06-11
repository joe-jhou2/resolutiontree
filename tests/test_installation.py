import scanpy as sc
import src as rt

print("Loading PBMC dataset...")
adata = sc.datasets.pbmc3k()


print("Running preprocessing...")
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

print("Testing cluster_resolution_finder...")
resolutions = [0.0, 0.2, 0.5, 1.0]
resolutions = [0.0, 0.2, 0.5, 1.0, 1.5, 2.0]

adata_new = rt.cluster_resolution_finder(adata,
                                         resolutions=resolutions,
                                         n_top_genes=3,
                                         min_cells=2,
                                         deg_mode="within_parent",
                                         inplace=False
                                         )

rt.cluster_decision_tree(adata_new, resolutions=resolutions, 
                        output_settings = {
                            "output_path": "examples/test_pypi_pbmc3k.png",
                            "draw": False,
                            "figsize": (12, 6),
                            "dpi": 300
                            },
                        node_style = {
                            "node_size": 500,
                            "node_colormap": None,
                            "node_label_fontsize": 12
                            },
                        edge_style = {
                            "edge_color": "parent",
                            "edge_curvature": 0.01,
                            "edge_threshold": 0.01,
                            "show_weight": True,
                            "edge_label_threshold": 0.05,
                            "edge_label_position": 0.8,
                            "edge_label_fontsize": 8
                            },
                        gene_label_settings = {
                            "show_gene_labels": True,
                            "n_top_genes": 2,
                            "gene_label_threshold": 0.001,
                            "gene_label_style": {"offset":0.5, "fontsize":8},
                            },
                        level_label_style = {
                            "level_label_offset": 15,
                            "level_label_fontsize": 12
                            },
                        title_style = {
                            "title": "Hierarchical Leiden Clustering",
                            "title_fontsize": 20
                            },
                        layout_settings = {
                            "node_spacing": 5.0,
                            "level_spacing": 1.5
                            },
                        clustering_settings = {
                            "prefix": "leiden_res_",
                            "edge_threshold": 0.05
                            }
                    )
print("All tests passed! 🎉")