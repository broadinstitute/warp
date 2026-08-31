# 1.1.2
2026-08-21 (Date of Last Commit)

* Fixed the confusion-matrix heatmaps (conf_Ttype_*.png, conf_ConsType_*.png) in 04_clusterability.py. The classification pickles store conf_mat as raw integer counts (0 to 860, rows summing to 164-546) but the heatmap colour scale is fixed to [0, 1], so every non-zero cell saturated to the darkest blue: a single misassigned cell rendered identically to a correct diagonal entry of 860. The figures were pure black-and-white noise and carried no information. Rows are now normalised to fractions before plotting, matching notebooks/3_evaluation.ipynb, which divides each row by its sum. On the K=89 run this turns 550 uniformly-dark cells into a diagonal averaging 0.94 against an off-diagonal averaging 0.0007.
* Widened the accuracy bar chart (classAcc_RF_K_*.png). At the previous width the three group labels ran together as "t-typesT CategoriesT Categories".
* Category annotations in 05_state_traversal.py are no longer ambiguous. _label_for_cat truncated the dominant t-type at its first space, so "L6 IT Car3" and "L6 CT Nxph2" both became "L6" and a ten-panel figure could be labelled L6, Pvalb, L6, L6, L4, Vip, L6, L6, Pvalb, Pvalb with no way to tell the categories apart. Labels now carry the category number, e.g. "c34 L6".
* State-scatter highlight colours no longer collide with the background. _default_palette used tab10, whose 8th entry is grey (#7f7f7f) against a #dbdbde background, so the 8th selected category's panel looked empty. Switched to husl, which is evenly spaced around the hue circle and returns no neutrals.
* All four are figure-rendering changes only. No model, no numeric output, and no upstream mmidas change: the MMIDAS package stays pinned at warp-v2. Re-running MMIDAS_Analyze is sufficient; MMIDAS_DataPrep and MMIDAS_Train outputs are unaffected.

# 1.1.1
2026-08-11 (Date of Last Commit)

* Documented the two unavoidable divergences from the reference notebooks in dashboard.md: 5_state_traversal.ipynb hardcodes a hand-picked selected_c list that a generic workflow cannot reproduce, and the reference train/test split is unseeded and therefore unrecoverable.
* Updated the Docker image to us.gcr.io/broad-gotc-prod/mmidas:1.0.0-0.1.0-1787578739 (MMIDAS pinned at warp-v2, which reverts the min_con pruning stop to upstream behaviour).

# 1.1.0
2026-08-06 (Date of Last Commit)

* Fixed KEGG pathway mapping, which silently produced zero pathways whenever a kegg_toml was supplied. mmidas.utils.data_tools.load_data read gene identifiers from adata.var.values (the column block, empty for an .h5ad written with no var columns) instead of the var index, so 03c_traversal_prep.py received an empty gene list, mapped no pathways, and exited 0. All pathway box plots in stage 05 were therefore missing without any error. load_data now reads adata.var.index.values and warns when the gene-identifier count does not match the number of expression columns.
* 05_state_traversal.py now selects categories by assigned-cell count, largest first, instead of taking the first n_selected_cats by category index. Surviving pruning is not the same as having cells assigned: on a collapsed model the index-ordered selection landed on ten categories that were all empty, and the ten resulting figures were pixel-identical apart from their titles. Empty categories are now excluded with a warning, and the task exits non-zero if no category has any cells at all.
* state_traversal_manifest.json gains selected_c_n_cells, n_active_categories, and n_populated_categories so a validation step can distinguish a real traversal from figures drawn over empty categories. Note n_selected_cats is now capped at the number of populated categories, so it can be below the requested n_selected_cats input.
* Fixed the silhouette figure (SC_K_*.png) in 04_clusterability.py. The sc_T arrays are shape (1, n_category); the code sorted them as 1-D and took len(), which is 1, so every point was plotted at x=1 as a separate broadcast line and each inherited the same legend label. The result was a 1507x9406 px figure with 134 legend entries and an x-axis spanning 0.96-1.04. The arrays are now flattened before sorting.
* Accuracy bar chart x-axis labels are no longer rotated and overlapping when n_arm > 1.
* The Docker image now installs the mmidas package from a pinned revision of the fork rather than a copy of a local directory, so any published image can be rebuilt from source. See "Docker image provenance" in dashboard.md. The load_data and train fixes above are candidates for upstreaming to AllenInstitute/MMIDAS.
* Updated the Docker image to us.gcr.io/broad-gotc-prod/mmidas:1.0.0-0.1.0-1786046379

# 1.0.0
2026-06-24 (Date of Last Commit)

* Initial release. Consumes the confirmed evaluation_results.json from MMIDAS_Train and runs RF classification + silhouette analysis (03b) and state traversal preparation (03c) in parallel, followed by clusterability figures (04) and state traversal figures (05).
