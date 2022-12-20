
dfmarker <- data.frame(marker_method_name = c("meanratio2", "wilcoxon", "ttest"), 
                       method_type1 = c("pairwise_comparisons", 
                                        "pairwise_comparisons",
                                        "pairwise_comparisons"), 
                       method_type2 = c("ratios", "differences", "differences"),
                       data_scale = c("lognorm_counts", "lognorm_counts", "lognorm_counts"),
                       package.name = c("DeconvoBuddies", "scran", "scran"), 
                       function_name = c("get_mean_ratios2", "findMarkers", "findMarkers"))
save(dfmarker, file = "lute_marker-method_transfer_learning.rda")