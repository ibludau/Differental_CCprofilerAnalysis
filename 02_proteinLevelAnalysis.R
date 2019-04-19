#' # Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rds")
design_matrix <- readRDS("design_matrix.rds")
calibration <- readRDS("calibration.rds")

#' # Integrate peptide traces across conditions
pepTracesIntegrated_sum <- integrateTraceIntensities(pepTracesList_filtered,
                                                     design_matrix = NULL,
                                                     integrate_within = NULL,
                                                     aggr_fun = "sum")

# pepTracesIntegrated_sum <- subset(pepTracesIntegrated_sum, unique(pepTracesIntegrated_sum$trace_annotation$protein_id)[1:30], trace_subset_type = "protein_id")

#' # Perform protein-centric analysis 
# Suppress messages for clarity
suppressMessages(
  proteinFeatures <-  findProteinFeatures(traces = pepTracesIntegrated_sum,
                                          corr_cutoff = 0.9,
                                          window_size = 7,
                                          parallelized = F,
                                          n_cores = 1,
                                          collapse_method = "apex_only",
                                          perturb_cutoff = "5%",
                                          rt_height = 2,
                                          smoothing_length = 7)
)
saveRDS(proteinFeatures,"proteinFeatures.rda")

#' Score protein features and filter for a detection FDR of 5%
proteinFeaturesFiltered <- scoreFeatures(proteinFeatures, FDR=0.05, PDF=TRUE, name="qvalueStats_proteinFeatures")
saveRDS(proteinFeaturesFiltered,"proteinFeaturesFiltered.rda")

#' Evaluate protein feature stats 
summarizeFeatures(proteinFeaturesFiltered,
                  plot=TRUE,
                  PDF=TRUE,
                  name="feature_summary_proteinFeaturesFiltered")

#' ## Extract feature values
proteinFeaturesFiltered_featureVals <- extractFeatureVals(traces = pepTracesList_filtered,
                                             features = proteinFeaturesFiltered,
                                             design_matrix = design_matrix,
                                             extract = "subunits_detected",
                                             imputeZero = T,
                                             verbose = F)
saveRDS(proteinFeaturesFiltered_featureVals, "proteinFeaturesFiltered_featureVals.rda")

#' ## Fill feature values
proteinFeaturesFiltered_featureValsFilled <- fillFeatureVals(featureVals = proteinFeaturesFiltered_featureVals,
                                                tracesList = pepTracesList_filtered,
                                                design_matrix = design_matrix)
saveRDS(proteinFeaturesFiltered_featureValsFilled, "proteinFeaturesFiltered_featureValsFilled.rda")

#' ## Perform peptide-level differential expression testing for all features
proteinFeaturesFiltered_DiffExprPep <- testDifferentialExpression(featureVals = proteinFeaturesFiltered_featureValsFilled,
                                                     compare_between = "Condition",
                                                     level = "peptide",
                                                     measuredOnly = FALSE)
saveRDS(proteinFeaturesFiltered_DiffExprPep, "proteinFeaturesFiltered_DiffExprPep.rda")

pdf("proteinFeaturesFiltered_DiffExprPep_stats.pdf", width = 4, height=4)
  hist(proteinFeaturesFiltered_DiffExprPep$pVal, breaks = 100)
  hist(proteinFeaturesFiltered_DiffExprPep$global_pVal, breaks = 100)
  hist(log(proteinFeaturesFiltered_DiffExprPep$qint1), breaks = 100)
  hist(log(proteinFeaturesFiltered_DiffExprPep$qint2), breaks = 100)
  hist(log(proteinFeaturesFiltered_DiffExprPep$global_int2_imp), breaks = 100)
  hist(log(proteinFeaturesFiltered_DiffExprPep$global_int1_imp), breaks = 100)
  hist(proteinFeaturesFiltered_DiffExprPep$log2FC, breaks = 50)
  hist(proteinFeaturesFiltered_DiffExprPep$global_log2FC, breaks = 50)
dev.off()

#' ## Aggregate differential expression results to the protein level
proteinFeaturesFiltered_DiffExprProtein <- aggregatePeptideTests(proteinFeaturesFiltered_DiffExprPep)
saveRDS(proteinFeaturesFiltered_DiffExprProtein, "proteinFeaturesFiltered_DiffExprProtein.rda")

pdf("proteinFeaturesFiltered_DiffExprProtein_stats.pdf", width = 4, height=4)
  hist(proteinFeaturesFiltered_DiffExprProtein$pVal, breaks = 100)
  hist(proteinFeaturesFiltered_DiffExprProtein$global_pVal, breaks = 100)
  hist(proteinFeaturesFiltered_DiffExprProtein$medianLog2FC, breaks = 50)
  hist(proteinFeaturesFiltered_DiffExprProtein$global_medianLog2FC, breaks = 50)
dev.off()

#' # Make volcano plots
#' General volcano plots
plotVolcano(proteinFeaturesFiltered_DiffExprPep, PDF = T, pBHadj_cutoff = 0.05, name = "proteinFeaturesFiltered_DiffExprPep")
plotVolcano(proteinFeaturesFiltered_DiffExprProtein, PDF = T, pBHadj_cutoff = 0.05, name = "proteinFeaturesFiltered_DiffExprProtein")

plotVolcano(proteinFeaturesFiltered_DiffExprPep, PDF = T, pBHadj_cutoff = 0.05, name = "proteinFeaturesFiltered_DiffExprPep", level = "global")
plotVolcano(proteinFeaturesFiltered_DiffExprProtein, PDF = T, pBHadj_cutoff = 0.05, name = "proteinFeaturesFiltered_DiffExprProtein", level = "global")

#' Get differential proteins
allProteins <- unique(proteinFeaturesFiltered_DiffExprProtein$feature_id)
diffProteins <- unique(proteinFeaturesFiltered_DiffExprProtein[pBHadj<0.05][abs(medianLog2FC)>1]$feature_id)
write.table(allProteins, "allProteins.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(diffProteins, "diffProteins.txt", col.names = F, row.names = F, quote = F, sep = "\t")
diffProteins_global <- unique(proteinFeaturesFiltered_DiffExprProtein[global_pBHadj<0.05][abs(global_medianLog2FC)>1]$feature_id)
write.table(diffProteins_global, "diffProteins_global.txt", col.names = F, row.names = F, quote = F, sep = "\t")

#' Volcanoplts highlighting specific proteins
library(ggrepel)
highlight_prot <- diffProteins[1]
plotVolcano(proteinFeaturesFiltered_DiffExprProtein, highlight=highlight_prot, PDF = T, pBHadj_cutoff = 0.05, name = "proteinFeaturesFiltered_DiffExprProtein_highlight")
plotVolcano(proteinFeaturesFiltered_DiffExprProtein, highlight=highlight_prot, PDF = T, pBHadj_cutoff = 0.05, name = "proteinFeaturesFiltered_DiffExprProtein_highlight", level = "global")

#' # Plot some differential protein features
pdf("diffProteinTraces.pdf", width=10,height=7)
for (id in diffProteins) {
  plotFeatures(feature_table = proteinFeaturesFiltered,
               traces = pepTracesList_filtered,
               feature_id = id,
               design_matrix=design_matrix,
               calibration=calibration,
               peak_area = T,
               legend = F,
               onlyBest = F, 
               monomer_MW=T)
}
dev.off()

pdf("diffProteinTraces_global.pdf", width=10,height=7)
for (id in diffProteins_global) {
  plotFeatures(feature_table = proteinFeaturesFiltered,
               traces = pepTracesList_filtered,
               feature_id = id,
               design_matrix=design_matrix,
               calibration=calibration,
               peak_area = T,
               legend = F,
               onlyBest = F, 
               monomer_MW=T)
}
dev.off()

##################################
##################################
##################################

testFilled_norm <- proteinFeaturesFiltered_featureValsFilled[feature_id %in% diffProteins]

library(betareg)
library(lmtest)

testFilled_norm <- testFilled_norm

local_vs_global_test <- testLocalVsGlobal(testFilled_norm)
saveRDS(local_vs_global_test,"local_vs_global_test.rds")

local_vs_global_test_protein <- aggregateLocalVsGlobal(local_vs_global_test)
saveRDS(local_vs_global_test_protein,"local_vs_global_test_protein.rds")

local_vs_global_test_sig <- unique(local_vs_global_test_protein[pBHadj<0.05][abs(medianDiff) > 0.3]$feature_id)
write.table(local_vs_global_test_sig, "local_vs_global_test_sig.txt", col.names = F, row.names = F, quote = F, sep = "\t")

pdf("local_vs_global_test_sig_features.pdf", width=10,height=7)
for (id in local_vs_global_test_sig) {
  plotFeatures(feature_table = proteinFeaturesFiltered,
               traces = pepTracesList_filtered,
               feature_id = id,
               design_matrix=design_matrix,
               calibration=calibration,
               peak_area = T,
               legend = F,
               onlyBest = F, 
               monomer_MW=T)
}
dev.off()

##################################
##################################
##################################
