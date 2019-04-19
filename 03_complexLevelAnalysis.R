#' # Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rds")
protein_traces_list <- readRDS("protein_traces_list.rds")
design_matrix <- readRDS("design_matrix.rds")
calibration <- readRDS("calibration.rds")

#' ## Summarize protein traces across conditions
protein_traces_sum <- integrateTraceIntensities(protein_traces_list,
                                                design_matrix = NULL,
                                                integrate_within = NULL,
                                                aggr_fun = "sum")

#' ## Save protein traces
saveRDS(protein_traces_sum,"protein_traces_sum.rds")

#corumHypotheses <- readRDS("../../../html/SECpaper/output_data/corum/complexTargetsPlusDecoys.rds")
#saveRDS(corumHypotheses,"corumHypotheses.rds")
corumHypotheses <- readRDS("corumHypotheses.rds")

#' ## Complex feature finding
complexFeatures <-  findComplexFeatures(traces = protein_traces_sum,
                                        complex_hypothesis = corumHypotheses,
                                        corr_cutoff = 0.9,
                                        window_size = 7,
                                        parallelized = F,
                                        n_cores = 1,
                                        collapse_method = "apex_network",
                                        perturb_cutoff = "5%",
                                        rt_height = 2,
                                        smoothing_length = 7)

saveRDS(complexFeatures,"corumComplexFeatures.rda")
# complexFeatures <- readRDS("corumComplexFeatures.rda")

complexFeatures <- filterFeatures(complexFeatures,
                                  complex_ids = NULL,
                                  protein_ids = NULL,
                                  min_feature_completeness = NULL,
                                  min_hypothesis_completeness = NULL,
                                  min_subunits = NULL,
                                  min_peak_corr = NULL,
                                  min_monomer_distance_factor = 2)

hypothesis <- "corum"

plotSummarizedMScoverage(hypotheses = corumHypotheses, protein_traces_sum, PDF = T)

complexFeaturesBest <- getBestFeatures(complexFeatures)
saveRDS(complexFeaturesBest,"complexFeaturesBest.rda")

complexFeaturesBestFiltered <- scoreFeatures(complexFeaturesBest, FDR=0.05, PDF=T, name="qvalueStats_complexFeatures")

scoredDataAll <- appendSecondaryComplexFeatures(scoredPrimaryFeatures = complexFeaturesBestFiltered, allFeatures = complexFeatures, peakCorr_cutoff = 0.7)
saveRDS(scoredDataAll,"scoredDataAll.rds")

plotSummarizedComplexes(scoredDataAll, corumHypotheses, protein_traces_sum, PDF=T, name="complex_completeness_pie")

summarizeFeatures(scoredDataAll,
                  plot=TRUE,
                  PDF=T,
                  name="feature_summary_scoredDataAll")

summarizeFeatures(complexFeaturesBestFiltered,
                  plot=TRUE,
                  PDF=T,
                  name="feature_summary_complexFeaturesBestFiltered")

#' ## Extract feature values
complex_featureVals <- extractFeatureVals(traces = pepTracesList_filtered,
                                          features = scoredDataAll,
                                          design_matrix = design_matrix,
                                          extract = "subunits_detected",
                                          imputeZero = T,
                                          verbose = F)
saveRDS(complex_featureVals, "complex_featureVals.rda")

#' ## Fill feature values 
complex_featureValsFilled <- fillFeatureVals(featureVals = complex_featureVals,
                                             tracesList = pepTracesList_filtered,
                                             design_matrix = design_matrix)
saveRDS(complex_featureValsFilled, "complex_featureValsFilled.rda")

#' ## Perform peptide-level differential expression testing for all features 
complex_DiffExprPep <- testDifferentialExpression(featureVals = complex_featureValsFilled,
                                                  compare_between = "Condition",
                                                  level = "peptide",
                                                  measuredOnly = FALSE)

saveRDS(complex_DiffExprPep, "complex_DiffExprPep.rda")
# complex_DiffExprPep <- readRDS("complex_DiffExprPep.rda")

pdf("complex_DiffExprPep_stats.pdf", width = 4, height=4)
  hist(complex_DiffExprPep$pVal, breaks = 100)
  hist(complex_DiffExprPep$global_pVal, breaks = 100)
  hist(log(complex_DiffExprPep$qint1), breaks = 100)
  hist(log(complex_DiffExprPep$qint2), breaks = 100)
  hist(log(complex_DiffExprPep$global_int2_imp), breaks = 100)
  hist(log(complex_DiffExprPep$global_int1_imp), breaks = 100)
  hist(complex_DiffExprPep$log2FC, breaks = 50)
  hist(complex_DiffExprPep$global_log2FC, breaks = 50)
dev.off()

#' ## Aggregate differential expression results to the protein level
complex_DiffExprProtein <- aggregatePeptideTests(complex_DiffExprPep)
saveRDS(complex_DiffExprProtein, "complex_DiffExprProtein.rda")

pdf("complex_DiffExprProtein_stats.pdf", width = 4, height=4)
  hist(complex_DiffExprProtein$pVal, breaks = 100)
  hist(complex_DiffExprProtein$global_pVal, breaks = 100)
  hist(complex_DiffExprProtein$medianLog2FC, breaks = 50)
  hist(complex_DiffExprProtein$global_medianLog2FC, breaks = 50)
dev.off()

#' ## Aggregate differential expression results to the complex level
complex_DiffExprComplex <- aggregateProteinTests(complex_DiffExprProtein)
saveRDS(complex_DiffExprComplex, "complex_DiffExprComplex.rda")

pdf("complex_DiffExprComplex_stats.pdf", width = 4, height=4)
  hist(complex_DiffExprComplex$pVal, breaks = 100)
  hist(complex_DiffExprComplex$global_pVal, breaks = 100)
  hist(complex_DiffExprComplex$medianLog2FC, breaks = 50)
  hist(complex_DiffExprComplex$global_medianLog2FC, breaks = 50)
dev.off()

#' # Make volcano plots
#' General volcano plots
plotVolcano(complex_DiffExprPep, PDF = T, pBHadj_cutoff = 0.05, name = "complex_DiffExprPep")
plotVolcano(complex_DiffExprProtein, PDF = T, pBHadj_cutoff = 0.05, name = "complex_DiffExprProtein")
plotVolcano(complex_DiffExprComplex, PDF = T, pBHadj_cutoff = 0.05, name = "complex_DiffExprComplex")

#' Get differential complexes
allComplexes <- unique(complex_DiffExprComplex$complex_id)
diffComplexes <- unique(complex_DiffExprComplex[pBHadj<0.05][abs(medianLog2FC)>1]$complex_id)
write.table(allComplexes, "allComplexes.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(diffComplexes, "diffComplexes.txt", col.names = F, row.names = F, quote = F, sep = "\t")

#' Volcanoplts highlighting different complexes
library(ggrepel)
highlight_complex <- diffComplexes[1]
plotVolcano(complex_DiffExprComplex, highlight=highlight_complex, PDF = T, pBHadj_cutoff = 0.05, name = paste0("complex_DiffExprComplex_",highlight_complex))

#' # Plot some differential protein features

pdf("diffComplexes.pdf", width=8, height=7)
for (id in diffComplexes) {
  plotFeatures(feature_table = scoredDataAll,
               traces = protein_traces_list,
               feature_id = id,
               design_matrix=design_matrix,
               calibration=calibration,
               peak_area = T,
               legend = F,
               onlyBest = F, 
               monomer_MW=T)
}
dev.off()

pdf("diffComplexes_entryName.pdf", width=8, height=7)
for (id in diffComplexes) {
  plotFeatures(feature_table = scoredDataAll,
               traces = protein_traces_list,
               feature_id = id,
               annotation_label = "Entry_name",
               design_matrix=design_matrix,
               calibration=calibration,
               peak_area = T,
               legend = F,
               onlyBest = F, 
               monomer_MW=T)
}
dev.off()


interestingComplexes <- c("COP9|proteasome|CCT|Septin|SMN|Prefoldin|snRNP|splice|Splice|DISC|CASP|EIF|RNA|Apop")
targets <- unique(scoredDataAll[grep(interestingComplexes,complex_name)]$complex_id)
pdf("generalComplexes_entryName.pdf", width=8, height=7)
for (id in targets) {
  plotFeatures(feature_table = scoredDataAll,
               traces = protein_traces_list,
               feature_id = id,
               annotation_label = "Entry_name",
               design_matrix=design_matrix,
               calibration=calibration,
               peak_area = T,
               legend = F,
               onlyBest = F, 
               monomer_MW=T)
}
dev.off()


