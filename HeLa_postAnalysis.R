# global stats
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rds")
protein_list <- lapply(pepTracesList_filtered, function(t){unique(t$trace_annotation$protein_id)})
peptide_list <- lapply(pepTracesList_filtered, function(t){unique(t$trace_annotation$id)})

length(Reduce(intersect, protein_list))
length(Reduce(intersect, peptide_list))

length(unique(unlist(protein_list)))
length(unique(unlist(peptide_list)))

# AMF analysis
protTraces <- readRDS("protein_traces_list.rds")
diffAssemblyState <- readRDS("diffAssemblyState.rds")

fold_change_cutoff=0.3

withAssemblyChanges <- subset(diffAssemblyState, (abs(meanDiff) > fold_change_cutoff) & (betaQval < 0.05))
withAssemblyChanges

trace_ann_sub <- subset(protTraces[[1]]$trace_annotation[protein_id %in% withAssemblyChanges$protein_id])
trace_ann_sub

length(unique(trace_ann_sub$protein_id)) == length(unique(withAssemblyChanges$protein_id))

info_summary <- subset(trace_ann_sub, select=c("protein_id","Entry_name","Protein_names","Gene_names"))
info_summary <- merge(info_summary, withAssemblyChanges, by="protein_id")

info_summary <- info_summary[order(meanDiff)]
info_summary$meanDiff <- format(info_summary$meanDiff, scientific = FALSE, digits = 2) 
info_summary$meanAMF1 <- format(info_summary$meanAMF1, scientific = FALSE, digits = 2) 
info_summary$meanAMF2 <- format(info_summary$meanAMF2, scientific = FALSE, digits = 2) 
info_summary$betaPval <- format(info_summary$betaPval, scientific = TRUE, digits = 2) 
info_summary$betaPval_BHadj <- format(info_summary$betaPval_BHadj, scientific = TRUE, digits = 2) 
info_summary$betaQval <- format(info_summary$betaQval, scientific = TRUE, digits = 2) 
fwrite(info_summary, "diffAMF_info_summary.tsv", sep="\t", quote = F, col.names = T, row.names = F)


DAVID_more_assembled_proteins <- fread("DAVID_more_assembled_proteins.txt")
plotEnrichment(DAVID_more_assembled_proteins,name="DAVID_more_assembled_proteins")
plotEnrichmentPval(DAVID_more_assembled_proteins,name="DAVID_more_assembled_proteins_pval")

# Protein features
proteinFeaturesFiltered_DiffExprProtein <- readRDS("proteinFeaturesFiltered_DiffExprProtein.rda")
protein_annotation <- unique(do.call(rbind,(lapply(protTraces, function(t){
  subset(t$trace_annotation, select=c("protein_id","protein_mw"))
}))))
proteinFeaturesFiltered_DiffExprProtein_ann <- merge(proteinFeaturesFiltered_DiffExprProtein, protein_annotation, 
                                                     by.x="feature_id", by.y="protein_id", all.x=T, all.y=F)
calibration <- readRDS("calibration.rds")
proteinFeaturesFiltered_DiffExprProtein_ann[,monomer_fraction:=calibration$MWtoFraction(2*protein_mw)]
proteinFeaturesFiltered_DiffExprProtein_ann[, monomeric := ifelse(monomer_fraction>apex, FALSE, TRUE)]

monomers_changing <- unique(proteinFeaturesFiltered_DiffExprProtein_ann[monomeric==T][pBHadj<0.05][abs(medianLog2FC)>1]$feature_id)
monomers_all <- unique(proteinFeaturesFiltered_DiffExprProtein_ann[monomeric==T]$feature_id)

assembled_changing <- unique(proteinFeaturesFiltered_DiffExprProtein_ann[monomeric==F][pBHadj<0.05][abs(medianLog2FC)>1]$feature_id)
assembled_all <- unique(proteinFeaturesFiltered_DiffExprProtein_ann[monomeric==F]$feature_id)

monomers_changing <- monomers_changing[! monomers_changing %in% assembled_all]
monomers_constant <- monomers_all[! monomers_all %in% c(monomers_changing,assembled_all)]

assembled_changing <- assembled_changing[! assembled_changing %in% monomers_all]
assembled_constant <- assembled_all[! assembled_all %in% c(assembled_changing,monomers_all)]     

fraction_monomers_changing <- length(monomers_changing)/length(monomers_constant)
fraction_assembled_changing <- length(assembled_changing)/length(assembled_constant)

differentialMonomers <- proteinFeaturesFiltered_DiffExprProtein_ann[monomeric==T][pBHadj<0.05][abs(medianLog2FC)>1]$feature_id
write.table(differentialMonomers, "differentialMonomers.txt", sep="\t", quote=F, col.names = F, row.names = F)

##############
# compare to yansheng
hela_proteome <- fread("Hela_selected_genes_Proteomics.txt")
setnames(hela_proteome,"V1","feature_id")
hela_proteome <- melt(hela_proteome, id.vars = "feature_id", value.name = "log10", variable.name = "cell_line")
proteome_annotation_table <- data.table(
  cell_line = c("Prot_01_log10","Prot_02_log10","Prot_03_log10",
                "Prot_04_log10","Prot_05_log10","Prot_06_log10",
                "Prot_07_log10","Prot_08_log10","Prot_09_log10",
                "Prot_10_log10","Prot_11_log10","Prot_12_log10",
                "Prot_13_log10","Prot_14_log10"),
  subtype = c("Kyoto","CCL2","Kyoto",
              "Kyoto","S3","CCL2",
              "CCL2","Kyoto","Kyoto",
              "Kyoto","CCL2","CCL2",
              "CCL2","CCL2")
)
hela_proteome <- merge(hela_proteome, proteome_annotation_table, by="cell_line")
hela_proteome <- subset(hela_proteome, subtype %in% c("Kyoto","CCL2"))
hela_proteome[, sample := gsub("Prot_","sample_",cell_line)]
hela_proteome[, sample := gsub("_log10","",sample)]
hela_proteome[, cell_line := NULL]

hela_transcriptome <- fread("Hela_selected_genes_mRNA.txt")
setnames(hela_transcriptome,"V1","feature_id")
hela_transcriptome <- melt(hela_transcriptome, id.vars = "feature_id", value.name = "log10", variable.name = "cell_line")
transcriptome_annotation_table <- data.table(
  cell_line = c("mRNA_01_log10","mRNA_02_log10","mRNA_03_log10",
                "mRNA_04_log10","mRNA_05_log10","mRNA_06_log10",
                "mRNA_07_log10","mRNA_08_log10","mRNA_09_log10",
                "mRNA_10_log10","mRNA_11_log10","mRNA_12_log10",
                "mRNA_13_log10","mRNA_14_log10"),
  subtype = c("Kyoto","CCL2","Kyoto",
              "Kyoto","S3","CCL2",
              "CCL2","Kyoto","Kyoto",
              "Kyoto","CCL2","CCL2",
              "CCL2","CCL2")
)
hela_transcriptome <- merge(hela_transcriptome, transcriptome_annotation_table, by="cell_line")
hela_transcriptome <- subset(hela_transcriptome, subtype %in% c("Kyoto","CCL2"))
hela_transcriptome[, sample := gsub("mRNA_","sample_",cell_line)]
hela_transcriptome[, sample := gsub("_log10","",sample)]
hela_transcriptome[, cell_line := NULL]
setnames(hela_transcriptome, "log10", "mRNA_log10")

hela_cnv <- fread("Hela_selected_genes_CNV.txt")
setnames(hela_cnv,"V1","feature_id")
hela_cnv <- melt(hela_cnv, id.vars = "feature_id", value.name = "log2", variable.name = "cell_line")
cnv_annotation_table <- data.table(
  cell_line = c("CNV_01","CNV_02","CNV_03",
                "CNV_04","CNV_05","CNV_06",
                "CNV_07","CNV_08","CNV_09",
                "CNV_10","CNV_11","CNV_12",
                "CNV_13","CNV_14"),
  subtype = c("Kyoto","CCL2","Kyoto",
              "Kyoto","S3","CCL2",
              "CCL2","Kyoto","Kyoto",
              "Kyoto","CCL2","CCL2",
              "CCL2","CCL2")
)
hela_cnv <- merge(hela_cnv, cnv_annotation_table, by="cell_line")
hela_cnv <- subset(hela_cnv, subtype %in% c("Kyoto","CCL2"))
hela_cnv[, sample := gsub("CNV_","sample_",cell_line)]
hela_cnv[, cell_line := NULL]

hela_merged <- merge(hela_proteome,hela_cnv, by=c("feature_id","sample","subtype"))
hela_merged <- merge(hela_merged,hela_transcriptome, by=c("feature_id","sample","subtype"))

hela_merged[,corr_prot_cnv := cor.test( ~ log10 + log2,
                               method = "spearman",
                               continuity = FALSE,
                               conf.level = 0.95)$estimate, by=c("feature_id")]
hela_merged[,corr_prot_cnv_pval := cor.test( ~ log10 + log2,
                                        method = "spearman",
                                        continuity = FALSE,
                                        conf.level = 0.95)$p.value, by=c("feature_id")]

hela_merged[,corr_prot_rna := cor.test( ~ log10 + mRNA_log10,
                               method = "spearman",
                               continuity = FALSE,
                               conf.level = 0.95)$estimate, by=c("feature_id")]
hela_merged[,corr_cnv_rna := cor.test( ~ log2 + mRNA_log10,
                                        method = "spearman",
                                        continuity = FALSE,
                                        conf.level = 0.95)$estimate, by=c("feature_id")]

hela_merged_stat <- unique(subset(hela_merged, select=c("feature_id","corr_prot_cnv","corr_prot_rna","corr_cnv_rna")))

proteins_assebled <- unique(proteinFeaturesFiltered_DiffExprProtein_ann[monomeric==F]$feature_id)
proteins_monomeric <- unique(proteinFeaturesFiltered_DiffExprProtein_ann[monomeric==T]$feature_id)

only_assembled <- proteins_assebled[! proteins_assebled %in% proteins_monomeric]
mixed_assembled <- intersect(proteins_assebled,proteins_monomeric)
only_monomeric <- proteins_monomeric[! proteins_monomeric %in% proteins_assebled]
hela_merged_stat[,complex_in:=NA]
hela_merged_stat[,complex_in:=ifelse(feature_id %in% mixed_assembled, "mixed", complex_in)]
hela_merged_stat[,complex_in:=ifelse(feature_id %in% only_assembled, "only_assembled", complex_in)]
hela_merged_stat[,complex_in:=ifelse(feature_id %in% only_monomeric, "only_monomeric", complex_in)]

wilcox_stats <- data.table(
    wilcox_monoVsComplex_CnvProt = wilcox.test(hela_merged_stat[complex_in=="only_monomeric"]$corr_prot_cnv, 
                                               hela_merged_stat[complex_in=="only_assembled"]$corr_prot_cnv,
                                               paired = F, alternative = "two.sided")$p.value,
    wilcox_monoVsMixed_CnvProt = wilcox.test(hela_merged_stat[complex_in=="only_monomeric"]$corr_prot_cnv, 
                                             hela_merged_stat[complex_in=="mixed"]$corr_prot_cnv,
                                             paired = F, alternative = "two.sided")$p.value,
    wilcox_mixedVsComplex_CnvProt = wilcox.test(hela_merged_stat[complex_in=="mixed"]$corr_prot_cnv, 
                                                hela_merged_stat[complex_in=="only_assembled"]$corr_prot_cnv,
                                                paired = F, alternative = "two.sided")$p.value,
    wilcox_monoVsComplex_RnaProt = wilcox.test(hela_merged_stat[complex_in=="only_monomeric"]$corr_prot_rna, 
                                               hela_merged_stat[complex_in=="only_assembled"]$corr_prot_rna,
                                               paired = F, alternative = "two.sided")$p.value,
    wilcox_monoVsMixed_RnaProt = wilcox.test(hela_merged_stat[complex_in=="only_monomeric"]$corr_prot_rna, 
                                             hela_merged_stat[complex_in=="mixed"]$corr_prot_rna,
                                             paired = F, alternative = "two.sided")$p.value,
    wilcox_mixedVsComplex_RnaProt = wilcox.test(hela_merged_stat[complex_in=="mixed"]$corr_prot_rna, 
                                                hela_merged_stat[complex_in=="only_assembled"]$corr_prot_rna,
                                                paired = F, alternative = "two.sided")$p.value,
    wilcox_monoVsComplex_CnvRna = wilcox.test(hela_merged_stat[complex_in=="only_monomeric"]$corr_cnv_rna, 
                                              hela_merged_stat[complex_in=="only_assembled"]$corr_cnv_rna,
                                              paired = F, alternative = "two.sided")$p.value,
    wilcox_monoVsMixed_CnvRna = wilcox.test(hela_merged_stat[complex_in=="only_monomeric"]$corr_cnv_rna, 
                                            hela_merged_stat[complex_in=="mixed"]$corr_cnv_rna,
                                            paired = F, alternative = "two.sided")$p.value,
    wilcox_mixedVsComplex_CnvRna = wilcox.test(hela_merged_stat[complex_in=="mixed"]$corr_cnv_rna, 
                                               hela_merged_stat[complex_in=="only_assembled"]$corr_cnv_rna,
                                               paired = F, alternative = "two.sided")$p.value
)

fwrite(wilcox_stats,"wilcox_stats_spearmanCorr_complexIn.txt", sep="\t", quote=F, row.names = F)

hela_merged_stat.m <- melt(hela_merged_stat, id.vars = c("feature_id","complex_in"))
hela_merged_stat.m <- subset(hela_merged_stat.m, ! is.na(complex_in))

hela_merged_stat.m$complex_in <- factor(hela_merged_stat.m$complex_in, levels=c("only_monomeric","mixed","only_assembled"))
#library("RColorBrewer")
pdf("spearmanCorr_complexIn_stats.pdf", width=10, height = 4)
  ggplot(hela_merged_stat.m, aes(x=complex_in, y=value, fill=complex_in)) + 
    #geom_violin() +
    geom_boxplot() + # width=0.1
    facet_wrap(. ~ variable) +
    theme_classic() +
    scale_fill_manual(values = c("#F7F7F7","#EF8A62","#67A9CF")) +
    theme(axis.title.x=element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())
dev.off()

hela_merged_P12814 <- hela_merged[feature_id=="P12814"]
pdf("HelaCNVprotCorr_P12814.pdf", width = 4, height = 3)
ggplot(hela_merged_P12814, aes(x=log2, y=log10)) + 
  geom_point(aes(colour = subtype)) +
  theme_classic() +
  labs(x="CNV (log2)",y="proteomics (log10)") +
  geom_smooth(method='lm', se = F, color='black', size=0.5) +
  annotate("text", x = 0.2, y = 5.1, label = paste0("p-value = ", round(unique(hela_merged_P12814$corr_prot_cnv_pval), digits = 5)))
dev.off()

hela_merged_Q969U7 <- hela_merged[feature_id=="Q969U7"]
pdf("HelaCNVprotCorr_Q969U7.pdf", width = 4, height = 3)
ggplot(hela_merged_Q969U7, aes(x=log2, y=log10)) + 
  geom_point(aes(colour = subtype)) +
  theme_classic() +
  labs(x="CNV (log2)",y="proteomics (log10)") +
  geom_smooth(method='lm', se = F, color='black', size=0.5) +
  annotate("text", x = 0.4, y = 4.6, label = paste0("p-value = ", round(unique(hela_merged_Q969U7$corr_prot_cnv_pval), digits = 5)))
dev.off()


local_vs_global_test_protein <- readRDS("local_vs_global_test_protein.rds")
protFeatureStat_Q969U7 <- proteinFeaturesFiltered_DiffExprProtein_ann[feature_id=="Q969U7"]
rearrangeFeatureStat_Q969U7 <- local_vs_global_test_protein[feature_id=="Q969U7"]
proteinStats_Q969U7 <- data.table(
  apex=c(protFeatureStat_Q969U7$apex),
  medianLog2FC=c(protFeatureStat_Q969U7$medianLog2FC),
  pBHadj=c(protFeatureStat_Q969U7$pBHadj),
  global_medianLog2FC_imp=c(protFeatureStat_Q969U7$global_medianLog2FC_imp),
  global_pBHadj=c(protFeatureStat_Q969U7$global_pBHadj),
  localVsGlobal_medianDiff=rearrangeFeatureStat_Q969U7$medianDiff,
  localVsGlobal_pBHadj=rearrangeFeatureStat_Q969U7$pBHadj
)

proteinStats_Q969U7$apex <- as.character(proteinStats_Q969U7$apex)

pdf("proteinStats_Q969U7.pdf", height = 3, width = 3)
ggplot(proteinStats_Q969U7, aes(x=apex, y=medianLog2FC, fill=-log10(pBHadj))) + 
  geom_col() +
  geom_text(aes(label=round(medianLog2FC, digits = 2)), vjust=-.5) + 
  theme_classic()+
  ylim(-3,3) +
  scale_fill_gradient(low = "white", high="blue", limits=c(0,12)) +
  theme(legend.position="bottom")

ggplot(proteinStats_Q969U7, aes(x="global",y=global_medianLog2FC_imp, fill=-log10(global_pBHadj))) + 
  geom_col() +
  geom_text(aes(label=round(global_medianLog2FC_imp, digits = 2)), vjust=-.5) + 
  theme_classic()+
  ylim(-3,3) +
  scale_fill_gradient(low = "white", high="blue", limits=c(0,12)) +
  theme(legend.position="bottom")

ggplot(proteinStats_Q969U7, aes(x=apex,y=localVsGlobal_medianDiff, fill=-log10(localVsGlobal_pBHadj))) + 
  geom_col() +
  geom_text(aes(label=round(localVsGlobal_medianDiff, digits = 2)), vjust=-.5) + 
  theme_classic() +
  ylim(-0.5,0.5) +
  scale_fill_gradient(low = "white", high="blue", limits=c(0,12)) +
  theme(legend.position="bottom")
dev.off()


# diffProtein feature enrichment
david_diffProteins <- fread("david_diffProteins.txt")
plotEnrichment(david_diffProteins,name="david_diffProteins")
plotEnrichmentPval(david_diffProteins,name="david_diffProteinss_pval")

# more beautiful summary plots
proteinFeaturesFiltered <- readRDS("proteinFeaturesFiltered.rda")
proteinFeaturesFiltered <- proteinFeaturesFiltered[grep("DECOY", protein_id, invert = T)]
proteinFeaturesFiltered[, n_subfeatures := .N, by = c("protein_id")]

proteinFeaturesFiltered_stat <- unique(subset(proteinFeaturesFiltered, select=c("protein_id","n_subfeatures")))

pdf("proteinFeaturesFiltered_nSubfeatures.pdf", width=3, height=4)
q <- ggplot(proteinFeaturesFiltered_stat,aes(x=n_subfeatures)) +
  stat_bin(binwidth=1) +
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5) +
  labs(x="N co-elution features",y="N proteins") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(q)
dev.off()

proteinFeatures <- readRDS("proteinFeatures.rda")
featuresScored <- calculateCoelutionScore(proteinFeatures)
pdf("proteinFeatures_niceStats.pdf", height=3,width=3)
  qvalueFeaturesScored <- calculateQvalue(featuresScored, name="qvalueStats", PDF=F)
  qvalueFeaturesScoredStats <- qvaluePositivesPlot(qvalueFeaturesScored, name="qvaluePositivesPlot",PDF=F)
dev.off()


# more beautiful summary plots on complex level
scoredDataAll <- readRDS("scoredDataAll.rds")
scoredDataAll <- scoredDataAll[grep("DECOY", complex_id, invert = T)]
scoredDataAll[, n_subfeatures := .N, by = c("complex_id")]

scoredDataAll_stat <- unique(subset(scoredDataAll, select=c("complex_id","n_subfeatures")))
scoredDataAll_stat[, n_subfeatures := ifelse(n_subfeatures>=7, 7, n_subfeatures)]

pdf("complexFeatures_nSubfeatures.pdf", width=3, height=4)
q <- ggplot(scoredDataAll_stat,aes(x=n_subfeatures)) +
  stat_bin(binwidth=1) +
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5) +
  labs(x="N co-elution features",y="N proteins") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_discrete(limits=seq(1,7,1), 
                   labels=c("1","2","3","4","5","6",expression(phantom(x) >= "7")))
print(q)
dev.off()

# Complex detection stats
complexFeaturesBest <- readRDS("complexFeaturesBest.rda")
complexFeaturesBestScored <- calculateCoelutionScore(complexFeaturesBest)
pdf("complexFeaturesBest_niceStats.pdf", height=3,width=3)
complexFeaturesBestqvalueFeaturesScored <- calculateQvalue(complexFeaturesBestScored, name="qvalueStats", PDF=F)
complexFeaturesBestqvalueFeaturesScoredStats <- qvaluePositivesPlot(complexFeaturesBestqvalueFeaturesScored, name="qvaluePositivesPlot",PDF=F)
dev.off()

# Table of differential complexes
complex_DiffExprComplex <- readRDS("complex_DiffExprComplex.rda")

complex_DiffExprComplex_sig <- complex_DiffExprComplex[pBHadj<0.05][abs(medianLog2FC)>1]
complex_DiffExprComplex_sig <- subset(complex_DiffExprComplex_sig, select=c("complex_id","apex","Nproteins","medianLog2FC","pVal","pBHadj"))

corumAnn <- readRDS("corumHypotheses.rds")
corumAnn <- unique(subset(corumAnn, select = c("complex_id","complex_name")))
complex_DiffExprComplex_sig <- merge(complex_DiffExprComplex_sig,corumAnn,by="complex_id", all.x=T, all.y=F)

complex_DiffExprComplex_sig <- complex_DiffExprComplex_sig[order(medianLog2FC)]
complex_DiffExprComplex_sig$medianLog2FC <- format(complex_DiffExprComplex_sig$medianLog2FC, scientific = FALSE, digits = 2) 
complex_DiffExprComplex_sig$pVal <- format(complex_DiffExprComplex_sig$pVal, scientific = TRUE, digits = 2) 
complex_DiffExprComplex_sig$pBHadj <- format(complex_DiffExprComplex_sig$pBHadj, scientific = TRUE, digits = 2) 

fwrite(complex_DiffExprComplex_sig, "complex_DiffExprComplex_sig_info_summary.tsv", sep="\t", quote = F, col.names = T, row.names = F)

library(ggrepel)
plotVolcano(complex_DiffExprComplex, highlight=complex_DiffExprComplex_sig$complex_id, PDF = T, pBHadj_cutoff = 0.05, name = paste0("complex_DiffExprComplex_allSignificant"))
plotVolcano(complex_DiffExprComplex, highlight="1345", PDF = T, pBHadj_cutoff = 0.05, name = paste0("complex_DiffExprComplex_Septin"))


# get diff proteins

proteinFeaturesFiltered_DiffExprProtein 
local_vs_global_test_protein <- readRDS("local_vs_global_test_protein.rds")

n_proteins <- length(unique(proteinFeaturesFiltered_DiffExprProtein$feature_id))
n_protein_features <- nrow(proteinFeaturesFiltered_DiffExprProtein)
n_differential_proteins <- length(unique(proteinFeaturesFiltered_DiffExprProtein[pBHadj<0.05][abs(medianLog2FC)>1]$feature_id))
n_differential_protein_features <- nrow(proteinFeaturesFiltered_DiffExprProtein[pBHadj<0.05][abs(medianLog2FC)>1])
n_proteins_rearranged <- length(unique(local_vs_global_test_protein[pBHadj<0.05][abs(medianDiff)>0.3]$feature_id))
n_protein_features_rearranged <- nrow(local_vs_global_test_protein[pBHadj<0.05][abs(medianDiff)>0.3])

proteinFeatureDiffStats <- data.table(
  level=c("proteins","proteins","proteins","features","features","features"),
  quant=c("n","n differential","n rearranged","n","n differential","n rearranged"),
  count=c(n_proteins, n_differential_proteins, n_proteins_rearranged,
          n_protein_features, n_differential_protein_features, n_protein_features_rearranged)
)

proteinFeatureDiffStats$level <- factor(proteinFeatureDiffStats$level, levels = c("proteins","features"))

pdf("proteinFeatureDiffStats.pdf", width=7, height=3)
ggplot(proteinFeatureDiffStats, aes(x=quant, y=count, fill=level, group=level)) + 
  geom_col() + 
  facet_wrap(. ~ level, scales = "free") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  geom_text(aes(label=count), vjust=-.1) + 
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))
dev.off()


#####


actin <- corumComplexHypotheses[complex_id=="5177"]$protein_id
actinTraces <- subset(protTraces, actin)
plot(actinTraces)
