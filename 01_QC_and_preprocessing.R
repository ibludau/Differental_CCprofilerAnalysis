#' # Load data
pepTraces_raw <- readRDS("pepTraces_raw.rds")
design_matrix <- readRDS("design_matrix.rds")

#' # QC
#' # Assess how well SEC traces align across all samples
alignTraces(pepTraces_raw, plot = T, PDF=T)
#' # Assess total intensity as a proxi for extraction efficiency across samples
plotGlobalIntensities(pepTraces_raw, plot = T, PDF=T)

#' Find missing values 
#' (defined as having identifications in left and right neigbouring fractions):
pepTracesMV <- findMissingValues(pepTraces_raw,
                                 bound_left = 1,
                                 bound_right = 1,
                                 consider_borders = TRUE)

#' Impute NA values by fitting a spline:
pepTracesImp <- imputeMissingVals(pepTracesMV, method = "spline")

#' Plot imputation summary:
plotImputationSummary(pepTracesMV, pepTracesImp, PDF = T, 
                      plot_traces = T, max_n_traces = 2)

#' ## Filter by consecutive IDs
pepTracesConsIds <- filterConsecutiveIdStretches(pepTracesImp,
                                                 min_stretch_length = 2,
                                                 remove_empty = T)
saveRDS(pepTracesConsIds,"pepTracesConsIds.rds")

pepTracesConsIds_multiPep <- filterSinglePeptideHits(pepTracesConsIds)
saveRDS(pepTracesConsIds_multiPep,"pepTracesConsIds_multiPep.rds")

pepTraces_maxCorr <- filterByMaxCorr(pepTracesConsIds_multiPep,
                                     cutoff = 0.5, plot = T, PDF=T, name="maxCorrHist")

#' ## Calculate Sibling Peptide Correlation
pepTraces_maxCorr_SPC <- calculateSibPepCorr(pepTraces_maxCorr,
                                             plot = T, PDF = T,
                                             name = "SibPepCorr_densityplot_afterMaxCorrFilter")

#' Update traces with additional metrics for each fraction:
pepTracesList_filtered <- updateTraces(pepTraces_maxCorr_SPC)

#' Inspect traces list:
summary(pepTracesList_filtered)

saveRDS(pepTracesList_filtered, "pepTracesList_filtered.rds")

#' # Plot some examplary peptide traces
examples <- unique(pepTracesList_filtered[[1]]$trace_annotation$protein_id)[1:3]
for(test_proteins in examples){
  pepTest <- subset(pepTracesList_filtered, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
  plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF=T, name = paste0("PeptideTraces_",test_proteins))
}

#' # Global stats
traces_names <- c(
  "pepTraces_raw",
  "pepTracesConsIds",
  "pepTracesConsIds_multiPep",
  "pepTracesList_filtered"
)

countDT <- data.table(evidence=character(), sample=character(), count=integer(), traces=character())

for (traces in traces_names) {
  #traces <- readRDS(paste0(traces_name,".rda"))
  traces_name <- traces
  traces <- eval(as.name(traces))
  
  getIDs <- lapply(traces, function(x){
    peptides <- unique(x$trace_annotation$id)
    proteins <- unique(x$trace_annotation$protein_id)
    list(peptides=length(peptides),proteins=length(proteins))
  })
  
  count_dt <- as.data.table(getIDs)
  count_dt[, evidence := c("peptides", "proteins")]
  count_dt <- melt(count_dt, id.vars = "evidence", variable.name = "sample", value.name = "count")
  count_dt[,traces:=traces_name]
  count_dt[, count := unlist(count)]
  countDT <- rbind(countDT,count_dt)
}

traces_names_dt <- data.table(
  traces=c("pepTraces_raw",
           "pepTracesConsIds",
           "pepTracesConsIds_multiPep",
           "pepTracesList_filtered"), 
  traces_name=c("raw",
                "consecutive ids",
                "no single peptides",
                "max. correlation")
)

countDT <- merge(countDT, traces_names_dt, by="traces")

countDT$traces_name <- factor(countDT$traces_name, levels=traces_names_dt$traces_name)

write.table(countDT, "traces_stat_counts.txt", sep="\t",quote = F, col.names = T, row.names = F)

pdf("traces_stat_counts.pdf",height=6,width=10)
g <- ggplot(countDT,aes(x=traces_name,y=count,fill=evidence, group=sample)) +
  geom_bar(stat="identity") +
  facet_wrap(sample ~ evidence, scales="free",nrow=2) +
  theme_classic() +
  theme(legend.position="bottom") + theme(legend.title = element_blank()) +
  labs(y = "count") +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(g)
dev.off()

#' ## Quantify on protein level
protein_traces_list <- proteinQuantification(pepTracesList_filtered ,
                                             quantLevel="protein_id",
                                             topN = 2,
                                             keep_less = TRUE,
                                             rm_decoys = TRUE,
                                             use_sibPepCorr = FALSE,
                                             use_repPepCorr = FALSE,
                                             full_intersect_only = FALSE,
                                             verbose = FALSE)

#' ## Update fraction annotation for protein traces
protein_traces_list <- updateTraces(protein_traces_list)
saveRDS(protein_traces_list,"protein_traces_list.rds")
