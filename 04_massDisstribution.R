#' Load data
design_matrix <- readRDS("design_matrix.rds") 
protTraces <- readRDS("protein_traces_list.rds")

# evaluate protein mass distribution
protTraces_assembly <- annotateMassDistribution(protTraces)
saveRDS(protTraces_assembly,"protTraces_assembly.rds")

library(betareg)
library(lmtest)

#protTraces_assembly <- subset(protTraces_assembly, unique(protTraces_assembly[[1]]$trace_annotation$protein_id)[1:500], trace_subset_type = "protein_id")

diffAssemblyState <- getMassAssemblyChange(protTraces_assembly, design_matrix, plot=T, PDF=T)
saveRDS(diffAssemblyState,"diffAssemblyState.rds")

fold_change_cutoff=0.3

write.table(diffAssemblyState,paste0("diffAssemblyState.txt"),sep="\t",quote=F,row.names = F,col.names=T)

if (length(unique(design_matrix$Replicate)) > 1) {
  all_proteins_independentOfAssemblyState <- unique(diffAssemblyState$protein_id)
  proteins_noAssemblyChanges <- subset(diffAssemblyState, (abs(meanDiff) <= fold_change_cutoff) & (betaPval_BHadj < 0.05))$protein_id
  more_assembled_proteins <- subset(diffAssemblyState, (meanDiff < -fold_change_cutoff) & (betaPval_BHadj < 0.05))$protein_id
  less_assembled_proteins <- subset(diffAssemblyState, (meanDiff > fold_change_cutoff) & (betaPval_BHadj < 0.05))$protein_id
  proteins_withAssemblyChanges <- subset(diffAssemblyState, (abs(meanDiff) > fold_change_cutoff) & (betaPval_BHadj < 0.05))$protein_id
} else {
  all_proteins_independentOfAssemblyState <- unique(diffAssemblyState$protein_id)
  proteins_noAssemblyChanges <- subset(diffAssemblyState, (abs(meanDiff) <= fold_change_cutoff))$protein_id
  more_assembled_proteins <- subset(diffAssemblyState, (meanDiff < -fold_change_cutoff))$protein_id
  less_assembled_proteins <- subset(diffAssemblyState, (meanDiff > fold_change_cutoff))$protein_id
  proteins_withAssemblyChanges <- subset(diffAssemblyState, (abs(meanDiff) > fold_change_cutoff))$protein_id
  
}

write.table(all_proteins_independentOfAssemblyState,paste0("all_proteins_independentOfAssemblyState.txt"),sep="\t",quote=F,row.names = F,col.names=F)
write.table(proteins_noAssemblyChanges,paste0("proteins_noAssemblyChanges.txt"),sep="\t",quote=F,row.names = F,col.names=F)
write.table(more_assembled_proteins,paste0("more_assembled_proteins.txt"),sep="\t",quote=F,row.names = F,col.names=F)
write.table(less_assembled_proteins,paste0("less_assembled_proteins.txt"),sep="\t",quote=F,row.names = F,col.names=F)
write.table(proteins_withAssemblyChanges,paste0("proteins_withAssemblyChanges.txt"),sep="\t",quote=F,row.names = F,col.names=F)

#' # Plot proteins_withAssemblyChanges
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rds")
examples <- proteins_withAssemblyChanges
pdf("proteins_withAssemblyChanges_pepTraces.pdf", width=10, height=6)
for(test_proteins in examples){
  pepTest <- subset(pepTracesList_filtered, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
  plot(pepTest, design_matrix = design_matrix, log = FALSE, monomer_MW=T, legend = F, PDF=F, name = paste0("PeptideTraces_",test_proteins))
}
dev.off()
pdf("proteins_withAssemblyChanges_protTraces.pdf", width=10, height=6)
for(test_proteins in examples){
  protTest <- subset(protTraces, trace_subset_ids = test_proteins)
  plot(protTest, design_matrix = design_matrix, log = FALSE, monomer_MW=T, legend = T, PDF=F, name = paste0("PeptideTraces_",test_proteins))
}
dev.off()


