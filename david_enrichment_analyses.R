#' # Do enrichment analysis in DAVID with all 3 GO direct terms

plotEnrichment <- function(david_enrichment, enrichment_cutoff = 1.5, Bonferroni_cutoff = 0.05, name="david_enrichment", PDF=T) {
  names(david_enrichment) <- gsub(" ", "_", names(david_enrichment))
  david_enrichment_sig <- subset(david_enrichment, Fold_Enrichment >= enrichment_cutoff & Bonferroni <= Bonferroni_cutoff)
  david_enrichment_sig[, neglog10Bonferroni := -log10(Bonferroni)]
  cat <- david_enrichment_sig$Term
  david_enrichment_sig$Term <- factor(david_enrichment_sig$Term, levels=cat)
  david_enrichment_sig[,name:=gsub("^GO:.*~","",Term)]
  david_enrichment_sig <- david_enrichment_sig[order(Fold_Enrichment)]
  david_enrichment_sig$name <- factor(david_enrichment_sig$name, levels = david_enrichment_sig$name)
  if(nrow(david_enrichment_sig) > 10) {
    david_enrichment_sig <- david_enrichment_sig[1:10]
  }
  
  q <- ggplot(data=david_enrichment_sig,aes(x=name,y=Fold_Enrichment, fill=neglog10Bonferroni)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.text.y  = element_text(angle = 30, hjust = 1)) +
    labs(fill='-log10Bonferroni',x='Category',y='Fold Enrichment') +
    coord_flip() +
    geom_text(aes(label=round(Fold_Enrichment,digits=1)), hjust=1.2,color="white") 
  if (PDF){
    pdf(paste0(name,"_enrichment_cutoff_",enrichment_cutoff,"_Bonferroni_cutoff_",Bonferroni_cutoff,".pdf"),height=4,width=6)
  }
  plot(q)
  if (PDF){
    dev.off()
  }
}

plotEnrichmentPval <- function(david_enrichment, enrichment_cutoff = 1.5, Bonferroni_cutoff = 0.05, name="david_enrichmentPval", PDF=T) {
  names(david_enrichment) <- gsub(" ", "_", names(david_enrichment))
  david_enrichment_sig <- subset(david_enrichment, Fold_Enrichment >= enrichment_cutoff & Bonferroni <= Bonferroni_cutoff)
  david_enrichment_sig[, neglog10Bonferroni := -log10(Bonferroni)]
  cat <- david_enrichment_sig$Term
  david_enrichment_sig$Term <- factor(david_enrichment_sig$Term, levels=cat)
  david_enrichment_sig[,name:=gsub("^GO:.*~","",Term)]
  david_enrichment_sig <- david_enrichment_sig[order(-Bonferroni)]
  david_enrichment_sig$name <- factor(david_enrichment_sig$name, levels = david_enrichment_sig$name)
  if(nrow(david_enrichment_sig) > 10) {
    david_enrichment_sig <- david_enrichment_sig[1:10]
  }
  
  q <- ggplot(data=david_enrichment_sig,aes(x=name,y=neglog10Bonferroni)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.text.y  = element_text(angle = 30, hjust = 1)) +
    labs(x='Category',y='-log10Bonferroni') +
    coord_flip() +
    geom_text(aes(label=round(neglog10Bonferroni,digits=0)), hjust=1.2,color="white") +
    geom_hline(yintercept=-log10(Bonferroni_cutoff), linetype="solid", size=0.5, colour="red")
  if (PDF){
    pdf(paste0(name,"_enrichment_cutoff_",enrichment_cutoff,"_Bonferroni_cutoff_",Bonferroni_cutoff,".pdf"),height=4,width=4)
  }
  plot(q)
  if (PDF){
    dev.off()
  }
}


david_diffProteinFeatures <- fread("david_diffProteinFeatures.txt")
david_diffProteinFeatures_global <- fread("david_diffProteinFeatures_global.txt")
david_allProteinFeatures_vs_human <- fread("david_allProteinFeatures_vs_human.txt")

plotEnrichment(david_diffProteinFeatures,name="david_diffProteinFeatures")
plotEnrichment(david_diffProteinFeatures_global,name="david_diffProteinFeatures_global")
plotEnrichment(david_allProteinFeatures_vs_human,name="david_allProteinFeatures_vs_human")


plotEnrichmentPval(david_diffProteinFeatures,name="davidPval_diffProteinFeatures")
plotEnrichmentPval(david_diffProteinFeatures_global,name="davidPval_diffProteinFeatures_global")
plotEnrichmentPval(david_allProteinFeatures_vs_human,name="davidPval_allProteinFeatures_vs_human")



