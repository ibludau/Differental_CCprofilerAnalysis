#' # Read Spectronaut output
data <-  fread("20181028_021159_SEC_HeLa_CCL2_Kyoto_peptides_swnormalized.tsv")
# change column names according to CCprofiler conventions
setnames(data,c("PG.ProteinAccessions","EG.ModifiedSequence","R.FileName","FG.Quantity"),c("protein_id","peptide_id","filename","intensity"))
# remove non-proteotypic peptides
data <- data[grep(";",data$protein_id, invert = TRUE)]
# reduce data to necessary columns
data_sub <- subset(data,select=c("protein_id","peptide_id","filename","intensity"))
# set negative values to zero > this is coming from data normalization > George
data_sub[, intensity := ifelse(intensity<0, 0, intensity)]

#' # Load annotation table
fraction_annotation <- fread("ccl2kyoto_sec_mw_spectronaut.csv")
setnames(fraction_annotation,c("run_id","sec_id"),c("filename","fraction_number"))
fraction_annotation[, sample := paste(c(condition_id, replicate_id), collapse = "_"), by=c("filename","fraction_number")]
fraction_annotation_sub <- subset(fraction_annotation, select=c("filename","fraction_number","sample"))

#' # Molecular weight calibration
molecularWeightCalibration <- unique(subset(fraction_annotation, select=c("sec_mw","fraction_number")))
setnames(molecularWeightCalibration,c("sec_mw","fraction_number"),c("std_weights_kDa","std_elu_fractions"))
molecularWeightCalibration[, std_weights_kDa := std_weights_kDa/1000]
calibration = calibrateMW(molecularWeightCalibration,
                          PDF=T,
                          plot=TRUE)
saveRDS(calibration,"calibration.rds")

#' # Create design matrix
design_matrix <- unique(subset(fraction_annotation, select=c("condition_id","replicate_id","sample")))
setnames(design_matrix, c("condition_id","replicate_id","sample"), c("Condition","Replicate","Sample_name"))
saveRDS(design_matrix,"design_matrix.rds")

#' # Import traces list
samples <- unique(fraction_annotation$sample)
# Import data as traces object for each sample
traces_list <- lapply(samples,function(x){
  message(x)
  ann <- subset(fraction_annotation_sub, sample==x)
  ann <- subset(ann, filename %in% data_sub$filename)
  ann <- subset(ann, select = c("filename","fraction_number"))
  setkey(ann,filename)
  data_sub_select <- subset(data_sub,filename %in% ann$filename)
  setkey(data_sub_select,filename)
  traces <- importPCPdata(input_data=data_sub_select,fraction_annotation=ann)
  return(traces)
})
names(traces_list) = samples
class(traces_list) <- "tracesList"
# remove input data to clean storage
rm(data,data_sub)
gc()

#' # Annotate traces with information from uniprot
pepTraces_raw <- annotateTraces(traces=traces_list,
                            trace_annotation=exampleTraceAnnotation,
                            traces_id_column = "protein_id",
                            trace_annotation_id_column = "Entry",
                            trace_annotation_mass_column = "Mass",
                            uniprot_mass_format = TRUE,
                            replace_whitespace = TRUE)

#' # Annotate traces with molecular weight calibration
pepTraces_raw <- annotateMolecularWeight(pepTraces_raw,
                                     calibration)

summary(pepTraces_raw)
saveRDS(pepTraces_raw, "pepTraces_raw.rds")
