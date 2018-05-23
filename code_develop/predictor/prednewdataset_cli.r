
library(R.utils)

#Rscript prednewdataset_cli.r --pretty_home . --species hsa --exp_file /fsclinic/common/analisis/validation_data/data/BRCA_genes_vals_BN.txt --design_file /fsclinic/common/analisis/validation_data/data/BRCA_Normal-Basal_ED.txt --design_type categorical --cond1 Tumor --cond2 Normal --model_file /home/aamadoz/001_pathipredRNAseq1015/example_results/model.RData --output_folder /home/aamadoz/001_pathipredRNAseq1015/example_results_prediction/

#### INPUT DATA

args <- commandArgs(trailingOnly = F, asValues = T,excludeEnvVars = F)
pretty_home <- paste0(dirname(normalizePath(args[["file"]])),"/")

species <- args[['species']]

exp_file <- args[['exp_file']]

design_file <- args[['design_file']]
design_type <- args[['design_type']]
cond1 <- args[['cond1']]
cond2 <- args[['cond2']]

model_file <- args[['model_file']]

decompose <- args[['decompose']]

output_folder <- args[['output_folder']]

# parse

if(is.null(decompose)) decompose <- F
if(design_type=="categorical") {
  if(is.null(cond1) | is.null(cond2)) {
    cond1 <- "Tumor"
    cond2 <- "Normal"
  }
}

cat("\nWELCOME TO PRETTYWAYS (perhaps CYGNUS)!!!\n")
cat("\nDetected params: \n")
cat("\tpretty_home: ",pretty_home,"\n")
cat("\tspecies: ",species,"\n")
cat("\texp_file: ",exp_file,"\n")
cat("\tdesign_file: ",design_file,"\n")
cat("\tdesign_type: ",design_type,"\n")
cat("\tcond1: ",cond1,"\n")
cat("\tcond2: ",cond2,"\n")
cat("\tmodel_file: ",model_file,"\n")
cat("\tdecompose: ",decompose,"\n")
cat("\toutput_folder: ",output_folder,"\n")
cat("\n\n")

# rm (list = ls ())
# pretty_home <- "/home/aamadoz/001_pathipredRNAseq1015/prettyways/"
# species <- "hsa"
# exp_file <- "/fsclinic/common/analisis/validation_data/data/BRCA_genes_vals_BN.txt"
# design_file <- "/fsclinic/common/analisis/validation_data/data/BRCA_Normal-Basal_ED.txt"
# cond1 <- "Tumor"
# cond2 <- "Normal"
# design_type <- "categorical"
# model_file <- "/home/aamadoz/001_pathipredRNAseq1015/example_results/model.RData"
# decompose <- F
# output_folder <- "/home/aamadoz/001_pathipredRNAseq1015/example_results/"

#### PREPARE DATA

source(paste0(pretty_home,"/prettyways.R"))
source(paste0(pretty_home,"/predict.r"))
source(paste0(pretty_home,"/stats.R"))
source(paste0(pretty_home,"/functions.r"))
source(paste0(pretty_home,"/web.r"))

dir.create(output_folder)
status <- function(value){
  write(value,file=paste0(output_folder,"/status.txt"))
}
status("0")

options(warn=-1)
load(paste0(pretty_home,"/files/pathigraphs_functions_AMEND.RData"))
fpathigraphs <- fpgs[c("hsa04014","hsa04015","hsa04010","hsa04012","hsa04310","hsa04350","hsa04390","hsa04370","hsa04630","hsa04068","hsa04071","hsa04024","hsa04151","hsa04150","hsa04110","hsa04210","hsa04115","hsa04620")]
pathigraph.genes <- all.needed.genes(fpathigraphs)
options(warn=0)
status("4")

#### LOAD DATA

cat("Loading data...\n")

exp <- read.table(exp_file,header=T,sep="\t",stringsAsFactors=F,row.names = 1)
des <- read.table(design_file,header=F,stringsAsFactors=F)
colnames(des) <- c("sample","group")
rownames(des) <- des$sample

sel_samples <- intersect(colnames(exp),rownames(des))
exp <- exp[,sel_samples]
des <- des[sel_samples,]

#### PREPROCESS DATA

cat("Scaling expression values...\n")

exp <- normalize.data(exp,by.quantiles = F,by.gene = F,percentil = F)
exp <- add.missing.genes(exp,genes=pathigraph.genes)

status("20")

#### RUN

# pretty results
cat("Propagating signaling...\n")
results <- prettyways(exp, fpathigraphs,verbose=F)

status("50")

#### ANALYSIS

if(decompose==T){
  path.vals <- results$all$path.vals
} else {
  path.vals <- results$all$effector.path.vals
}
n <- ncol(path.vals)

## predictor
cat("Predicting...\n")
load(model_file)
predres <- predict.newdataset(model[["model"]], t(path.vals), des, design_type) 

cat("Calculating prediction model statistics...\n")
pred_stats <- get.stats(predres[["svm_pred"]], model, design_type)

status("95")

#### SAVE RESULTS

cat("Creating report...\n")

save.prednewdataset.res(predres[["pred_res"]], pred_stats,results,path.vals,fpathigraphs,output_folder)

################## new report

results <- init.result("prettypred")

results <- add.input.param(results,"species",species)

results <- add.section(results,"Input params")
results <- add.download(results,"Expression file",basename(exp_file))
results <- add.download(results,"Design file",basename(design_file))
results <- add.param(results,"Design type",design_type)
results <- add.param(results,"Sample size",n)

results <- add.section(results,"Path values")
results <- add.download(results,"Path values","path_vals.txt")

results <- add.section(results,"Prediction")
results <- add.download(results,"Prediction results","prediction_results.txt")
results <- add.download(results,"Statistics","prediction_stats.txt")

write(render.xlm(results),file=paste0(output_folder,"/report.xml"))

unlink(paste0(output_folder,"/sifs4CellMaps"),recursive = T)

cat("[Finished]\n")

status("100")

