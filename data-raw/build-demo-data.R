library(data.table)
library(phyloseq)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript data-raw/build-demo-data.R /path/to/MicroFun_demonstration_data_filtered /path/to/MicroFun")
source_dir <- normalizePath(args[1], mustWork = TRUE)
package_dir <- normalizePath(args[2], mustWork = TRUE)
dir.create(file.path(package_dir, "data"), showWarnings = FALSE)
manifest <- read.csv(file.path(source_dir, "manifest.csv"), stringsAsFactors = FALSE)
for (i in seq_len(nrow(manifest))) {
  path <- file.path(source_dir, manifest$path[i])
  if (!file.exists(path) || file.info(path)$size != manifest$bytes[i] ||
      unname(tools::md5sum(path)) != manifest$md5[i]) stop("Source checksum mismatch: ", manifest$path[i])
}
e <- new.env()
load(file.path(source_dir,"phyloseq_asv.RData"),envir=e)
# fread's gzip support is available in the build environment.
read_input <- function(path) as.data.frame(fread(path,sep="\t",check.names=FALSE,showProgress=FALSE))
core <- as.list(e)
core$standard <- lapply(setNames(c("KO","EC","KEGG_2026-11-04","MetaCyc"),c("KO","EC","KEGG","MetaCyc")),function(n) read_input(file.path(source_dir,"standard",paste0(n,".tsv.gz"))))
core$selections <- list(
  KO=read.csv(file.path(source_dir,"MicroFun_selected_KOs.csv"),stringsAsFactors=FALSE)$KO,
  standard_KO=read.csv(file.path(source_dir,"standard/KO_figure_features.csv"),stringsAsFactors=FALSE)$feature,
  standard_KEGG=read.csv(file.path(source_dir,"standard/KEGG_figure_features.csv"),stringsAsFactors=FALSE)$feature)
# Genus selection from the source tutorial at the recorded commit.
core$selections$genera <- c(
  "Saccharopolyspora", "Collinsella", "Pseudoxanthomonas", "Streptomyces",
  "Shouchella", "Flavobacterium", "Pseudomonas", "Cloacibacterium",
  "Pedobacter", "Dyadobacter", "Mycovorax", "Devosia",
  "Schlegelella", "Cellulomonas", "Brachybacterium", "Thermobispora",
  "Geobacillus", "Brevibacillus", "Mammaliicoccus", "Weissella",
  "Bacteroides", "Faecalibacterium", "Tyzzerella", "Veillonella",
  "Thomasclavelia", "Streptococcus", "Lactobacillus"
)
core$KO_display_labels <- read.csv(file.path(source_dir,"KO_display_labels.csv"),stringsAsFactors=FALSE)
core$filter_summary <- read.csv(file.path(source_dir,"filter_summary.csv"),stringsAsFactors=FALSE)
core$removed_ASVs <- read.csv(file.path(source_dir,"removed_ASVs.csv"),stringsAsFactors=FALSE)
core$source_manifest <- read.csv(file.path(source_dir,"manifest.csv"),stringsAsFactors=FALSE)
core$provenance <- list(source_study_doi="10.1016/j.cell.2019.07.008",source_repository="https://github.com/bioscinema/MicroFun_website_tutorials",source_commit="998e0473ec492e1a5f155dd67a419354e12e8558",sample_ids="S001-S068",KEGG_snapshot="2026-11-04",packaging_date="2026-11-22",PICRUSt2_version=NA_character_,reference_database=NA_character_,data_license="Not established in the supplied source; confirm before redistribution")
microfun_demo <- core
save(microfun_demo,file=file.path(package_dir,"data/microfun_demo.rda"),compress="xz",compression_level=6)
rm(core,microfun_demo);gc()
for (kind in c("KO","EC","KEGG","MetaCyc")) {
 filename <- if(kind=="KEGG") "KEGG_2026-11-04" else kind
 name <- paste0("microfun_demo_",kind)
 obj <- read_input(file.path(source_dir,"stratified",paste0(filename,"_contrib.tsv.gz")))
 assign(name,obj)
 save(list=name,file=file.path(package_dir,"data",paste0(name,".rda")),compress="xz",compression_level=6)
 cat(name,nrow(obj),"rows;",file.info(file.path(package_dir,"data",paste0(name,".rda")))$size,"bytes\n");flush.console()
 rm(list=name);rm(obj);gc()
}
