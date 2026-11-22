library(MicroFun)
reader <- getFromNamespace(".microfun_read_contributions", "MicroFun")
expect_error <- function(expr, text) {
  err <- tryCatch({force(expr);NULL}, error=identity)
  stopifnot(inherits(err,"error"),grepl(text,conditionMessage(err),fixed=TRUE))
}
expect_error(reader(),"exactly one")
expect_error(reader("missing.tsv",data.frame()),"exactly one")
expect_error(reader(stratified_data=matrix(1)),"data frame")
expect_error(reader(c("a","b")),"existing contribution file")
x <- data.table::data.table(sample="S001",taxon="ASV1",`function`="K00001",taxon_function_abun=2.5)
before <- data.table::copy(x)
y <- reader(stratified_data=x)
y[, taxon_function_abun := 9]
stopifnot(identical(x,before))
path <- tempfile(fileext=".tsv")
data.table::fwrite(x,path,sep="\t")
stopifnot(identical(reader(path),reader(stratified_data=as.data.frame(x))))
unlink(path)
# Run the same small differential analysis through both input routes.
samples <- paste0("S",seq_len(8))
asvs <- c("a","b")
meta <- data.frame(group=rep(c("Control","Treatment"),each=4),row.names=samples)
tax <- matrix(c("GenusA","GenusB"),ncol=1,dimnames=list(asvs,"Genus"))
ps <- phyloseq::phyloseq(
  phyloseq::otu_table(matrix(1L,2,8,dimnames=list(asvs,samples)),taxa_are_rows=TRUE),
  phyloseq::tax_table(tax),phyloseq::sample_data(meta))
contrib <- expand.grid(sample=samples,taxon=asvs,stringsAsFactors=FALSE)
contrib[["function"]] <- "K00001"
contrib$taxon_function_abun <- c(10,11,12,13,80,90,100,110,40,41,39,42,40,42,41,39)
path <- tempfile(fileext=".tsv")
data.table::fwrite(contrib,path,sep="\t")
args <- list(physeq=ps,group="group",direction="function",pathway="KO",
             reference="Control",daa_method="Maaslin2",threshold=1,verbose=FALSE)
a <- do.call(pathway_sdaa,c(list(stratified_path=path),args))
b <- do.call(pathway_sdaa,c(list(stratified_data=contrib),args))
stopifnot(isTRUE(all.equal(a,b)),length(a)>0L)
unlink(path)
