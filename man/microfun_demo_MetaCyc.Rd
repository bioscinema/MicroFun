\name{microfun_demo_MetaCyc}
\alias{microfun_demo_MetaCyc}
\docType{data}
\title{Complete MetaCyc demonstration contributions}
\description{The filtered MicroFun demonstration derived from Riquelme et al. (2019). Sample IDs are consistently renamed S001 to S068. Genus Incertae sedis ASVs were removed in the source release.}
\usage{microfun_demo_MetaCyc}
\format{A data frame with 1,591,512 rows and 8 columns: \code{sample}, \code{function}, \code{taxon}, \code{taxon_abun}, \code{taxon_rel_abun}, \code{genome_function_count}, \code{taxon_function_abun}, \code{taxon_rel_function_abun}. All 68 samples and 2317 contributing ASVs are represented. Identifier columns are character; abundance columns are numeric. KO IDs retain their source \code{ko:} prefix. Additional normalized columns are upstream values and were not recalculated after filtering.}
\details{The MDA cohort contains 22 STS and 21 LTS samples; JHU contains 10 STS and 15 LTS samples. Species annotations are entirely missing. The KEGG tables preserve the 2026-11-04 export and differ from the bundled KO-to-KEGG mapping. Predictions are fractional functional estimates, not measured read counts. The original PICRUSt2 version and database are not available. Confirm data redistribution terms before publishing; the software MIT license does not establish rights to all data. Full tables consume more memory than their compressed package files. Load namespaces independently.}
\source{\url{https://github.com/bioscinema/MicroFun_website_tutorials}, commit 998e0473ec492e1a5f155dd67a419354e12e8558, assets/downloads/study filtered data archives. Riquelme et al. (2019), Cell 178(4):795-806.e12, \doi{10.1016/j.cell.2019.07.008}.}
\author{Ruitao Liu, Chenyu Liu, Yaoran Yu, Zilin Li and Liangliang Zhang}
\seealso{\code{\link{pathway_sdaa}}}
\examples{
data("microfun_demo_MetaCyc", package = "MicroFun")
head(microfun_demo_MetaCyc)
}
\keyword{datasets}
