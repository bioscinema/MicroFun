MicroFun demonstration: loading and provenance

Use data("microfun_demo", package="MicroFun") for the core list.
Use data("microfun_demo_KO", package="MicroFun") for full KO contributions.
EC, KEGG and MetaCyc have independent datasets named microfun_demo_EC,
microfun_demo_KEGG and microfun_demo_MetaCyc.
All contribution datasets are ordinary data frames, retaining all source columns.
Pass one table using pathway_sdaa(stratified_data=table, physeq=ps, ...).

Source study: Riquelme E, Zhang Y, Zhang L, et al. (2019).
Tumor Microbiome Diversity and Composition Influence Pancreatic Cancer Outcomes.
Cell 178(4):795-806.e12. https://doi.org/10.1016/j.cell.2019.07.008
PICRUSt2 method: Douglas GM et al. (2020), Nature Biotechnology 38:685-688.
https://doi.org/10.1038/s41587-020-0548-6
Phyloseq: McMurdie PJ and Holmes S (2013), PLOS ONE 8(4):e61217.
https://doi.org/10.1371/journal.pone.0061217

68 renamed sample IDs S001-S068; MDA 22 STS / 21 LTS; JHU 10 STS / 15 LTS.
2,487 retained ASVs after removal of 290 Genus Incertae sedis ASVs.
2,317 ASVs occur in each contribution table; no reason for the other 170
ASVs' absence can be established from these exports alone.
Genus is unassigned for 725 retained ASVs; all Species values are unassigned.
Do not interpret missing contribution rows as measured absence of a function.
Predicted contributions can be fractional; they are not shotgun read counts.
Only sample, function, taxon and taxon_function_abun are needed by this workflow.
Additional relative/normalized fields are preserved from upstream and were not
renormalized after ASV filtering. Prefer raw taxon_function_abun for analysis.

Standard KO / EC / KEGG / MetaCyc dimensions: 8881 / 2898 / 234 / 542 by 68.
The dated KEGG_2026-11-04 export has different mapping coverage than the
package's KO-to-KEGG reference. That underlying dated mapping is not supplied.
The original PICRUSt2 version, reference database and full prediction commands
are not available in these processed exports. Raw reads are not bundled.
Data are for demonstration and do not establish a validated clinical signature.

Core fields: ps, metadata, otu, taxonomy, phylo_tree, standard, selections,
KO_display_labels, filter_summary, removed_ASVs, source_manifest, provenance.
The source_manifest lists original TSV.GZ checksums, not RDA checksums.
The package's data-raw directory contains a rebuilding script and source details.

Dataset authors: Ruitao Liu; Chenyu Liu; Yaoran Yu; Zilin Li; Liangliang Zhang.
Data rights are not established by the software MIT license. Confirm original
study data and third-party annotation redistribution terms before publishing.
