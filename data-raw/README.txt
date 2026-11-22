MicroFun bundled demonstration data build

Source: https://github.com/bioscinema/MicroFun_website_tutorials
Source commit: 998e0473ec492e1a5f155dd67a419354e12e8558
Source archives: assets/downloads/study/MicroFun_filtered_core.zip,
MicroFun_filtered_KO_part01.zip, MicroFun_filtered_KO_part02.zip,
MicroFun_filtered_EC_stratified.zip,
MicroFun_filtered_KEGG_2026-11-04_stratified.zip,
MicroFun_filtered_MetaCyc_stratified.zip.
Also use assets/downloads/study/MicroFun_selected_KOs.csv.
Extract all six ZIPs into one directory and run the original assemble-data.R
with that directory as working directory to join the two KO fragments.
Copy MicroFun_selected_KOs.csv into the assembled directory.

From the package root run:
Rscript data-raw/build-demo-data.R /absolute/path/to/MicroFun_demonstration_data_filtered /absolute/path/to/MicroFun

Build dependencies: R, phyloseq, data.table and R.utils (for gzip input).
Data are serialized as ordinary data frames and a core list with xz compression.
No original study values or row orders are changed during packaging.
The source MD5 and byte sizes are verified before serialization.
The removed-ASV list is provenance and must not be reintroduced to analyses.

Dataset authors, in order:
Ruitao Liu; Chenyu Liu; Yaoran Yu; Zilin Li; Liangliang Zhang.
Data license: not established by supplied source files; confirm before redistribution.
