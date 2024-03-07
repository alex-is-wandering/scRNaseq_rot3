SCRIPTS.

1. FOLDER: converting H5ADtoSeurat
	The first part is a python script, requiring h5ad-type data, outputting filesl which go into R - the second part of the script. "convertingH5AD" was a brief attempt at using seurat-disk to just read the H5AD directly - it didn't work as well. 
	Notes are included in the R script re: what most lines should be doing.
	
2. "cosgrove2020_Seuratconstruction.R"
	Sourced from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143704), this pulls together the raw data file, the metadata.txt.gz, and more metadata (sex, age, sample location) from the series matrix file. 
	
3. "cosgrove2020_skeletalmuscle.R"
	IN PROGRESS.
	ECM component synthesis/Matrisome association/etc of the above data with the ECM. So far, only contains command for QC plots (features, counts, mito%)
	
4. "ECMprocessing_heatmaps_dotplots.R"
	A compilation of assorted snippets, not particular to one dataset  - everything in it should also be more organized and contained within other files, but I haven't cleaned it up piece by piece yet.
	
5. FOLDER: elmentaite2021_gutdatatesting
	Contains 2 files: "gutatlas_rawdata" and "gutatlas_comparisons". "rawdata" pulls together the downloaded data from https://www.gutcellatlas.org/ and runs compares gene expression to Matrisome categories etc, after QC and processed/raw files were compared in "_comparisons". However, due to odd annotations we're leaving that one be.
	
6. "ReactomePathways_reformatting.R"
	ReactomePathways (this pathway specifically: https://maayanlab.cloud/Harmonizome/gene_set/Extracellular+matrix+organization/Reactome+Pathways) had text formatted oddly, and this script applies packages purrr and jsonlite to translate that into a functional list.

7. FOLDER: subramanian2021_ECMprocessing
	I started out making individual scripts for each Seurat object made from cell-type subsets from the original whole-kidney data; this includes dotplots, violins, addmodules scores for pathways (ECM organization) and the Matrisome Project categories. Pulling that together into a more organized script ("Compiled_ECMrelated_kidneyanalysis").
	"ReactomePathway_ECMOrganization_processing" is similarly mostly redundant due to the section in "Compiled" about module scoring against ECMOrg and ECM Receptor Interaction. 
	
8. "tutorialdataprocessing.R"
	From the Seurat tutorial (), following commands on a data set that had already been QC'd - no results here, disregard. 
	
9. "filtering_scaledata_bycluster.R"
	Filtering the Seurat objects lakemenon2023 and subramanian2021 for genes of relevance, i.e. 30% of cells >=2 expression in scale.data for at least one cluster, for a gene to be selected. Produces a table with "gene", "matrisome category", "cell types that meet threshold"
