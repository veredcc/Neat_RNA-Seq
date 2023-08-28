## Global settings for Neat_RNA-Seq script

# About

Update this file with the parameters and input files of your current project, before executing "Neat_RNA-Seq_1.0.Rmd".

# Global parameters
  
```{r global params}

## General

set.seed(111)

project_name = "HFD_Ad"

analysis_round = "01a_Midaged_devel"  # options, for example:  Kidney, Heart, or any name for a round of an analysis        

## Gene annotation

#modify gene ID?
CHOP_GENE_ID_BY_DELIMITER = TRUE #FALSE: to leave gene ID as is, e.g. in trinity results.
#TRUE: to remove the part after the delimiter, e.g. when gene IDs have the form e.g. ENSMUSG00000027351_Spred1
GENE_ID_DELIMITER = "_"   

#annotation source
#required to modify annotation column headers and edit their content
#first column must be the key for joining to the RSEM counts data frame (e.g. to contain Ensembl ID in the same format as in RSEM results)
#currently only implemented for ensembl annotation downloaded from Ensembl/biomart which includes Gene stable ID, Gene name, Gene type, Gene description
ANNOT_SOURCE = "Ensembl" #options: Ensembl, Trinotate (currently Trinotate does not have any effect)

## Statistical analysis for differential expression

#set design formula for DESeq2

CALC_CONTRASTS   = TRUE   #whether or not to calculate contrasts (TURE, FALSE)
CALC_INTERACTION = FALSE   #whether or not to calculate interaction (TURE, FALSE)

#set DESIGN for contrasts

DESIGN = "~ Diet_per_age"  #examples: "~ Stage", "~ Batch + Diet"

#set DESIGN for interaction

DESIGN_INTERACTION = "~ Age + Diet + Age:Diet"
LRT                = "~ Age + Diet"

#use expanded or simple model in DESeq function (default - expanded)
#simple model is when the control of all treatment groups is in the first treatment group
#for more complex designs, use expanded model

use_expanded_model = FALSE

#set cutoffs for DE contrasts

LINEAR_FC_CUTOFF  = 1     #e.g. 1.3
PADJ_CUTOFF       = 0.05  #e.g. 0.05
#PVAL_CUTOFF = 0.01       #pay attention, if you use this variable, you have to change the code in the "compute_contrasts" function in the functions file !!!
DESEQ_PADJ_CUTOFF = 0.05  #e.g. 0.05

#set cutoffs for DE interaction

PADJ_CUTOFF_INTERACTION       = 0.1
#PVAL_CUTOFF_INTERACTION      = 0.01       #pay attention, if you use this variable, you have to change the code in the "compute_contrasts" function in the functions file !!!
DESEQ_PADJ_CUTOFF_INTERACTION = 0.1

## Visualization

#set factors for PCA, heatmap, and excel top rows
#specify col names from col_data, starting from the factor with the biggest effect

EFFECTS = c("Batch", "Diet_per_age")

#normalization method (for visualization)
NORM_METHOD = 'VSD'  #'VSD' or 'RLOG'

#batch correction

REMOVE_BATCH_EFFECT = F
BATCH_EFFECTS = c('Batch')  #batch effects as columns of col_data. sva can use only one. limma can use up to two.
BATCH_CORR_METHOD = 'sva' #sva or limma

## Clustering

#set group factor for the binary pattern calculations and for partition clustering
#(this is usually the factor that was used for the contrasts)

GROUP = "Diet_per_age"

#set correlation cutoff for binary pattern calculation

CORR_CUTOFF = 0.8

#set min counts cutoff for 1s in the binary pattern calculation

COUNTS_CUTOFF = 500 

#show only top n DE genes in hierarchical clustering?

SHOW_TOP_GENES_IN_HEATMAP = FALSE   #TRUE: show only top DE genes in heatmap. FALSE: show all DE genes in heatmap
NR_TOP_GENES = 5000

#partition clustering

K_FIXED = 12 #no. of requested clusters. if NA, K_MAX will be used
K_MAX = 20

GROUP1 = "Diet_per_age"   #usually GROUP will be the main treatment (e.g. diet), and GROUP1 will be another biological factor, e.g. age

#manual clustering
PERFORM_MANUAL_CLUSTERING = T  #TRUE or FALSE

## Enrichment analysis

FILTER_KEGG_PATHWAYS_BY_TAXON = NA  #provide KEGG taxon (name or ID?) to filter the pathways. for no filtering use NA
PLOT_GO_TERM_OVERLAPS = F  #plot heatmaps and create excels of gene2term and term2term overlaps

#passed to clusterProfiler::enrichr
ENRICHMENT_PVAL_CUTOFF = 0.05
ENRICHMENT_PADJ_METHOD = 'fdr' #options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#passed to enrichplot::dotplot as argument showCategory
MAX_TERMS_IN_DOTPLOT = 30  #number of enriched terms to display in dotplot

```

# Input files and directories

```{r infiles}

#the three following should be prepared by the user, and are analysis-specific

rsem_files_locations     = "Count_file_location_all_wo_Young_H8.txt"

experiment_design_file   = "Experiment_design_all_wo_Young_H8.txt"

contrasts_file           = "Contrasts_midaged.txt"

#the following files are per organism (or per de novo transcriptome assembly in non-model organisms)

#external annotation file (for organisms with a known genome, that have data in public databases such as ensembl)
#e.g. download from Ensembl/BioMart: Gene stable ID, Gene name, Gene type, Gene description (by default the file is called "mart_export.txt")
#or, use the data retrieved from AnnotationHub using the Retrieve_annotation.Rmd script ("Func_annot_data/Annotation.tab"), but this has not been tested yet
#the annotation data will appear as first cols in the Excel results file
annotation_file          = "mart_export.txt"  #or 

#file prepared by trinotate (for non-model organisms)
trinotate_file           = "Embryos_assembly_Tom.trino_anno_rep.xls"

#directory with data required for enrichment analyses
#kegg data are taken from Vered Perl script.
#GO data and the Annotation.tab file are from the R script.
#The perl script is located on Vered’s PC, at
#Programming/Perl_workspace/Assaf_Rudich_03_HFD_KEGG_from_Vered_script/KEGG_2_ensembl_conversion8.pl
functional_annot_dir     = "Func_annot_data"  



```