## Functions for Vered_RNA-Seq_1.0 script
## by Vered Chalifa-Caspi

##### Preparations #####

import_rsem = function (rsem_files_locations, CHOP_GENE_ID_BY_DELIMITER=F, GENE_ID_DELIMITER=NA) {
  
  #read file location info
  file_table        = read.delim(rsem_files_locations, stringsAsFactors = F)
  file_table$Sample = make.names(file_table$Sample) #make the sample names Syntactically Valid names (it adds X before 0)
  
  #create vector of file names with sample names as element names
  file_vec          = file_table$File
  names(file_vec)   = file_table$Sample
  
  #tximport - import RSEM files
  txi.rsem <- tximport(file_vec, 
                       type = "rsem")
  
  head(txi.rsem$counts) 

  if (CHOP_GENE_ID_BY_DELIMITER) {
    rownames(txi.rsem$abundance)=lapply(X =rownames(txi.rsem$abundance),FUN = function(x)     unlist(stringi::stri_split(str = x,regex = GENE_ID_DELIMITER))[1] )
    rownames(txi.rsem$counts)=lapply(X =rownames(txi.rsem$counts),FUN = function(x) unlist(stringi::stri_split(str = x,regex = GENE_ID_DELIMITER))[1] )
    rownames(txi.rsem$length)=lapply(X =rownames(txi.rsem$length),FUN = function(x) unlist(stringi::stri_split(str = x,regex = GENE_ID_DELIMITER))[1] )
  }
    
  return (txi.rsem)
}

prepare_col_data = function (experiment_design_file) {

  #read column data from file
  col_data = read.delim(experiment_design_file,
                        header       = T,
                        comment.char = "#")

  #make sure the first column header is SampleID
  colnames(col_data)[1] <- 'SampleID'
    
  #make the sample names Syntactically Valid names (it adds X before 0)
  col_data$SampleID = col_data$SampleID %>%  
    as.character %>% 
    make.names()
  
  #re-factoring all col_data factors
  col_data = data.frame(lapply(X = col_data, FUN = factor))
  
  #set row names to SampleID
  rownames(col_data) <- col_data$SampleID

  #sampleID order in col_data
  identical(rownames(col_data),   
            colnames(txi.rsem$counts))   
  
  return (col_data)

}

remove_samples_from_coldata = function (col_data, sample_list) {
    col_data[rownames(col_data) %in% setdiff(col_data$SampleID, sample_list),]
}

filter_coldata_by_tissue = function (col_data, col_attribute, tissue) {
  col_attribute <- rlang::sym(col_attribute)
  filter (col_data, !!col_attribute == tissue)
}

get_replicate_struct_from_coldata = function (col_data, GROUP) {
  #get a vector of the number of replicates in each treatment group, in the order that they appear in col_data
  #GROUP specifies the column name of the treatment group in col_data
  col_data[[GROUP]] %>% fct_inorder %>% table %>% as.numeric
}

match_txi_rsem_to_coldata = function (txi.rsem, col_data) {
  
  SampleIDs = as.character(col_data$SampleID)
  for(tbl in c("abundance", "counts", "length")) {
    txi.rsem[[tbl]] <- txi.rsem[[tbl]][,SampleIDs]
  }
  
  # Check col_data and txi.rsem agree in column headers:
  for(tbl in c("abundance", "counts", "length")) {
    print(tbl)
    identical(rownames(col_data),   # SampleID order in col_data
              colnames(txi.rsem[[tbl]])) %>% print
  }
  
  return (txi.rsem)
}

process_annotation_data = function (annotation_file, ANNOT_SOURCE="") {
  
  ## Read annotation data
  
  annot = read.delim(file = annotation_file,
                     stringsAsFactors = F)
  
  #process annotation data

  if (ANNOT_SOURCE == "Ensembl") {
    #define new column names (this method avoids error in case column(s) do not exist)
    lookup <- c(Ensembl_ID = "Gene.stable.ID", Gene_Symbol = "Gene.name", Gene_Type = "Gene.type", Gene_Title = "Gene.description")
    #process
    annot = 
      annot %>%
      rename(any_of(lookup)) %>%                                   #rename columns according to lookup
      dplyr::select (Ensembl_ID, Gene_Symbol, Gene_Title, Gene_Type) %>%  #reorder columns
      mutate(Gene_Title = str_replace(string = Gene_Title,         #edit the gene description column: remove everything after " ["
                                      pattern = " \\[.*",
                                      replacement = ""))   
   
  }
  return (annot)
  
}

process_trinotate_data = function (trinotate_file) {
  
  ## Read Trinotate data
  
  trinotate <- read.delim(file = trinotate_file,
                          colClasses = c(X.gene_id="factor",
                                         transcript_id="factor",
                                         sprot_Top_BLASTX_hit="character",
                                         infernal="character", 
                                         prot_id="character", 
                                         prot_coords="character", 
                                         sprot_Top_BLASTP_hit="character",
                                         Pfam="character", 
                                         SignalP="character", 
                                         TmHMM="character", 
                                         eggnog="character",
                                         Kegg="character", 
                                         gene_ontology_BLASTX="character",
                                         gene_ontology_BLASTP="character",
                                         gene_ontology_Pfam="character", 
                                         transcript="character", 
                                         peptide="character"))
  
  ## Parse Trinotate data
  
  #get all gene IDs (to be used on the final merge)
  all_genes = trinotate %>%
    select (X.gene_id) %>%
    unique
  
  #get RNA predictions from infernal
  trinotate_main_infernal =
    trinotate[,c("X.gene_id",
                 "infernal")] %>% 
    as_tibble %>% 
    # Remove unidetified
    filter(infernal!="." ) %>% 
    unique %>%
    mutate(RNA_Infernal = str_replace(string = infernal,
                                      pattern = "`",
                                      replacement = "|")) %>%
    select ("X.gene_id", "RNA_Infernal")
  
  #get sprot_BLASTX results  
  trinotate_main_sprot_BLASTX <-
    trinotate[,c("X.gene_id",
                 "sprot_Top_BLASTX_hit")] %>% 
    as_tibble %>% 
    # Remove unidetified
    filter( sprot_Top_BLASTX_hit!=".") %>% 
    # For sprot blastx:
    # Hits are encoded in a table where field separator is '^' and record separator is '`' (back tick)
    # Keep only first hit
    mutate(sprot_edit = str_replace(string = sprot_Top_BLASTX_hit,
                                    pattern = "\\`.*",
                                    replacement = "")) %>% 
    # Remove curly braces
    mutate(sprot_edit = str_replace(string = sprot_edit,
                                    pattern = "\\{.*\\}",
                                    replacement = "")) %>%  
    select(X.gene_id,sprot_edit) %>% 
    # Separate 'sprot' data into columns
    separate(sprot_edit,
             into = paste("sprot",
                          c("Name","Acc","Pos","percID","eval","RecName","Lineage"),
                          sep="_"),
             sep = "\\^") %>%  
    # Remove RecName extra characters
    mutate(sprot_RecName = str_replace(string = sprot_RecName,
                                       pattern = "^RecName: Full=",
                                       replacement = "")) %>% 
    unique %>% 
    group_by(X.gene_id) %>% 
    summarise(sprot_Name=paste(sprot_Name,collapse="|"),
              sprot_Acc=paste(sprot_Acc,collapse="|"),
              sprot_Pos=paste(sprot_Pos,collapse="|"),
              sprot_percID=paste(sprot_percID,collapse="|"),
              sprot_eval=paste(sprot_eval,collapse="|"),
              sprot_RecName=paste(sprot_RecName,collapse="|"),
              sprot_Lineage=paste(sprot_Lineage,collapse="|")) %>%
    ungroup() %>% 
    # Select required columns
    select("X.gene_id","sprot_Name","sprot_RecName","sprot_percID","sprot_eval") %>% 
    unique
  
  # Doing the same for the protein BLASTP results
  trinotate_main_sprot_BLASTP <-
    trinotate[,c("X.gene_id",
                 "sprot_Top_BLASTP_hit")] %>% 
    as_tibble %>% 
    # Remove unidetified
    filter( sprot_Top_BLASTP_hit!="." ) %>% 
    # For sprot blastx
    # Keep only first hit
    mutate(sprot_edit = str_replace(string = sprot_Top_BLASTP_hit,
                                    pattern = "\\`.*",
                                    replacement = "")) %>% 
    # Remove curly braces
    mutate(sprot_edit = str_replace(string = sprot_edit,
                                    pattern = "\\{.*\\}",
                                    replacement = "")) %>%  
    select(X.gene_id,sprot_edit) %>% 
    # Separate into columns
    separate(sprot_edit,
             into = paste("sprot",
                          c("Name","Acc","Pos","percID","eval","RecName","Lineage"),
                          sep="_"),
             sep = "\\^") %>%  
    # Remove RecName extra characters
    mutate(sprot_RecName = str_replace(string = sprot_RecName,
                                       pattern = "^RecName: Full=",
                                       replacement = "")) %>% 
    # Select required columns
    select("X.gene_id","sprot_Name","sprot_RecName","sprot_percID","sprot_eval") %>%
    # Agglomerating multiple values per record
    group_by(X.gene_id) %>%
    summarise(blastp_Name = paste0(sprot_Name, collapse = "|"),
              blastp_RecName = paste0(sprot_RecName, collapse = "|"),
              blastp_percID = paste0(sprot_percID, collapse = "|"), 
              blastp_eval = paste0(sprot_eval, collapse = "|")) %>% 
    unique
  
  
  trinotate_main_pfam <-
    trinotate[,c("X.gene_id",
                 "Pfam")] %>% 
    filter(Pfam != ".") %>%  # Remove unidetified
    unique

  trinotate_main_pfam$pfam_domains = sapply(X = trinotate_main_pfam[,"Pfam"],
                                            FUN = parse_pfam)

  trinotate_main_pfam1 = 
    trinotate_main_pfam %>%
    select ("X.gene_id", "pfam_domains") %>%
    group_by(X.gene_id) %>% 
    summarise(pfam_all = paste(pfam_domains,collapse="|"))
  
  #merge all tables
  
  trinotate_main <- 
    all_genes %>%
    left_join(trinotate_main_sprot_BLASTX,by="X.gene_id") %>%
    left_join(trinotate_main_infernal,by="X.gene_id") %>%
    left_join(trinotate_main_sprot_BLASTP,by="X.gene_id") %>%
    left_join(trinotate_main_pfam1,by="X.gene_id") %>% 
    unique
  
  rm (trinotate_main_sprot_BLASTX)
  rm(trinotate_main_infernal)
  rm(trinotate_main_sprot_BLASTP)
  rm(trinotate_main_pfam)
  rm(trinotate_main_pfam1)
  
  gc()
  save.image()
  
  ## Statistical analysis of trinotate data:
  
  
  cat(sprintf("Number of genes for which annotation exists: %s\n",
              trinotate_main$X.gene_id  %>% unique %>% length))
  
  for (coln in  c( "sprot_Name", "sprot_RecName", "sprot_percID", "sprot_eval", "RNA_Infernal", "blastp_Name", "blastp_RecName", "blastp_percID", "blastp_eval", "pfam_all")) {
    cat(sprintf("Number of genes positive in column %s: %s\n",coln,
                trinotate_main[,coln] %>% is.na %>% not  %>% sum))
    
  }
  
  return (trinotate_main)
  
}

parse_pfam = function (pfam_data) {
  #example:
  #PF13540.9^RCC1_2^NULL^945-970^E:8e-07`PF00415.21^RCC1^NULL^957-1006^E:1.8e-14`PF13540.9^RCC1_2^NULL^993-1022^E:1.4e-10`PF00415.21^RCC1^NULL^1009-1058^E:1.4e-12`PF13540.9^RCC1_2^NULL^1051-1074^E:1.2e-06`PF00415.21^RCC1^NULL^1061-1108^E:2.7e-11`PF13540.9^RCC1_2^NULL^1351-1380^E:6e-08`PF00415.21^RCC1^NULL^1368-1414^E:1.1e-14`PF13540.9^RCC1_2^NULL^1401-1430^E:5.5e-09
  pfam_items = strsplit(pfam_data, "`", fixed=T) %>% unlist
  pfam_items1 = c()
  for (item in pfam_items) {
    elems = strsplit(item, "^", fixed=T) %>% unlist
    pfam_id = str_replace(string=elems[1], pattern="\\..*", replacement="")
    pfam_short_name = elems[2]
    pfam_join = paste(pfam_id, pfam_short_name)
    pfam_items1 = c(pfam_items1, pfam_join)
  }
  pfam_summary = pfam_items1 %>% unique %>% paste0(collapse=", ")
  return(pfam_summary)
}

##### Statistical analysis #####

run_deseq2 = function (txi.rsem, col_data, DESIGN, use_expanded_model=TRUE) {
  
  dds0 <- DESeqDataSetFromTximport(txi = txi.rsem,
                                   colData = col_data,
                                   design =  as.formula(DESIGN))
  
  if (use_expanded_model) {
    dds <- DESeq(dds0,
                 #parallel = TRUE,
                 betaPrior = TRUE,  
                 modelMatrixType = "expanded")
  } else {
    dds <- DESeq(dds0,
                 betaPrior = FALSE)
  }
  
  resultsNames(dds)
  
  return (dds)
}

run_deseq2_interaction = function (txi.rsem, col_data, DESIGN_INTERACTION, LRT, use_expanded_model=TRUE) {
  
  dds_interaction0 <- DESeqDataSetFromTximport(txi = txi.rsem,
                                               colData = col_data,
                                               design =  as.formula(DESIGN_INTERACTION))
  
  if (use_expanded_model) {
    dds_interaction <- DESeq(dds_interaction0,             #I have not tested that with LRT (26.3.2023)
                       test="LRT",
                       reduced = as.formula(LRT),
                       betaPrior = TRUE,  
                       modelMatrixType = "expanded")
  } else {
    dds_interaction <- DESeq(dds_interaction0,
                       test="LRT",
                       reduced = as.formula(LRT),
                       betaPrior = FALSE)
  }
  
  resultsNames(dds_interaction)
  
  return (dds_interaction)
}

compute_contrasts = function (dds, contrasts_data) {

  #initialize stats_df data frame with gene name as first column
  stats_df = data.frame(gene = rownames(rowData(dds)))
  
  #compute stats for each contrast and add to stats_df data frame
  for (i in 1:nrow(contrasts_data)) {
    contrast = c(contrasts_data[i,'Factor'], contrasts_data[i,'Numerator'], contrasts_data[i,'Denominator'])
    contrast_name = contrasts_data[i,'Contrast_name']
    print (contrast_name)
    stats_df = 
      results(dds,
              contrast=contrast,
              alpha = DESEQ_PADJ_CUTOFF) %>%
      as.data.frame %>% 
      as_tibble %>% 
      mutate(linearFC = ifelse(is.na(log2FoldChange),
                               yes = NA,
                               no = ifelse(log2FoldChange>0,
                                           yes = 2^log2FoldChange,
                                           no = -1/(2^log2FoldChange))%>% 
                                 signif(digits = 3))) %>% 
      mutate(pvalue = ifelse(test = is.na(pvalue),
                             yes = NA,
                             no = signif(pvalue,digits = 3))) %>% 
      mutate(padj = ifelse(test = is.na(padj),
                           yes = NA,
                           no = signif(padj,digits = 3))) %>% 
      mutate(pass = ifelse(test = abs(as.numeric(linearFC)) >= LINEAR_FC_CUTOFF & 
                             padj <= PADJ_CUTOFF & #uncomment to use FDR cutoff !!! 
                             #pvalue <= PVAL_CUTOFF &  #uncomment to use unadjusted p-value cutoff !!!            
                             !is.na(padj),
                           yes = ifelse(test = as.numeric(linearFC)>0,
                                        yes="up",
                                        no="down"),
                           no = "")) %>%
      mutate(manual_cutoffs = manual_cutoff_formula) %>%
      select(linearFC,pvalue,padj,pass,manual_cutoffs) %>%
      # select(linearFC,pvalue,padj,pass) %>%                #without manual_cutoffs formula cols
      dplyr::rename(!!paste0("linearFC.",      contrast_name) := linearFC) %>%
      dplyr::rename(!!paste0("pvalue.",        contrast_name) := pvalue) %>%
      dplyr::rename(!!paste0("padj.",          contrast_name) := padj) %>%
      dplyr::rename(!!paste0("pass.",          contrast_name) := pass) %>%
      dplyr::rename(!!paste0("manual_cutoffs.",contrast_name) := manual_cutoffs) %>%
      as.data.frame %>% 
      cbind(stats_df,.)
  }
  
  # Add pass any field
  
  if (sum(str_detect(string = names(stats_df),
                     pattern = "pass")) >1 ) {
    stats_df$pass_any <-
      apply(X = stats_df[,str_detect(string = names(stats_df),
                                     pattern = "pass")],
            MARGIN = 1,
            FUN = function(x) any(unlist(x)!=""))
    # Convert TRUE/FALSE to 1/""
    stats_df$pass_any <- ifelse(test = !is.na(stats_df$pass_any) & stats_df$pass_any==TRUE,
                                yes = 1,
                                no = "")
  } else {
    stats_df$pass_any <- ifelse(test = stats_df[,str_detect(string = names(stats_df),pattern = "pass")]!="",
                                yes = 1,
                                no = "")
  }
  
  names(stats_df) <- make.names(names(stats_df),
                                unique = T) 
  
  return (stats_df)
}

compute_interaction = function (dds_interaction) {
  
  stats_df_interaction0 = results(dds_interaction, alpha = DESEQ_PADJ_CUTOFF_INTERACTION)

  #initialize stats_df_interaction data frame with gene name as first column
  stats_df_interaction = data.frame(gene = rownames(rowData(dds_interaction)))
    
  stats_df_interaction = 
    stats_df_interaction0 %>%
    as.data.frame %>% 
    as_tibble %>% 
    mutate(pvalue.interaction = ifelse(test = is.na(pvalue),
                           yes = NA,
                           no = signif(pvalue,digits = 3))) %>% 
    mutate(padj.interaction = ifelse(test = is.na(padj),
                         yes = NA,
                         no = signif(padj,digits = 3))) %>% 
    mutate(pass.interaction = ifelse(test = 
                           padj <= PADJ_CUTOFF_INTERACTION & #uncomment to use FDR cutoff !!! 
                           #pvalue <= PVAL_CUTOFF &          #uncomment to use unadjusted p-value cutoff !!!            
                           !is.na(padj),
                         yes = "interaction",
                         no = "")) %>%
    select(pvalue.interaction,padj.interaction,pass.interaction) %>%
    as.data.frame %>%
    cbind(stats_df_interaction,.)
  
  names(stats_df_interaction) <- make.names(names(stats_df_interaction),
                                unique = T) 
  
  return (stats_df_interaction)
}

combine_stat_dfs = function (stats_df, stats_df_interaction) {
  stats_df_combined = cbind(stats_df, select(stats_df_interaction, -gene))
  
  stats_df_combined$pass_combined <- ifelse(
    test = stats_df_combined$pass_any=="1" | stats_df_combined$pass.interaction=="interaction",
    yes = 1,
    no = "")
  return (stats_df_combined)
}

de_summary_stats = function (stats_df, DE_genes_stats_file) {
  
  DE_genes_stats =
    apply(X = stats_df[,str_detect(string = names(stats_df),
                                   pattern = "pass")],
          MARGIN = 2,
          FUN = function(x) c(up=sum(x=="up"),down=sum(x=="down"),any=sum(x!=""))) %>%
    data.frame %>%
    t
  
  export_table (DE_genes_stats, DE_genes_stats_file)
  
  return (DE_genes_stats)
}

##### Filtering #####

get_DE_genes_list = function (stats_df) {
  DE_genes <-
    stats_df %>%
    filter(pass_any=="1") %>% 
    dplyr::select(gene) %>%
    unlist
  return (DE_genes)
}

get_DE_genes_list_interaction = function (stats_df_interaction) {
  DE_genes_interaction <-
    stats_df_interaction %>%
    filter(pass.interaction=="interaction") %>% 
    dplyr::select(gene) %>%
    unlist
  return (DE_genes_interaction)
}

get_DE_genes_list_combined = function (stats_df_combined) {
  DE_genes_combined <-
    stats_df_combined %>%
    filter(pass_combined=="1") %>% 
    dplyr::select(gene) %>%
    unlist
  return (DE_genes_combined)
}

get_top_DE_genes = function (stats_df, n=2000) {

  #get top n DE genes based on min padj in any comparison, sorted by this padj
  
  padj_cols <- str_detect(string = names(stats_df),
                          pattern = "padj")
  
  top_DE_genes <-
    stats_df %>%
    # Add column with minimum padj
    cbind(min_padj=apply(X = stats_df[,padj_cols] ,
                         MARGIN=1,
                         FUN = function(x) ifelse(test = all(is.na(x)),
                                                  yes = NA,
                                                  no = min(x,na.rm = T)))) %>%
    as_tibble %>% 
    filter(pass_any=="1") %>% 
    # Select top n genes by descending padj
    top_n(n = n,wt = desc(min_padj)) %>% 
    dplyr::select(gene) %>%
    unlist

  #if, for example, the n+1 gene has the same min padj as the n gene, the top_n function will return n+1 genes.
  #therefore, I add the code below to make sure that only n genes will be returned.
  
  if(length(top_DE_genes)>n) {
    top_DE_genes = top_DE_genes[1:n]
  }
  return (top_DE_genes)
}

get_DE_genes_list_per_contrast = function (contrasts_data, stats_df) { #Vered 13.8.2023
  DE_genes_lists = list()
  for (contrast in contrasts_data$Contrast_name) {
    pass_col = paste0("pass.", contrast)
    up   = stats_df %>% filter (get({{pass_col}}) == "up") %>% pull (gene)
    down = stats_df %>% filter (get({{pass_col}}) == "down") %>% pull (gene)
    all = c(up, down)
    DE_genes_lists[[paste0(contrast, "_up")]]   = up
    DE_genes_lists[[paste0(contrast, "_down")]] = down
    DE_genes_lists[[paste0(contrast, "_all")]]  = all
  }
  return (DE_genes_lists)
}

get_genes_from_file = function (file_name, has_header=FALSE) {
  
  # get gene list from file, assuming genes are in the first column
  
  data = read.delim(file_name, as.is=T, header=has_header)
  gene_list = data[[1]]
  
  return(gene_list)
}

get_genes_with_corr_to_pattern = function (corrs, pattern, CORR_CUTOFF, sort_by_corr=T) {
  corr_col_name = paste0("_", pattern)
  corr_col_name <- rlang::sym(corr_col_name)
  corrs_f = corrs %>% rownames_to_column("gene") %>% filter (!!corr_col_name > CORR_CUTOFF)
  if (sort_by_corr == TRUE) {
    corrs_f = arrange(corrs_f, desc(!!corr_col_name))
  }
  return (corrs_f$gene)
}

get_genes_with_corr_to_best_pattern = function (corrs2, pattern, sort_by_corr=T) {
  pattern_name = paste0("_", pattern)
  corrs2_f = corrs2 %>% rownames_to_column("gene") %>% filter (best_pattern == pattern_name)
  if (sort_by_corr == TRUE) {
    pattern_name <- rlang::sym(pattern_name)
    corrs2_f = arrange(corrs2_f, desc(!!pattern_name))
  }
  return(corrs2_f$gene)
}

filter_expression_matrix_by_gene_list = function (expr_data, gene_list) {
  #assuming row names are gene names
  genes <- rownames(expr_data) %in% gene_list
  return (expr_data[genes,])
}

##### Binary patterns #####

zero_one_sequences = function(lenght){
  return(zero_one_sequences_helper(lenght, "", c()))
}

zero_one_sequences_helper = function(lenght, sequence, vec){
  if (nchar(sequence) == lenght){
    vec <- c(vec, sequence)
    return(vec)
  }
  
  vec <- zero_one_sequences_helper(lenght, paste(sequence, "0", sep = "", collapse = ""), vec)
  vec <- zero_one_sequences_helper(lenght, paste(sequence, "1", sep = "", collapse = ""), vec)
  return(vec)
}

get_all_possible_binary_patterns = function (pattern_length) {
  patterns <- zero_one_sequences(pattern_length)
  cat ("No. of binary patterns: ", length(patterns), " (", patterns[1:4], "... )") 
  return (patterns)
}

calc_pattern_1s_pass_min_counts = function (expr_data, patterns, replicates_structure, count_cutoff) {
  
  #create a matrix basen on the norm counts.
  #every gene is marked "TRUE" for a pattern if all the samples that have to be 1 according to the pattern, have counts above a cutoff
  
  #count_cutoff = 500           #for testing - to run the function directly
  #expr_data = norm_counts_DE   #for testing - to run the function directly
  
  #create matrix with all results
  #in the matrix, an _ is appended to pattern name, otherwise it does not show well in Excel when exported
  
  pattern_length = nchar(patterns[1])
  
  #initialize results table by just taking row names
  passes = expr_data %>% as.data.frame %>% select(starts_with("999"))   
  
  calc_pass = function (x) {
    bool_vec = x > count_cutoff   #vector showing for each element of the expression pattern, is it above count cutoff
    indexes = which (pat1 %in% 1) #indexes of samples which has to be 1 according to the pattern
    all (bool_vec[indexes])       #if all the samples which has to be 1 have true in the vector that checks expression above cutoff
  }
  
  for (pattern in patterns) {
    # neglect patterns  "0000" and "1111"
    if (pattern == paste(rep(0, pattern_length), collapse ="") || pattern == paste(rep(1, pattern_length), collapse ="")) next
    # convert pattern to vector, e.g. 1, 0, 1, 1
    pat = str_split(pattern, "") %>% unlist %>% strtoi  
    # match to replicate structure, e.g. 1 1 1 0 0 0 1 1 1 1 1 1
    pat1 = rep(pat, replicates_structure)                               
    # apply
    a = apply(expr_data, 1, calc_pass) %>% as.data.frame
    passes_col_name = paste0("_", pattern) 
    names(a) = passes_col_name
    passes = cbind(passes, a)
  }
  
  return (passes)  
  
}

calc_correlations_to_expression_patterns = function (expr_data, patterns, replicates_structure) {

  #create matrix with all correlation results
  #in the matrix, an _ is appended to pattern name, otherwise it does not show well in Excel when exported
  
  pattern_length = nchar(patterns[1])
  
  #initialize results table by just taking row names
  corrs = expr_data %>% as.data.frame %>% select(starts_with("999"))   
  
  for (pattern in patterns) {
    # neglect patterns  "0000" and "1111"
    if (pattern == paste(rep(0, pattern_length), collapse ="") || pattern == paste(rep(1, pattern_length), collapse ="")) next
    # convert pattern to vector, e.g. 1, 0, 1, 1
    pat = str_split(pattern, "") %>% unlist %>% strtoi  
    # match to replicate structure, e.g. 1 1 1 0 0 0 1 1 1 1 1 1
    pat1 = rep(pat, replicates_structure)                               
    # apply
    a = apply(expr_data, 1, cor, pat1) %>% as.data.frame
    corr_col_name = paste0("_", pattern) 
    names(a) = corr_col_name
    corrs = cbind(corrs, a)
  }
  
  return (corrs)

}

calc_best_pattern_per_gene = function (corrs, CORR_CUTOFF=0.8, result_file_name, stats_file_name, passes) {

  #create matrix indicating only correlations above cutoff, no. of patterns above cutoff, and the best pattern
  #best pattern is the pattern with the max. correlation, if the correlation is above the cutoff
  
  #best pattern pass counts cutoff shows the best pattern if counts in all samples considered 1 in the pattern have more than a certain no. of counts,
  #based on the "passes" object. The "passes" object is supposed to contain the same genes in the same order as in corrs
  
  get_pattern = function (a) {
    if(max(a)>CORR_CUTOFF){
      col_no = which.max(a)
      col_name = names(a)[col_no]
      col_name1 = str_remove(col_name, "Corr_with")
      return (col_name1)
    } else {
      return (NA)
    }
  }
  
  check_pass = function (a) {
    best_p = as.character (a["best_pattern"])
    if (!is.na (best_p) && isTRUE(as.logical(a[best_p]))) {
      return (best_p)
    } else {
      return (NA)
    }
  }
  
  corrs1 = as.matrix(corrs)
  corrs2 = as.data.frame (ifelse(corrs1 > CORR_CUTOFF, corrs1, NA))
  corrs2$nr_patterns  = apply(corrs, 1, function (x) sum(x>CORR_CUTOFF))
  corrs2$best_pattern = apply(corrs, 1, get_pattern)
  corrs2$best_pattern_pass_counts_cutoff = apply(cbind(passes, best_pattern=corrs2$best_pattern), 1, check_pass)

  #print corrs2
  
  export_table(corrs2, result_file_name)
    
  #show and print stats

  corr2_stats = table(corrs2$best_pattern)
  cat ("No. genes per pattern")
  print (corr2_stats)
  
  corr2_stats_with_counts_cutoff = table(corrs2$best_pattern_pass_counts_cutoff)
  cat ("No. genes passing counts cutoff")
  print (corr2_stats_with_counts_cutoff)
  
  corr2_stats_merged = as.data.frame(cbind(corr2_stats, corr2_stats_with_counts_cutoff))
  
  corr2_stats_merged$Pattern = rownames(corr2_stats_merged)
  
  corr2_stats_merged = corr2_stats_merged[,c(3,1,2)]

  colnames(corr2_stats_merged) = c("Pattern", "No._genes", "No._genes_passing_counts_cutoff")
  
  write.table(x = corr2_stats_merged,
              file = stats_file_name,
              quote = F,
              sep = "\t",
              row.names = F)

  return (list(corrs2=corrs2, binary_patterns_stats=corr2_stats_merged))
}

##### Clustering #####

direction = function (FC, cutoff) {
  if (FC>cutoff) {
    return ("up")
  } else if (FC<-cutoff) {
    return ("down")
  } else {
    return ("nc")
  }
}

##### Graphics #####

draw_sample_correlation_matrix = function (norm_log_counts, results_dir) {

  result_file = sprintf("%s/sample_correlation_heatmap.png",results_dir)
  
  norm_log_counts = remove_first_x_from_all_colnames_if_exist(norm_log_counts)
  
  sampleDists <- dist(t(norm_log_counts))

  sampleDistMatrix <- as.matrix(sampleDists)

  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  png(filename = result_file)
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  dev.off()
  
}

draw_pca = function (norm_log_counts, results_dir, pc_x=1, pc_y=2, shape, color) {

  axis1 = paste0("PC", pc_x)
  axis2 = paste0("PC", pc_y)
  
  shape_factor <- shape
  color_factor <- color

  result_file = sprintf("%s/PCA_%s.vs.%s.png",results_dir,axis1,axis2)
    
  #perform PCA analysis
  norm_log_counts_pca <- prcomp(t(norm_log_counts))
  
  norm_log_counts_pca_percentVar <- norm_log_counts_pca$sdev^2/sum(norm_log_counts_pca$sdev^2)
  names(norm_log_counts_pca_percentVar) <- colnames(norm_log_counts_pca$x)
  
  pca2plot <- data.frame(norm_log_counts_pca$x[,c(axis1,axis2)],
                         Shape=col_data[,shape_factor],
                         Color=col_data[,color_factor],
                         Sample=col_data$SampleID)
  
  
  png(filename = result_file)
  
  plot = ggplot(pca2plot,aes(x=get(axis1),
                      y=get(axis2),
                      col=Color,
                      shape=Shape)) +
    geom_point(size=5) +
    # geom_label_repel(aes(label=Sample)) +
    # geom_path(aes(color=Type,
    #               group=Type)) +
    xlab(sprintf("%s: %1.2f%% of variance", axis1, norm_log_counts_pca_percentVar[axis1]*100)) + 
    ylab(sprintf("%s: %1.2f%% of variance", axis2, norm_log_counts_pca_percentVar[axis2]*100)) + 
    coord_fixed()
  print(plot)
  dev.off()
  
  return (norm_log_counts_pca)
}

draw_pca_3d = function (norm_log_counts_pca, pca_3d_file, shape, color) {

  shape_factor <- shape
  color_factor <- color
  
  result_file = pca_3d_file  #I did not manage to save it to a file within a directory

  norm_log_counts_pca_percentVar <- norm_log_counts_pca$sdev^2/sum(norm_log_counts_pca$sdev^2)
  names(norm_log_counts_pca_percentVar) <- colnames(norm_log_counts_pca$x)
      
  pca2plot3d <- data.frame(norm_log_counts_pca$x[,c("PC1","PC2","PC3")],
                           Color=col_data[,color_factor],
                           Shape=col_data[,shape_factor],
                           Sample=str_replace(string = col_data$SampleID,
                                              pattern = "^X",
                                              replacement = ""))
  
  symbols2use <- c('diamond','triangle-down','square','circle','x','o')
  
  p<-plot_ly(pca2plot3d) %>% 
    add_markers(x = ~PC1, y = ~PC2, z = ~PC3, 
                colors = rainbow(4),     #c("red","blue"),
                text=~Sample, 
                color = ~Color,
                symbol = ~Shape,
                symbols = symbols2use)  %>% 

    layout(scene = list(xaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC1",
                                                     norm_log_counts_pca_percentVar["PC1"]*100)),
                        yaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC2",
                                                     norm_log_counts_pca_percentVar["PC2"]*100)),
                        zaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC3",
                                                     norm_log_counts_pca_percentVar["PC3"]*100))))
  
  saveWidget(p, file = result_file, selfcontained = T, title=sprintf("Principal Component Analysis"))
}

draw_ma_and_volcano_plots = function (dds, contrasts_data, results_dir) {

  #compute stats for each contrast and draw plots
  for (i in 1:nrow(contrasts_data)) {
    contrast = c(contrasts_data[i,'Factor'], contrasts_data[i,'Numerator'], contrasts_data[i,'Denominator'])
    contrast_name = contrasts_data[i,'Contrast_name']
    print (contrast_name)

    compar_res = results(dds,
                          contrast=contrast,
                          alpha = DESEQ_PADJ_CUTOFF)
    
    MA_plot_file = sprintf("%s/MAPlot_%s.png",results_dir,contrast_name)
    volcano_plot_file = sprintf("%s/VolcanoPlot_%s.png",results_dir,contrast_name)
    
    # Draw MA plot
    
    png(filename = MA_plot_file)
    plotMA(compar_res,
           alpha = PADJ_CUTOFF,
           main=paste0("MAplot ",contrast_name),
           ylim = range(compar_res$log2FoldChange,na.rm = T)
    )
    abline(h = c(LOG_FC_CUTOFF,-LOG_FC_CUTOFF),col="green")
    dev.off()

    # Draw volcano plot
    # code adopted from: https://support.bioconductor.org/p/78962/

      data2plot <- compar_res %>% 
      as_tibble %>% 
      mutate(highFC=(abs(log2FoldChange)>LOG_FC_CUTOFF),
             lowPADJ=(padj<=PADJ_CUTOFF),
             signif=((abs(log2FoldChange)>LOG_FC_CUTOFF) & (padj<=PADJ_CUTOFF)),
             neglog10=-log10(pvalue)) 
    
    volcano_plot <- 
      ggplot(data2plot, aes(x=log2FoldChange, y=neglog10)) +
      geom_point(aes(colour = signif), 
                 size=2.5) +
      scale_colour_manual(values = c("TRUE"= "red", "FALSE"= "black")) + 
      xlab(expression(log[2]~fold~change)) +
      ylab(expression(-log[10]~pvalue)) +
      geom_vline(xintercept = c(-LOG_FC_CUTOFF, LOG_FC_CUTOFF),color = "blue", lty = 2) +
      ggtitle(paste0("Volcano plot ",contrast_name)) 
    ggsave(plot = volcano_plot, 
           filename = volcano_plot_file)
    
  }
  
}

draw_ma_and_volcano_plots_pdf = function (dds, contrasts_data, results_dir) {

  #this function works (for volcano), butthe PDF file it creates is too big and makes problems after opening it.
    
  pdf_file_name = sprintf("%s/MA_and_volcano_Plot.pdf",results_dir)
  p = list()  #list data structure to hold all graphs
   
  #compute stats for each contrast and draw plots
  for (i in 1:nrow(contrasts_data)) {
    contrast = c(contrasts_data[i,'Factor'], contrasts_data[i,'Numerator'], contrasts_data[i,'Denominator'])
    contrast_name = contrasts_data[i,'Contrast_name']
    print (contrast_name)
    
    compar_res = results(dds,
                         contrast=contrast,
                         alpha = DESEQ_PADJ_CUTOFF)

    # Draw MA plot:
    
    # plot_name = paste0 ("MA plot ", contrast_name)
    # 
    # p[[plot_name]] = plotMA(compar_res,
    #                         alpha = PADJ_CUTOFF,
    #                         main=paste0("MAplot ",contrast_name),
    #                         ylim = range(compar_res$log2FoldChange,na.rm = T)
    # )
    # #abline(h = c(LOG_FC_CUTOFF,-LOG_FC_CUTOFF),col="green")
    
    # Draw volcano plot:
    
    # plot_name = paste0 ("Volcano plot", contrast_name)
    plot_name = contrast_name
    
    data2plot <- compar_res %>% 
      as_tibble %>% 
      # filter(!is.na(padj) & padj != 1) %>%
      mutate(highFC=(abs(log2FoldChange)>LOG_FC_CUTOFF),
             lowPADJ=(padj<=PADJ_CUTOFF),
             signif=((abs(log2FoldChange)>LOG_FC_CUTOFF) & (padj<=PADJ_CUTOFF)),
             neglog10=-log10(pvalue)) 
    
    p[[plot_name]] = 
      ggplot(data2plot, aes(x=log2FoldChange, y=neglog10)) +
      geom_point(aes(colour = signif), 
                 size=1.5) +
      scale_colour_manual(values = c("TRUE"= "red", "FALSE"= "black")) + 
      xlab(expression(log[2]~fold~change)) +
      ylab(expression(-log[10]~pvalue)) +
      geom_vline(xintercept = c(-LOG_FC_CUTOFF, LOG_FC_CUTOFF),color = "blue", lty = 2) +
      ggtitle(paste0("Volcano plot ",contrast_name)) 

  
  }
  
  pdf(pdf_file_name)
  marrangeGrob(p, ncol=2, nrow=2)

  dev.off() 
  
}

plot_expression_heatmap = function (mat2plot, EFFECTS, stats_df, col_data, file_name, color_range="red2green", row_distance_measure="correlation", plot_title = NA) {

  # heatmap row annotation
  row_annotation <- stats_df[, str_detect(names(stats_df),"pass") & !names(stats_df)=="pass_any" & !names(stats_df)=="pass_combined",drop=F]
  colnames(row_annotation) <- colnames(row_annotation) %>% str_replace("pass.","")
  rownames(row_annotation) <- stats_df$gene
  row_annotation <- row_annotation[row.names(mat2plot),,drop=F]
  row_annotation[row_annotation==""] <- NA
  #remove row_annotation columns (a.k.a. comparisons) with all NA (0 DE genes), otherwise pheatmap throws error (Vered)
  row_annotation = row_annotation %>% select_if (function(x) any(!is.na(x)))

  # heatmap column annotation
  col_annotation = col_data %>% select(all_of(EFFECTS))

  col_annotation$sample <- NULL
  
      
  # Z-scoring: scale expression data by row
  mat2plot <- z_score(mat2plot)

  #Vered 24.9.2020: remove rows which contain NaN values (otherwise pheatmap throws an error)
  mat2plot <- mat2plot[complete.cases(mat2plot), ]
  
  #make color scale
  hmcol = switch (color_range,
                  "red2green"  =colorRampPalette(c("red", "black", "green"))(100),
                  "red2green.via.gray"  =colorRampPalette(c("brown2", "grey", "seagreen4"))(100),
                  "yellow2blue"=colorRampPalette(viridis(n=3) %>% rev)(255),
                  "red2blue"   =colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100),
                  "red2blue_W" =colorRampPalette(c("red", "white", "blue"))(100)
  )
  
  png(filename = file_name,
      width = 28,
      height = 20,
      units = "cm",
      res=300)
  pheatmap_data <-
    pheatmap(mat2plot,
             clustering_method = "complete",
             clustering_distance_rows = row_distance_measure,
             col = rev(hmcol),
             cluster_cols = F,
             annotation_row = row_annotation,
             annotation_col = col_annotation,
             main = plot_title,
             show_rownames = F,
             width = 15, height = 10                       # Comment these two lines for stdout
    )
  dev.off()
  
  return (pheatmap_data)
  
}

plot_expression_heatmap_wo_contrasts = function (mat2plot, EFFECTS, stats_df_interaction, col_data, file_name, color_range="red2green", row_distance_measure="correlation", plot_title = NA) {
  
  # heatmap row annotation
  row_annotation <- stats_df_interaction[, str_detect(names(stats_df_interaction),"pass") ,drop=F]
  colnames(row_annotation) <- colnames(row_annotation) %>% str_replace("pass.","")
  rownames(row_annotation) <- stats_df$gene
  row_annotation <- row_annotation[row.names(mat2plot),,drop=F]
  row_annotation[row_annotation==""] <- NA
  #remove row_annotation columns (a.k.a. comparisons) with all NA (0 DE genes), otherwise pheatmap throws error (Vered)
  row_annotation = row_annotation %>% select_if (function(x) any(!is.na(x)))
  
  # heatmap column annotation
  col_annotation = col_data %>% select(all_of(EFFECTS))
  
  col_annotation$sample <- NULL
  
  
  # Z-scoring: scale expression data by row
  mat2plot = z_score(mat2plot)
  
  #Vered 24.9.2020: remove rows which contain NaN values (otherwise pheatmap throws an error)
  mat2plot <- mat2plot[complete.cases(mat2plot), ]
  
  #make color scale
  hmcol = switch (color_range,
                  "red2green"  =colorRampPalette(c("red", "black", "green"))(100),
                  "red2green.via.gray"  =colorRampPalette(c("brown2", "grey", "seagreen4"))(100),
                  "yellow2blue"=colorRampPalette(viridis(n=3) %>% rev)(255),
                  "red2blue"   =colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100),
                  "red2blue_W" =colorRampPalette(c("red", "white", "blue"))(100)
  )
  
  png(filename = file_name,
      width = 28,
      height = 20,
      units = "cm",
      res=300)
  pheatmap_data <-
    pheatmap(mat2plot,
             clustering_method = "complete",
             clustering_distance_rows = row_distance_measure,
             col = rev(hmcol),
             cluster_cols = F,
             #annotation_row = row_annotation,
             annotation_col = col_annotation,
             main = plot_title,
             show_rownames = F,
             width = 15, height = 10                       # Comment these two lines for stdout
    )
  dev.off()
  
  return (pheatmap_data)
  
}

plot_expression_heatmap_for_precomputed_clusters = function (mat2plot, clusters, EFFECTS, GROUP, GROUP1, col_data, heatmap_file_name, row_distance_measure="correlation", plot_title = NA) {

  #code adjusted from Liron's DeSeq2 module
    
  New_clusters=c()
  for (num in sort(unique(clusters))){  #for each cluster
    #get a list of the genes in this cluster
    res_red=unique(names(clusters))
    Genes=res_red[res_red  %in% names(clusters[clusters %in% num])]
    
    if (length(Genes)>1){
      
      heatmap_title = paste0("Cluster ", num, " (", length(Genes), " genes)")
      
      #create a heatmap of the genes in the cluster
      heat_map=pheatmap(mat = mat2plot[Genes,],
                        cluster_rows = T,
                        show_rownames = F,
                        clustering_distance_rows = "correlation" ,
                        clustering_method = "ward.D2",
                        main = heatmap_title,
                        cluster_cols = F,
                        silent = F,   
                        scale = "row")
      
      
      file_per_cluster = str_replace(heatmap_file_name, ".png", paste0("_cluster", num, ".png"))
      ggsave(file_per_cluster,heat_map$gtable,dpi = 600, width = 15, height = 15, units = "cm",scale = 1.5)
      
      #add row order as computed by the hierarchical clustering
      #an array where element names are the genes and values are cluster no.
      New_clusters=c(New_clusters,clusters[heat_map$tree_row$labels[heat_map$tree_row$order]])
    }else{
      New_clusters=c(New_clusters,clusters[Genes])
    }
  }
  
  #prepare for final clustering  
  
  clusters=New_clusters
  
  #order the expression matrix by the order from the clustering (each cluster after the other,
  #and within each cluster, the order from the sample-wise hierarchical clustering)
  
  mat2plot_ordered = mat2plot[names(clusters),]
  
  #Vered 24.9.2020: remove rows which contain NaN values (otherwise pheatmap throws an error)
  mat2plot_ordered <- mat2plot_ordered[complete.cases(mat2plot_ordered), ]
  
  #create column annotation (sample names) for the heatmap
  #creates a data frame where row names are sample names and there is one column called X_AXIS with the biol groups
  
  #choose one of the options. if GROUP is not equal to GROUP1 you may need to change the order
  annotation_col = col_data[GROUP] #to show only the main effect (e.g. Diet)
  annotation_col = col_data %>% select(all_of(EFFECTS))  #to show all effects
  
  #create row annotation (cluster no.) for the heatmap
  #create a data frame where row names are genes and there is one column called Clusters with the cluster no. clusters are factors
  annotation_row = data.frame("Clusters" = factor(clusters))
  
  #define row gaps
  row_gaps = cumsum(aggregate.data.frame(x = as.data.frame(clusters),by =list(clusters),FUN = length )[,"clusters"] )
  
  #define col gaps
  if (GROUP == GROUP1) {
    col_gaps = cumsum(aggregate.data.frame(x = annotation_col[GROUP],by =annotation_col[GROUP],FUN = length)[,-1])
  } else {
    #check this! may need to switch between GROUP and GROUP1. Note, this required that the samples are originally ordered as expected
    #(otherwise, see code from Liron how to reorder them according the needs of this hiererchical clustering)
    col_gaps = c(cumsum(rev(aggregate.data.frame(x = annotation_col[GROUP1],by =annotation_col[c(GROUP1,GROUP)],FUN = length))[,1]),
                 rep(cumsum(rev(aggregate.data.frame(x = annotation_col[GROUP],by =annotation_col[c(GROUP)],FUN = length))[,1]),2))    
  }
  
  #do clustering
  
  #create heatmap
  heat_map=pheatmap(mat = mat2plot_ordered ,
                    cluster_rows= F,   #don't cluster rows. they are already ordered.
                    show_rownames=F,
                    cluster_cols=F,
                    scale = "row",
                    silent = F,
                    annotation_row=annotation_row,
                    #annotation_colors = annotation_row_COLORS,
                    #annotation_names_row=T,
                    annotation_col=annotation_col,
                    show_colnames = T,
                    gaps_row = row_gaps,
                    gaps_col = col_gaps,
                    border_color=NA
  )
  
  ggsave(heatmap_file_name,heat_map$gtable,dpi = 600, width = 15, height = 15, units = "cm",scale = 1.5)
  
  return(New_clusters)
}

plot_cluster_profiles<-function(clusters,mat2plot,col_data, color_group=c('Type','Time'),titles=c("Type","Time","Normalized counts"),split_by=c(),smooth=T,X_AXIS_ORDER=NA, out_file = ""){

  #code adjusted from Liron's DeSeq2 module
  
  #arguments hard coded (for testing)
  
  # color_group = c(GROUP1,GROUP)
  # titles = c(GROUP1,GROUP,"Normalized counts")
  # split_by=c()
  # smooth=T
  # X_AXIS_ORDER=unique(col_data[GROUP])
  # out_file = cluster_profiles_file
  
  vs_ano = cbind(col_data[c(GROUP1, GROUP)], t(mat2plot))  #should be col_data[color_group]
  
  plot_list  = list()
  title_list = list()
  count=1
  for (i in sort(unique(clusters))){
    Df=data.frame()
    genes=names(clusters[clusters==i])
    for (j in genes){
      temp=subset.data.frame(x =vs_ano,,c(color_group,j))
      colnames(temp)=c("Type","Time","EXP")
      Df=rbind(Df,temp)
    }
    if (color_group[1]==color_group[2]){
      Df$Type='Trend'
    }
    if (length(X_AXIS_ORDER)>1){
      Df$Time <-factor(Df$Time,levels=intersect(X_AXIS_ORDER,unique(Df$Time)))
    }
    if (length(split_by)>0){
      title_list[count]=list(c(i,length(genes)) ) 
      plot_list[count]=list(ggplot(data=Df, aes(x=Time, y=EXP, group=Type)) + 
                              facet_wrap(~split_by, ncol=1,scales = "free_y",strip.position = c("right"))+
                              theme(strip.text.y = element_text(size = 7, colour = "black"))+#,face="bold" ,angle = 90))+
                              #theme(strip.background = element_blank())+
                              ggtitle( paste(paste("Cluster ", i  ,sep="") ,length(genes) ,sep="\n")) +
                              #theme(plot.title = element_text(colour = "black", size = 12)) + 
                              theme(legend.position  ="none")+
                              xlab(titles[2]) +
                              ylab(titles[3]) +
                              theme(axis.title   = element_text(colour = "black", size = 10) )+
                              theme(legend.title = element_text(colour = "black", size = 10) )+
                              theme(legend.text  = element_text(colour = "black", size = 8) )+
                              theme(axis.text.y  = element_text(colour = "black", size = 6))+
                              theme(axis.text.x  = element_text(colour = "black", size = 6))+
                              theme(plot.margin  = unit(c(0.5,0.2,0.5,0.2),"cm")  )+
                              #scale_color_discrete(name=titles[1])+
                              #scale_linetype_manual(name=titles[1],values=c(1,5))+
                              scale_x_discrete(limits=unique(Df$Time))+
                              theme(legend.key.width =unit(3,"line")) 
                            #scale_y_continuous(breaks=c(seq(0,10,by=2)) )
      )
    }else{
      title_list[count]=list(c(i,length(genes)) ) 
      plot_list[count]=list(ggplot(data=Df, aes(x=Time, y=EXP, group=Type)) + 
                              ggtitle( paste(paste("Cluster ", i  ,sep="") ,length(genes) ,sep="\n")) +
                              #theme(plot.title = element_text(colour = "black", size = 12)) + 
                              theme(legend.position  ="none")+
                              xlab(titles[2]) +
                              ylab(titles[3]) +
                              theme(axis.title   = element_text(colour = "black", size = 10) )+
                              theme(legend.title = element_text(colour = "black", size = 10) )+
                              theme(legend.text  = element_text(colour = "black", size = 8) )+
                              theme(axis.text.y  = element_text(colour = "black", size = 6))+
                              theme(axis.text.x  = element_text(colour = "black", size = 6))+
                              theme(plot.margin  = unit(c(0.5,0.2,0.5,0.2),"cm")  )+
                              #scale_color_discrete(name=titles[1])+
                              #scale_linetype_manual(name=titles[1],values=c(1,5))+
                              scale_x_discrete(limits=unique(Df$Time))+
                              theme(legend.key.width = unit(3,"line")) 
                            #scale_y_continuous(breaks=c(seq(0,10,by=2)) )
      )
    }
    
    
    if (smooth){
      if ( (dim(Df)[1]<30000) && (length(unique(Df$Time))>3)){
        old_plot = plot_list[[count]] 
        res <- try( list(plot_list[[count]] + 
                           stat_smooth(method="loess",
                                       fullrange=F,
                                       size=0.5,
                                       span=1,
                                       aes(linetype=Type,color=Type)) ),
                    silent = T)
        if (inherits(res,"try-error")){
          plot_list[count]=list(old_plot +
                                  stat_smooth(method='glm',
                                              formula = y ~ poly(x, length(unique(Df$Time))-1),
                                              fullrange=F,
                                              size=0.5,
                                              aes(linetype=Type,color=Type)) 
          )
          # plot_list[count]=list(plot_list[[count]] + 
          # stat_smooth(method="loess", span = 0.1 ,fullrange=F, size=0.5 ,aes(linetype=Type,color=Type)) 
          # )
        }else{
          plot_list[count]=res
          
        }
      }else{
        plot_list[count]=list(plot_list[[count]]+
                                stat_smooth(method='glm',
                                            formula = y ~ poly(x, length(unique(Df$Time))-1),
                                            fullrange=F,
                                            size=0.5,
                                            aes(linetype=Type,color=Type)) 
        )
      }
      #below is Liron's code which gave bad results for me
      # gp2=plot_list[count]
      # res<-try(get_legend(gp2[[1]] + theme(legend.position  ="right")+guides(color=guide_legend(title=titles[1]),linetype=guide_legend(title=titles[1]))) ,silent = T)
      # if (inherits(res,"try-error")){
      #   plot_list[count]=list(plot_list[[count]]+
      #                           stat_summary(fun.y=mean, geom="point", size = 3,aes(color=Type))+
      #                           stat_summary(fun.data = "mean_se", geom = "errorbar", width = .3,size = 1,aes(color=Type),show.legend = F)
      #   )
      # }
    }else{
      plot_list[count]=list(plot_list[[count]]+
                              stat_summary(fun.y=mean, geom="line", size = 1.3,aes(linetype=Type,color=Type))+
                              stat_summary(fun.data = "mean_se", geom = "errorbar", width = .3,size = 1,aes(color=Type),show.legend = F)
      )
    }
    
    out_file1 = str_replace(out_file, ".pdf", paste0("_cluster", i, ".png"))
    ggsave(out_file1,plot_list[[count]],dpi = 600, width = 3, height =3, units = "in")
    
    count=count+1  
  }
  
  gp2=plot_list[1]
  num_of_plots=length(plot_list)
  if (num_of_plots <= 2) {
    grob_rows=1
    grob_cols=num_of_plots
  } else if (num_of_plots <= 4) {
    grob_rows = grob_cols = 2
  } else {
    grob_rows=2
    grob_cols=3
  }
  if (color_group[1]==color_group[2]){
    Clustering_Plot<-marrangeGrob(plot_list, nrow=grob_rows, ncol=grob_cols, top=NULL)
  } else {
    legend <- get_legend(gp2[[1]] + theme(legend.position  ="right")+guides(color=guide_legend(title=titles[1]),linetype=guide_legend(title=titles[1])))
    Clustering_Plot<-marrangeGrob(plot_list, nrow=grob_rows, ncol=grob_cols,right=legend, top=NULL)
  }
  
  
  #Clustering_Plot = clusters_plots[[1]]
  
  ggsave(out_file,Clustering_Plot,dpi = 600, width = 6.99, height =3.99, units = "in")
  
  #return(list(Clustering_Plot))
  return(Clustering_Plot)
}

##### Results to Excel #####

write_all_genes_to_Excel = function (res_df) {
  

  
  # index in grouping_header for which to produce column info above the table (2=counts data)
  grouping_header_meta <- c(2)
  
  # Remove leading Xs from sample names
  res_df2write = remove_first_x_from_colnames_starting_with_xdd (res_df)
  
  res_df_formula_columns <-
    res_df2write %>%
    names %>%
    str_detect(pattern = "manual_cutoffs") %>%
    which
  
  for (coli in res_df_formula_columns) {
    class(res_df2write[,coli]) <- c(class(res_df2write[,coli]), "formula")
  }
  
  # Set to FALSE if there will be more additions to the workbook in further sections
  save_workbook <- TRUE
  
  # Writing data to an excel workbook
  wb <- createWorkbook()
  
  # Add first sheet: cutoffs form
  addWorksheet(wb, 
               sheetName = "Cutoffs", 
               gridLines = TRUE)
  writeData(wb, "Cutoffs", 
            x = c("p-value","Adjusted pvalue (FDR)","linear Fold Change (linearFC)"), 
            startCol = 1, 
            startRow = 4)
  writeData(wb, "Cutoffs", 
            x = c(1,0.05,1.3), 
            startCol = 2, 
            startRow = 4)
  
  createNamedRegion(wb, "Cutoffs", cols=2, rows=4, "PVAL_CO")
  createNamedRegion(wb, "Cutoffs", cols=2, rows=5, "FDR_CO")
  createNamedRegion(wb, "Cutoffs", cols=2, rows=6, "LFC_CO")
  
  style_COs <- 
    createStyle(border = "TopBottomLeftRight",
                borderStyle = "thick", 
                fgFill = "yellow")
  addStyle(wb, "Cutoffs", style_COs, rows=4:6, cols=2, gridExpand = FALSE, stack = FALSE)
  
  
  # Add second sheet: results for all genes:
  sheetname = "All genes"
  
  addWorksheet(wb, 
               sheetName = sheetname, 
               gridLines = TRUE)
  
  
  writeDataTable(wb, 
                 sheet = sheetname, 
                 x = res_df2write,# add this if for testing: %>% head(n=2150),
                 startRow = metadata_rows+1,
                 colNames = TRUE,
                 rowNames = FALSE,
                 tableStyle = "TableStyleLight9")
  
  
  s_num_cs_round <- createStyle(numFmt = "#,##0")
  
  addStyle(wb, 
           sheet = sheetname, 
           style =  s_num_cs_round, 
           rows = (metadata_rows+2):(metadata_rows+1+dim(res_df2write)[1]), 
           cols = (grouping_header[1]+1):sum(grouping_header[1:2]), 
           gridExpand = TRUE, stack = FALSE)
  
  # Adding additional column information above column headers
  writeData(wb, 
            sheet = sheetname, 
            x = contrasts_header %>% as.matrix %>% t,
            startRow = metadata_rows,
            colNames = F,
            rowNames = F)
  
  # Adding sample meta-data (col_data) where required:
  for (i in cumsum(grouping_header)[grouping_header_meta-1]){
    writeData(wb, 
              sheet = sheetname, 
              x = col_data[,metadata_cols] %>% t , 
              startRow = 1,
              startCol = i,
              colNames = FALSE,
              rowNames = TRUE)
  }
  
  
  #
  # 1. Create style
  # 2. Add style to location
  # Creating `colorscheme` with background colors and appropriate text color, based on luminence (brightness)
  
  mycolors = rep(brewer.pal(12, name = "Paired"),3)
  colorscheme <- data.frame(colnum = grouping_header,
                            fgcol= mycolors[1:length(grouping_header)],
                            stringsAsFactors = F)
  
  colorscheme$textcol <- ifelse((as(hex2RGB(colorscheme$fgcol),"polarLUV"))@coords[,1] > 65, 
                                yes = "dimgray",
                                no = "white")
  
  
  
  for (i in 1:length(grouping_header)) {
    cat(sprintf("%s of %s\n",i,length(grouping_header)))
    # Add title color:
    hs1 <- createStyle(fgFill = colorscheme$fgcol[i], 
                       fontColour = colorscheme$textcol[i])
    addStyle(wb, 
             sheet = sheetname, 
             style =  hs1, 
             rows = metadata_rows+1, 
             cols = ifelse(i==1,1,cumsum(grouping_header)[i-1]+1):cumsum(grouping_header)[i], 
             gridExpand = FALSE, stack = TRUE)
    # Add group separator:
    hs2 <- createStyle(border = "Right",
                       borderColour = "black",
                       borderStyle = "thick")
    addStyle(wb = wb, 
             sheet = sheetname, 
             style = hs2, 
             rows = (metadata_rows+2):(metadata_rows+1+dim(res_df2write)[1]), 
             cols = cumsum(grouping_header)[i], 
             gridExpand = FALSE, stack = TRUE)
    
  }
  
  if (save_workbook) {
    
      # Writing workbook:
      print(sprintf("Writing file: %s", final_results_ALL_file))
      # openXL(wb)
      saveWorkbook(wb, 
                   file = final_results_ALL_file, 
                   overwrite = TRUE) ## save to working directory
    
  }
  
  return (res_df2write)
  
}

write_DE_genes_to_Excel = function (res_df2write) {
  
  # Get only differentially expressed:
  res_df2write_DE <- res_df2write[(!is.na(res_df2write$pass_any) & res_df2write$pass_any=="1"),]
  
  # Order of genes in heatmap.

  #note, I think that left join is used when adding order to res_df2write_DE and not cbind 
  #because top_DE_genes may include only part of the DE genes
  
  top_diff_genes_order <- 
    cbind(gene=top_DE_genes[pheatmap_data_DE_genes$tree_row$order],
          order = seq(from=1,by=1,along.with = top_DE_genes)) %>% 
    as_tibble()
  
  res_df2write_DE <-
    res_df2write_DE %>% 
    as_tibble %>% 
    left_join(y = top_diff_genes_order,by = c("gene"))
  
  #add cluster number and order of DE genes according to clustering using DeSeq2 module method

  top_diff_genes_order1 = data.frame(gene=names(New_clusters),
                                     cluster=New_clusters,
                                     order_clusters = seq(from=1,by=1,along.with = New_clusters)) %>%
    as_tibble()
  
  res_df2write_DE <-
    res_df2write_DE %>% 
    as_tibble %>% 
    left_join(y = top_diff_genes_order1,by = c("gene"))
 
  #add cluster number and order of DE genes according to manual clustering
  
  if (PERFORM_MANUAL_CLUSTERING) {
    
    top_diff_genes_order2 = data.frame(gene=names(New_clusters_man),
                                       man_cluster=New_clusters_man,
                                       order_man_clusters = seq(from=1,by=1,along.with = New_clusters_man)) %>%
      as_tibble()
    
    res_df2write_DE <-
      res_df2write_DE %>% 
      as_tibble %>% 
      left_join(y = top_diff_genes_order2,by = c("gene"))    
  }
   
  # Add Z-score data
  
  zscore = z_score(norm_log_counts)
  colnames(zscore) <- zscore %>% colnames %>% paste0(.,".zscore")
  
  zscore <- zscore %>% as_tibble(rownames = "gene") 
  
  res_df2write_DE <- 
    res_df2write_DE %>% 
    left_join(zscore,by="gene") 
  
  dim(res_df2write_DE)
  res_df2write_DE <- res_df2write_DE %>% unique()
  dim(res_df2write_DE)
  
  ### Output to excel
 
  # Set to FALSE if there will be more additions to the workbook in further sections
  save_workbook <- TRUE
  
  sheetname = "DE genes"
  
  grouping_header_meta_DE <- c(2,length(grouping_header_DE))
  
  
  # Remove leading Xs from sample names
  names(res_df2write_DE) <- 
    res_df2write_DE %>% 
    names %>% 
    str_replace(pattern = "^X(\\d)",
                replacement = "\\1")
  
  # Writing data to an excel workbook:
  
  if (results_all_with_DE) {
    wb1 = wb
  } else {
    wb1 = createWorkbook()
  }
  
  # Add worksheet and data:
  addWorksheet(wb1, 
               sheetName = sheetname, 
               gridLines = TRUE)
  
  
  writeDataTable(wb1, 
                 sheet = sheetname, 
                 x = res_df2write_DE ,     #%>% head(n=2150),
                 startRow = metadata_rows+1,
                 colNames = TRUE,
                 rowNames = FALSE,
                 tableStyle = "TableStyleLight9")
  
  
  s_num_cs_round <- createStyle(numFmt = "#,##0")
  
  addStyle(wb1, 
           sheet = sheetname, 
           style =  s_num_cs_round, 
           rows = (metadata_rows+2):(metadata_rows+1+dim(res_df2write_DE)[1]), 
           cols = (grouping_header_DE[1]+1):sum(grouping_header_DE[1:2]), 
           gridExpand = TRUE, stack = FALSE)
  
  writeData(wb1, 
            sheet = sheetname, 
            x = contrasts_header_DE %>% as.matrix %>% t,
            startRow = metadata_rows,
            colNames = F,
            rowNames = F)
  
  
  # Add col_data where required
  for (i in cumsum(grouping_header_DE)[grouping_header_meta_DE-1]){
    writeData(wb1, 
              sheet = sheetname, 
              x = col_data[,metadata_cols] %>% t , 
              startRow = 1,
              startCol = i,
              colNames = FALSE,
              rowNames = TRUE)
  }
  
  
  #
  # 1. Create style
  # 2. Add style to location
  # Creating `colorscheme` with background colors and appropriate text color, based on luminence (brightness)
  
  mycolors = rep(brewer.pal(12, name = "Paired"),3)
  colorscheme <- data.frame(colnum = grouping_header_DE,
                            fgcol= mycolors[1:length(grouping_header_DE)],
                            stringsAsFactors = F)
  
  colorscheme$textcol <- ifelse((as(hex2RGB(colorscheme$fgcol),"polarLUV"))@coords[,1] > 65, 
                                yes = "dimgray",
                                no = "white")

  
  for (i in 1:length(grouping_header_DE)) {
    cat(sprintf("%s of %s\n",i,length(grouping_header_DE)))
    # Add title color:
    hs1 <- createStyle(fgFill = colorscheme$fgcol[i], 
                       fontColour = colorscheme$textcol[i])
    addStyle(wb1, 
             sheet = sheetname, 
             style =  hs1, 
             rows = metadata_rows+1, 
             cols = ifelse(i==1,1,cumsum(grouping_header_DE)[i-1]+1):cumsum(grouping_header_DE)[i], 
             gridExpand = FALSE, stack = TRUE)
    # Add group separator:
    hs2 <- createStyle(border = "Right",
                       borderColour = "black",
                       borderStyle = "thick")
    addStyle(wb = wb1, 
             sheet = sheetname, 
             style = hs2, 
             rows = (metadata_rows+2):(metadata_rows+1+dim(res_df2write_DE)[1]), 
             cols = cumsum(grouping_header_DE)[i], 
             gridExpand = FALSE, stack = TRUE)
    
  }
  
  if (save_workbook) {

      print(sprintf("Writing file: %s", final_results_DE_file))
      saveWorkbook(wb1, 
                   file = final_results_DE_file, 
                   overwrite = TRUE) ## save to working directory

  }
  

}

write_interaction_genes_to_Excel = function (res_df2write) {
  
  # Get only genes with significant interaction
  res_df2write_DE <- res_df2write[(!is.na(res_df2write$pass.interaction) & res_df2write$pass.interaction=="interaction"),]
  
  
  # Order of genes in heatmap
  
  top_diff_genes_order <- 
    cbind(gene=top_DE_genes_interaction[pheatmap_data_interaction_genes$tree_row$order],
          order = seq(from=1,by=1,along.with = top_DE_genes_interaction)) %>% 
    as_tibble()
  
  res_df2write_DE <-
    res_df2write_DE %>% 
    as_tibble %>% 
    left_join(y = top_diff_genes_order,by = c("gene"))
  
  # Adding Z-score data
  zscore = z_score(norm_log_counts)
  colnames(zscore) <- zscore %>% colnames %>% paste0(.,".zscore")
  
  zscore <- zscore %>% as_tibble(rownames = "gene") 
  
  res_df2write_DE <- 
    res_df2write_DE %>% 
    left_join(zscore,by="gene") 
  
  dim(res_df2write_DE)
  res_df2write_DE <- res_df2write_DE %>% unique()
  dim(res_df2write_DE)
  
  ### Output to excel
  
  # Set to FALSE if there will be more additions to the workbook in further sections
  save_workbook <- TRUE
  
  sheetname = "Interaction genes"
  
  grouping_header_meta_DE <- c(2,length(grouping_header_DE))
  
  # Remove leading Xs from sample names
  names(res_df2write_DE) <- 
    res_df2write_DE %>% 
    names %>% 
    str_replace(pattern = "^X(\\d)",
                replacement = "\\1")
  
  # Writing data to an excel workbook
  if (results_all_with_DE) {
    wb1 = wb
  } else {
    wb1 = createWorkbook()
  }
  
  # Add worksheet and data
  addWorksheet(wb1, 
               sheetName = sheetname, 
               gridLines = TRUE)
  
  writeDataTable(wb1, 
                 sheet = sheetname, 
                 x = res_df2write_DE ,
                 startRow = metadata_rows+1,
                 colNames = TRUE,
                 rowNames = FALSE,
                 tableStyle = "TableStyleLight9")
  
  
  s_num_cs_round <- createStyle(numFmt = "#,##0")
  
  addStyle(wb1, 
           sheet = sheetname, 
           style =  s_num_cs_round, 
           rows = (metadata_rows+2):(metadata_rows+1+dim(res_df2write_DE)[1]), 
           cols = (grouping_header_DE[1]+1):sum(grouping_header_DE[1:2]), 
           gridExpand = TRUE, stack = FALSE)
  
  writeData(wb1, 
            sheet = sheetname, 
            x = contrasts_header_DE %>% as.matrix %>% t,
            startRow = metadata_rows,
            colNames = F,
            rowNames = F)
  
  
  # Adding col_data where required:
  for (i in cumsum(grouping_header_DE)[grouping_header_meta_DE-1]){
    writeData(wb1, 
              sheet = sheetname, 
              x = col_data[,metadata_cols] %>% t , 
              startRow = 1,
              startCol = i,
              colNames = FALSE,
              rowNames = TRUE)
  }
  
  # 1. Create style
  # 2. Add style to location
  # Creating `colorscheme` with background colors and appropriate text color, based on luminence (brightness)
  
  mycolors = rep(brewer.pal(12, name = "Paired"),3)
  colorscheme <- data.frame(colnum = grouping_header_DE,
                            fgcol= mycolors[1:length(grouping_header_DE)],
                            stringsAsFactors = F)
  
  colorscheme$textcol <- ifelse((as(hex2RGB(colorscheme$fgcol),"polarLUV"))@coords[,1] > 65, 
                                yes = "dimgray",
                                no = "white")
  
  for (i in 1:length(grouping_header_DE)) {
    cat(sprintf("%s of %s\n",i,length(grouping_header_DE)))
    # Add title color:
    hs1 <- createStyle(fgFill = colorscheme$fgcol[i], 
                       fontColour = colorscheme$textcol[i])
    addStyle(wb1, 
             sheet = sheetname, 
             style =  hs1, 
             rows = metadata_rows+1, 
             cols = ifelse(i==1,1,cumsum(grouping_header_DE)[i-1]+1):cumsum(grouping_header_DE)[i], 
             gridExpand = FALSE, stack = TRUE)
    # Add group separator:
    hs2 <- createStyle(border = "Right",
                       borderColour = "black",
                       borderStyle = "thick")
    addStyle(wb = wb1, 
             sheet = sheetname, 
             style = hs2, 
             rows = (metadata_rows+2):(metadata_rows+1+dim(res_df2write_DE)[1]), 
             cols = cumsum(grouping_header_DE)[i], 
             gridExpand = FALSE, stack = TRUE)
    
  }
  
  if (save_workbook) {
    

      print(sprintf("Writing file: %s", final_results_Interaction_file))
      saveWorkbook(wb1, 
                   file = final_results_Interaction_file, 
                   overwrite = TRUE) ## save to working directory
  }
  
}

write_combined_genes_to_Excel = function (res_df2write) {
  
  # Get only differentially expressed:
  res_df2write_DE <- res_df2write[(!is.na(res_df2write$pass_combined) & res_df2write$pass_combined=="1"),]
  
  # Order of genes in heatmap
  
  top_diff_genes_order <- 
    cbind(gene=top_DE_genes_combined[pheatmap_data_combined_genes$tree_row$order],
          order = seq(from=1,by=1,along.with = top_DE_genes_combined)) %>% 
    as_tibble()
  
  res_df2write_DE <-
    res_df2write_DE %>% 
    as_tibble %>% 
    left_join(y = top_diff_genes_order,by = c("gene"))
  
  # Adding Z-score data
  zscore = z_score(norm_log_counts)
  colnames(zscore) <- zscore %>% colnames %>% paste0(.,".zscore")
  
  zscore <- zscore %>% as_tibble(rownames = "gene") 
  
  res_df2write_DE <- 
    res_df2write_DE %>% 
    left_join(zscore,by="gene") 

  
  dim(res_df2write_DE)
  res_df2write_DE <- res_df2write_DE %>% unique()
  dim(res_df2write_DE)
  
  ### Output to excel
  
  # Set to FALSE if there will be more additions to the workbook in further sections
  save_workbook <- TRUE
  
  sheetname = "Combined"
  
  grouping_header_meta_DE <- c(2,length(grouping_header_DE))
  
  
  # Remove leading Xs from sample names
  names(res_df2write_DE) <- 
    res_df2write_DE %>% 
    names %>% 
    str_replace(pattern = "^X(\\d)",
                replacement = "\\1")
  
  if (results_all_with_DE) {
    wb1 = wb
  } else {
    wb1 = createWorkbook()
  }
  
  # Add worksheet and data:
  addWorksheet(wb1, 
               sheetName = sheetname, 
               gridLines = TRUE)
  
  writeDataTable(wb1, 
                 sheet = sheetname, 
                 x = res_df2write_DE,
                 startRow = metadata_rows+1,
                 colNames = TRUE,
                 rowNames = FALSE,
                 tableStyle = "TableStyleLight9")
  
  
  s_num_cs_round <- createStyle(numFmt = "#,##0")
  
  addStyle(wb1, 
           sheet = sheetname, 
           style =  s_num_cs_round, 
           rows = (metadata_rows+2):(metadata_rows+1+dim(res_df2write_DE)[1]), 
           cols = (grouping_header_DE[1]+1):sum(grouping_header_DE[1:2]), 
           gridExpand = TRUE, stack = FALSE)
  
  writeData(wb1, 
            sheet = sheetname, 
            x = contrasts_header_DE %>% as.matrix %>% t,
            startRow = metadata_rows,
            colNames = F,
            rowNames = F)
  
  # Add col_data where required
  for (i in cumsum(grouping_header_DE)[grouping_header_meta_DE-1]){
    writeData(wb1, 
              sheet = sheetname, 
              x = col_data[,metadata_cols] %>% t , 
              startRow = 1,
              startCol = i,
              colNames = FALSE,
              rowNames = TRUE)
  }
  

  # 1. Create style
  # 2. Add style to location
  # Creating `colorscheme` with background colors and appropriate text color, based on luminence (brightness)
  
  mycolors = rep(brewer.pal(12, name = "Paired"),3)
  colorscheme <- data.frame(colnum = grouping_header_DE,
                            fgcol= mycolors[1:length(grouping_header_DE)],
                            stringsAsFactors = F)
  
  colorscheme$textcol <- ifelse((as(hex2RGB(colorscheme$fgcol),"polarLUV"))@coords[,1] > 65, 
                                yes = "dimgray",
                                no = "white")

  
  for (i in 1:length(grouping_header_DE)) {
    cat(sprintf("%s of %s\n",i,length(grouping_header_DE)))
    # Add title color
    hs1 <- createStyle(fgFill = colorscheme$fgcol[i], 
                       fontColour = colorscheme$textcol[i])
    addStyle(wb1, 
             sheet = sheetname, 
             style =  hs1, 
             rows = metadata_rows+1, 
             cols = ifelse(i==1,1,cumsum(grouping_header_DE)[i-1]+1):cumsum(grouping_header_DE)[i], 
             gridExpand = FALSE, stack = TRUE)
    # Add group separator
    hs2 <- createStyle(border = "Right",
                       borderColour = "black",
                       borderStyle = "thick")
    addStyle(wb = wb1, 
             sheet = sheetname, 
             style = hs2, 
             rows = (metadata_rows+2):(metadata_rows+1+dim(res_df2write_DE)[1]), 
             cols = cumsum(grouping_header_DE)[i], 
             gridExpand = FALSE, stack = TRUE)
    
  }
  
  if (save_workbook) {
      print(sprintf("Writing file: %s", final_results_Combined_file))
      saveWorkbook(wb1, 
                   file = final_results_Combined_file, 
                   overwrite = TRUE) ## save to working directory
  }
}


##### Enrichment analysis #####

Filter_Pathways_By_Taxa <- function(taxon){
  ko_paths_names_all = KEGGREST::keggList("pathway", "ko")
  kegg_organisms     = KEGGREST::keggList("organism")
  kegg_organisms     = as.data.frame(kegg_organisms)
  organisms          = kegg_organisms[str_detect(kegg_organisms$phylogeny,taxon),'organism']
  
  paths = c()  #initialize paths vector
  for (org in organisms) {
    a = unique(KEGGREST::keggLink("pathway", org))  #get all pathways per species, make the list nonredundant
    print (paste (org, length(a), "pathways"))
    a = str_replace(a, "path:[a-z]+", "path:ko")  #replace species-specific with ko path name
    paths = c(paths, a)  #add to paths vector
  }
  paths_unique = unique(paths)  #remove redundancy
  print (paste ("No of unique KO pathways retrieved:", length(paths_unique)))
  paths_unique = str_replace(paths_unique, "ko", "map")
  return(paths_unique)
}

Clusters_Enrichment_Test=function(outDir,clusters,TERM2NAME,TERM2GENE,file_name,Type,pAdjustMethod='fdr',pvalueCutoff=0.05,gene2ko=FALSE,maxCategory=1000){
  
  #start
  
  allRes=list()
  original_allRes = list()
  cluster_names=list()
  count=1
  for (i in sort(unique(clusters))){
    temp=General_Enrichment_Test(clusters,c(i),TERM2NAME,TERM2GENE,pAdjustMethod,pvalueCutoff )
    if (length(temp)>0){
        allRes[[count]]<-temp
      cluster_names[count]=i
      count=count+1
    }
  } 
  names(allRes)=cluster_names
  allRes=clusterProfiler::merge_result(enrichResultList =allRes)
  if (Type=="KO"){
    allRes<-generat_urls(allRes,gene2ko)
  }
  
  enrichment_table = process_clusterprofiler_results_table (allRes@compareClusterResult)
  
  write.csv(x = enrichment_table,
            file =file.path(outDir,paste(file_name,'.csv',collapse = "")) ,
            quote = TRUE,
            row.names = TRUE)
  
  original_allRes = allRes
  Simplify_title = c()
  if (Type=="GO"){
    allRes@fun<-"enrichGO"
    clusterProfiler::simplify(allRes, cutoff=0.7, by="p.adjust", select_fun=min)
    file_name=paste("Simplify",file_name,collapse = "_")
    
    enrichment_table_simplify = process_clusterprofiler_results_table (allRes@compareClusterResult)
    
    write.csv(x = enrichment_table_simplify,
              file =file.path(outDir,paste(file_name,'.csv',sep='',collapse = "")) ,
              quote = FALSE,
              row.names = TRUE)
    
  }else{
    allRes@fun<-"enrichKEGG"
  }
  if (dim(allRes@compareClusterResult)[1]>50){
    font.size=4
  }else{
    font.size=9
  }
  
  #prepare enrichment table(s) object for delivery to RShiny
  enrich_tables = list()
  if(exists("enrichment_table")) {
    enrich_tables$all = enrichment_table
  }
  if (exists("enrichment_table_simplify")) {
    enrich_tables$simplify = enrichment_table_simplify
  }
  #draw a dot plot and return the relevant objects and graph
  
   if (dim(allRes@compareClusterResult)[1]>0){
     x=enrichplot::dotplot(allRes,showCategory=maxCategory,font.size=font.size)
     ggplot2::ggsave(filename = file.path(outDir,paste(file_name,".pdf",sep='',collapse = "")),dpi = 600,device = "pdf",width = 20,height = 20)
     #Vered 26.7.2023 - added enrichment table and enrichment_table_simplify for delivery to RShiny
     return(list(original_allRes,list(htmltools::h4(paste(unlist(stringi::stri_split(str = file_name,fixed='_')),collapse = " ")),plotly::ggplotly(x),Simplify_title), enrich_tables))
   }
  #Vered 26.7.2023 - added enrichment table and enrichment_table_simplify for delivery to RShiny
  return(list(original_allRes,list(htmltools::h4(paste(unlist(stringi::stri_split(str = file_name,fixed='_')),collapse = " "))), enrich_tables))															
}

Contrasts_Enrichment_Test=function(outDir,DE_genes_lists,TERM2NAME,TERM2GENE,file_name,Type,pAdjustMethod='fdr',pvalueCutoff=0.05,gene2ko=FALSE,maxCategory=1000){

  #start
  
  allRes=list()
  original_allRes = list()
  cluster_names=list()
  count=1
  for (key in names(DE_genes_lists)) {
    gene_list = DE_genes_lists[[key]]
    if (length(gene_list)>0) {
      temp=General_Enrichment_Test1(gene_list,TERM2NAME,TERM2GENE,pAdjustMethod,pvalueCutoff )
      if (length(temp)>0){
        allRes[[key]]<-temp
      }
    }

  }
  
  allRes=clusterProfiler::merge_result(enrichResultList =allRes)
  if (Type=="KO"){
    allRes<-generat_urls(allRes,gene2ko)
  }
  
  enrichment_table = process_clusterprofiler_results_table (allRes@compareClusterResult)
  
  write.csv(x = enrichment_table,
            file =file.path(outDir,paste(file_name,'.csv',collapse = "")) ,
            quote = TRUE,
            row.names = TRUE)
  
  original_allRes = allRes
  Simplify_title = c()
  if (Type=="GO"){
    allRes@fun<-"enrichGO"
    clusterProfiler::simplify(allRes, cutoff=0.7, by="p.adjust", select_fun=min)
    file_name=paste("Simplify",file_name,collapse = "_")
    
    enrichment_table = process_clusterprofiler_results_table (allRes@compareClusterResult)
    
    write.csv(x = enrichment_table,
              file =file.path(outDir,paste(file_name,'.csv',sep='',collapse = "")) ,
              quote = FALSE,
              row.names = TRUE)
    
  }else{
    allRes@fun<-"enrichKEGG"
  }
  if (dim(allRes@compareClusterResult)[1]>50){
    font.size=4
  }else{
    font.size=9
  }
  
  #draw a dot plot and return the relevant objects and graph
  
  if (dim(allRes@compareClusterResult)[1]>0){
    x=enrichplot::dotplot(allRes,showCategory=maxCategory,font.size=font.size)
    ggplot2::ggsave(filename = file.path(outDir,paste(file_name,".pdf",sep='',collapse = "")),dpi = 600,device = "pdf",width = 20,height = 20)
    return(list(original_allRes,list(htmltools::h4(paste(unlist(stringi::stri_split(str = file_name,fixed='_')),collapse = " ")),plotly::ggplotly(x),Simplify_title)))
  }
  return(list(original_allRes,list(htmltools::h4(paste(unlist(stringi::stri_split(str = file_name,fixed='_')),collapse = " ")))))
}

process_clusterprofiler_results_table = function(clusterprofiler_results_table) {
  
  MAX_NR_GENES_TO_SHOW = 1000
  text_to_show = paste0("Too many to show (>", MAX_NR_GENES_TO_SHOW, ")")
  
  enrichment_table = clusterprofiler_results_table %>%
    separate(col = GeneRatio, into = c("in_cluster_in_term", "in_cluster")) %>%
    separate(col = BgRatio,   into = c("in_term", "in_genome"))  %>%
    mutate_at (c('in_cluster_in_term', 'in_cluster', 'in_term', 'in_genome'), as.numeric) %>%
    mutate (Fold_enrichment = signif(((in_cluster_in_term/in_cluster)/(in_term/in_genome)),digits=2), .before = Count) %>%
    mutate (geneID = ifelse (test = Count <= MAX_NR_GENES_TO_SHOW,
                             yes  = geneID,
                             no   = text_to_show))
  return(enrichment_table)
}

General_Enrichment_Test<-function(clusters,num,TERM2NAME,TERM2GENE,pAdjustMethod='fdr',pvalueCutoff=0.05){
  
  res_red=unique(TERM2GENE[,2])
  Genes=res_red[res_red  %in% names(clusters[clusters %in% num])]
  
  res=enricher(Genes, TERM2GENE=TERM2GENE, 
               TERM2NAME=TERM2NAME,
               minGSSize     = 0,
               #maxGSSize     = length(res_red),
               maxGSSize     = 10000,
               pAdjustMethod = pAdjustMethod,
               pvalueCutoff  = pvalueCutoff,
               qvalueCutoff  = 1  #Vered added 10.8.2023
  ) 
  
  return(res)
}

General_Enrichment_Test1<-function(gene_list,TERM2NAME,TERM2GENE,pAdjustMethod='fdr',pvalueCutoff=0.05){  
  
  res_red=unique(TERM2GENE[,2])
  Genes = intersect(gene_list, res_red)
  
  res=enricher(Genes, TERM2GENE=TERM2GENE, 
               TERM2NAME=TERM2NAME,
               minGSSize     = 0,
               #maxGSSize     = length(res_red),
               maxGSSize     = 10000,
               pAdjustMethod = pAdjustMethod,
               pvalueCutoff  = pvalueCutoff,
               qvalueCutoff  = 1  #Vered added 10.8.2023
  ) 
  
  return(res)
}

plot_shared_genes<-function(allRes,outDir,file_name){
  
  heatmaps=list()
  if (length(allRes)>1){
    for (cluster in unique(allRes$Cluster)){
      TEMP = na.omit(allRes[allRes$Cluster==cluster,])   
      genes = unique(unlist(stringi::stri_split(paste(TEMP$geneID,sep = ,collapse = '/'),fixed = '/')))
      
      if ((length(TEMP$ID)>1) &&((length(genes)>1))){
        mat  = matrix(0, nrow = length(TEMP$ID), ncol = length(genes))
        rownames(mat) = TEMP$ID
        colnames(mat) = genes
        
        
        for (x in TEMP$ID){
          for (y in genes){
            mat[x,y] = length( unlist(intersect( unlist(stringi::stri_split(TEMP[TEMP$ID==x,'geneID'],fixed = '/')),y )))
            
          }
        }
        
        rownames(mat) = TEMP$Description
        
        if ((length(TEMP$ID)>20) || ((length(genes)>200))){
          heat_map = pheatmap(mat = mat,cluster_rows = T,
                              cluster_cols = T,
                              silent = T,
                              legend = F)
          
        }else{
          if (length(unique(as.vector(mat)))==2 ){
            color=colorRampPalette(c('white','pink'))(2)
            mybreaks=NA
          }else{
            color=colorRampPalette(c('pink'))(1)
            mybreaks=c(0,1)
          }
          
          heatmap_file = file.path(outDir,paste0('Cluster_',cluster,'_genes2term_',file_name,".pdf"))
          
          heat_map = pheatmap(mat = mat,cluster_rows = T,
                              treeheight_col = 0,
                              treeheight_row = 0,
                              width = 5,
                              height =5,
                              cellheight = 10,
                              fontsize_row =5,
                              fontsize_col = 1,
                              main = paste('Cluster ',cluster),
                              border_color = 'white',
                              cluster_cols = T,
                              breaks = mybreaks,
                              filename = heatmap_file,
                              color = color,
                              legend = F)
          
        }
        
        csv_file = file.path(outDir,paste0('Cluster_',cluster,'_genes2term_',file_name,".csv"))
        
        write.csv(x = mat[heat_map$tree_row$order,heat_map$tree_col$order],
                  quote = T,
                  row.names = T,
                  file = csv_file)
        
        
        mat  = matrix(0, nrow = length(TEMP$ID), ncol = length(TEMP$ID))
        rownames(mat) = TEMP$ID
        colnames(mat) = TEMP$ID
        
        for (x in TEMP$ID){
          for (y in TEMP$ID){
            mat[x,y] =100*(2*length( unlist(intersect( unlist(stringi::stri_split(TEMP[TEMP$ID==x,'geneID'],fixed = '/')),unlist(stringi::stri_split(TEMP[TEMP$ID==y,'geneID'],fixed = '/'))))))/(length(unlist(stringi::stri_split(TEMP[TEMP$ID==x,'geneID'],fixed = '/'))) +  length(unlist(stringi::stri_split(TEMP[TEMP$ID==y,'geneID'],fixed = '/'))))
            
          }
        }
        
        rownames(mat) = TEMP$Description   
        colnames(mat) = TEMP$Description   
        if (length(unique(as.vector(mat)))>1 ){
          color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu") ))( 100)
          mybreaks=seq(0,1,0.01)
        }else{
          color=colorRampPalette(c('red'))(1)
          mybreaks=c(0,1)
        }
        
        heatmap_file1 = file.path(outDir,paste0('Cluster_',cluster,'_term2term_',file_name,".pdf"))
        
        if (length(TEMP$ID)>20){
          heat_map = pheatmap(mat = mat,cluster_rows = T,
                              main = paste('Cluster ',cluster),
                              treeheight_col = 0,
                              treeheight_row = 0,
                              cluster_cols = T,
                              silent = T,
                              legend = T)
          
        }else{
          heat_map = pheatmap(mat = mat,cluster_rows = T,
                              width = 5,
                              height =5,
                              treeheight_col = 0,
                              treeheight_row = 0,
                              breaks = mybreaks,
                              fontsize_row =5,
                              fontsize_col = 5,
                              main = paste('Cluster ',cluster),
                              border_color = 'white',
                              cluster_cols = T,
                              color = color,
                              filename = heatmap_file1,
                              legend = T)
        }
        
        csv_file1 = file.path(outDir,paste0('Cluster_',cluster,'_term2term_',file_name,".csv"))
        
        write.csv(x = mat[heat_map$tree_row$order,heat_map$tree_col$order],
                  quote = T,
                  row.names = T,
                  file = csv_file1)
        
        heatmap_title=htmltools::h4(paste(unlist(stringi::stri_split(str = paste('Cluster',cluster,file_name,' Terms Genes Overlap',collapse = ""),fixed='_')),collapse = " "))
        heatmaps = c(heatmaps ,list(heatmap_title))
        mat = mat[heat_map$tree_row$order,heat_map$tree_col$order]
        heatmaps = c(heatmaps ,list(plot_ly(z = mat,
                                            x=colnames(mat),
                                            y=rownames(mat),
                                            colors = color,
                                            type = "heatmap")
        )
        )
      }
    }
    
  }
  return(heatmaps)
}

generat_urls<-function(allRes,gene2ko){
  
  #not sure this function is necessary (but it is used in other functions above)
  
  temp_table=allRes@compareClusterResult
  temp_table$URL=apply(X = temp_table,MARGIN = 1,FUN = function(x) paste(c("http://www.kegg.jp/kegg-bin/show_pathway?",stringi::stri_replace_all(str = x["ID"],replacement = "",regex = "path:",collapse = ""),"/", paste(sapply( unlist(stringi::stri_split(str = x["geneID"],regex = "/")),FUN = function(x) gene2ko[gene2ko$V1==x,"V2"]),collapse = "+") ),collapse = ""))
  allRes@compareClusterResult=temp_table
  return(allRes)
}

##### Utilities #####

remove_first_n_chars_from_col_names = function (x, n=1) {
  #by default, removes the first character from the name of each column
  i = n+1
  colnames(x) = str_sub(colnames(x), i, -1)
  return (x)
}

remove_first_x_from_all_colnames_if_exist = function (x) {
  #if all col names have leading X, remove it
  if (all(str_detect(colnames(x), "^X"))) {
    colnames(x) = str_sub(colnames(x), 2, -1)
  }
  return (x)
}

remove_first_x_from_colnames_starting_with_xdd = function (x) {
  #if colname starts with "x" follwed by two digits (i.e. it is a sample name), remove the X
  fixed_colnames = str_replace(colnames(x), pattern = "^X(\\d)", replacement = "\\1")
  colnames(x) = fixed_colnames
  return (x)
}

make_data_frame_from_clusters = function (clusters) {  #vered Aug2023
  clusters_df = data.frame(gene=names(clusters), cluster=clusters)
  return (clusters_df)
}

z_score = function (expr_matrix) {
  expr_matrix %>% t %>% scale %>% t
}

make_list_from_clusters = function (clusters) {
  clusters_list = list()
  for (cluster_name in sort(unique(clusters))){  #for each cluster get a list of the genes in this cluster  #vered Aug2023
    res_red=unique(names(clusters))
    Genes=res_red[res_red  %in% names(clusters[clusters %in% cluster_name])]
    clusters_list[[cluster_name]] = Genes
  }
  return(clusters_list)
}

export_table = function (data, file_name) {
  data <- cbind(Name = rownames(data), data)
  write.table(x = data,
              file = file_name,
              quote = F,
              sep = "\t",
              na = "",
              row.names = F,
              col.names = T)
}