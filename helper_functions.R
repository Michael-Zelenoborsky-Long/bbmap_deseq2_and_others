# Michael Gordon, 10 Oct 2020
# Functions for sgefix_combined4deduplicated_public.Rmd

# dependencies are the same as in Rmd
require(tidyverse)
require("ggbeeswarm")
require("magrittr")
require(tidyr)
require(dplyr)
library("DESeq2")


# functions require paths to data

if (exists("path_to_annot") & exists("mercator_categories_colours") &
    exists("PsT_fix_metadata_path") & exists("path_to_expr_matrix")){
  print("Paths to data are defined.")
} else {
  stop("Check paths to data.")
}


# func to attach annotation
load(path_to_annot)

addAnnotation <- function(deseq_result_object, 
                          removeDuplicated=FALSE, 
                          annot=annot_for_genome_based, 
                          al=alpha,
                          column_with_id="transcriptID") {
  #current_DESeq_result <- as.data.frame( subset(deseq_result_object, padj < al ))
  current_DESeq_result <- as.data.frame( deseq_result_object)
  current_DESeq_result$transcriptID <- row.names(current_DESeq_result)
  mergedDF <- dplyr::left_join(x = annot, y = current_DESeq_result, by= c("transcriptID"=column_with_id),  
                               copy=TRUE)
  mergedDF %<>%  dplyr::select("transcriptID", "log2FoldChange", "padj", dplyr::everything())
  if (!removeDuplicated) {
    return(subset(mergedDF, !is.na(padj) & padj < al ))
  } else {
    return( subset(mergedDF, !is.na(padj) & padj < al & !duplicated(mergedDF$transcriptID) ) )
  }
}


load(mercator_categories_colours) # loads mercatorCatsColoursDf with russian translation

mercatorTopCats <- function(annotated_df, report_empty_categ=FALSE, categories_df=mercatorCatsColoursDf, 
                            ret_raw_freq=F, ret_tmp_df=F){
  # make frequency table for an annotated dataframe
  # returns dataframe
  require(dplyr)
  
  tmp_df <- annotated_df %>% separate(Mercator_Name, into = c("mercMainCat", NA),
                                      sep='[.]', remove = F, extra = "drop", fill = "right") %>% 
    separate(Mercator_Bincode, into = c("merc1digit", NA),
             sep='[.]', remove = F, extra = "drop", fill = "right") 
  tmp_df$merc1digit <- as.numeric( tmp_df$merc1digit)
  
  tmp_df <- dplyr::full_join(tmp_df, categories_df, by=c("merc1digit"="mercatorFirstDigit") )
  tmp_df %<>% select("merc1digit", "mercMainCat", "mercatorCategory", "mercatorColours")
  
  freq_table_df <- as.data.frame(table(tmp_df$mercMainCat))
  
  # debug
  if(ret_raw_freq){ return(freq_table_df) }
  if(ret_tmp_df){ return(tmp_df) }
  
  
  # avoid warnings
  tmp_df$mercMainCat <- as.character(tmp_df$mercMainCat)
  freq_table_df$Var1 <- as.character(freq_table_df$Var1)
  
  freq_table_df <- left_join(unique(tmp_df), freq_table_df, by=c("mercMainCat"="Var1"))
  freq_table_df <- freq_table_df[order(freq_table_df$merc1digit),]
  
  if(report_empty_categ){
    freq_table_df %<>% mutate_all(~replace(., is.na(.), 0)) # mark empty categories as 0
  }
  freq_table_df %<>% drop_na()
  
  sum_freq <- sum(freq_table_df$Freq)
  freq_table_df$Percent <- freq_table_df$Freq/sum_freq*100
  sum_counts_annotated <- sum(freq_table_df$Freq[ freq_table_df$merc1digit != 35 ])
  freq_table_df$Percent_of_annotated <- freq_table_df$Freq/sum_counts_annotated*100
  
  freq_table_df %<>% select("mercatorCategory", "Freq", "Percent", "Percent_of_annotated", everything())
  
  return(freq_table_df)
}


mercatorPie <- function(categories_freq_table, russian=FALSE, predefined_colors=TRUE, axis_labels=FALSE){
  # takes freq table with mercator categories
  # returns barplot in polar coordinates
  
  require(ggplot2)
  require(RColorBrewer)
  
  if (predefined_colors) {
    categories_freq_table$mercatorColours <- as.character(categories_freq_table$mercatorColours)
    # color is saved as factor. To use text-encoded colors coersion to character is necessary
  }
  
  plot <- ggplot(categories_freq_table)+
    geom_bar(aes(x="", 
                 y=Freq,
                 fill=mercatorCategory
    ), stat = "identity"
    #, position ="dodge"
    )+
    
    scale_fill_manual(values = setNames(object = categories_freq_table$mercatorColours,
                                        nm = categories_freq_table$mercatorCategory)) +
    coord_polar("y")
  
  if(axis_labels){
    plot <- plot + geom_text(aes(x="", y=Freq, fill=mercatorCategory, label=Freq, color="white"), 
                             position = position_stack(vjust = 0.1))
  }
  print(plot)
}


# mercator top categories
mercatorBar <- function(categories_freq_table, russian=FALSE, predefined_colors=TRUE, axis_labels=FALSE){
  # takes freq table with mercator categories
  # returns barplot in polar coordinates
  require(ggplot2)
  require(RColorBrewer)
  
  if (predefined_colors) {
    categories_freq_table$mercatorColours <- as.character(categories_freq_table$mercatorColours)
    # color is saved as factor. To use text-encoded colors coersion to character is necessary
  }
  
  plot <- ggplot(categories_freq_table)+
    geom_bar(aes(x="", 
                 y=Freq,
                 fill=mercatorCategory
    ), stat = "identity"
    , position ="dodge"
    )+
    
    scale_fill_manual(values = setNames(object = categories_freq_table$mercatorColours,
                                        nm = categories_freq_table$mercatorCategory)) #+
  #coord_polar("y")
  
  if(axis_labels){
    plot <- plot + geom_text(aes(x="", y=Freq, fill=mercatorCategory, label=Freq, color="white"), 
                             position = position_stack(vjust = 0.1))
  }
  print(plot)
}


# Annotation stats
print_annot_stats <- function(annotated_df) {
  mercator_annot_count <- nrow(
    annotated_df %>% dplyr::filter(Mercator_Bincode != "35.2")
  )
  mercator_assigned_count <- nrow(
    annotated_df %>% dplyr::filter(Mercator_Bincode != "35.2" & Mercator_Bincode != "35.1")
  )
  eggnog_annotated <- nrow(
    annotated_df %>% dplyr::filter( !is.na(seed_eggNOG_ortholog) )
  )
  has_go_terms <- nrow(
    annotated_df %>% dplyr::filter( !is.na(Gene_Ontology_terms) )
  )
  has_no_mercator_nor_eggnog <- nrow(
    annotated_df %>% dplyr::filter( Mercator_Bincode == "35.2" & is.na(seed_eggNOG_ortholog) )
  )
  has_mirna_cleavage_site <- nrow(
    annotated_df %>% dplyr::filter( !is.na(SiteID) )
  )
  total_n_rows <- nrow(annotated_df)
  unique_n_rows <- nrow(dplyr::distinct(annotated_df, transcriptID))
  
  res_matrix <- matrix(
    c(total_n_rows, 100,
      unique_n_rows, unique_n_rows/total_n_rows*100,
      mercator_annot_count, mercator_annot_count/total_n_rows*100,
      mercator_assigned_count, mercator_assigned_count/total_n_rows*100,
      eggnog_annotated,  eggnog_annotated/total_n_rows*100,
      has_go_terms,  has_go_terms/total_n_rows*100,
      has_no_mercator_nor_eggnog, has_no_mercator_nor_eggnog/total_n_rows*100,
      has_mirna_cleavage_site, has_mirna_cleavage_site/total_n_rows*100
    ),
    nrow = 8,
    byrow = T
  )
  
  rownames(res_matrix) <- c("Total # of transcripts", "Unique tr IDs", "Mercator annotated", "Mercator categorised", 
                            "eggNOG annotated", "Has GO terms", "No eggnog, no mercator", "Has miRNA cleavage site")
  colnames(res_matrix) <- c("Count", "Percent")
  
  
  return(res_matrix)
}


# print counts boxplot
print_pretty_counts <- function(gene_IDs_vector,
                                d=dds,
                                deseq_normalized=TRUE,
                                annot_source=annot_for_genome_based) {
  
  require(ggbeeswarm)
  
  for (i in gene_IDs_vector) {
    geneCounts <- plotCounts(d, gene = i, intgroup = c("group"),
                             returnData = TRUE, normalized = deseq_normalized)
    plot <- ggplot(geneCounts, aes(x = group, y = count, color = group)) +
      
      ggtitle(i)+
      labs(subtitle = annot_source %>% filter(transcriptID == i) %>% select(Mercator_Name) %>% as.character(),
           caption = annot_source %>% filter(transcriptID == i) %>% select(Mercator_Description) %>% as.character())+
      geom_boxplot()+
      theme_bw()+
      geom_beeswarm(cex = 6, groupOnX=FALSE)+
      theme(axis.text.x = element_blank())
    
    #element_text(angle = 90, hjust = 1))
    print(plot)
  } 
}

print_counts_boxplots <- function(DEresult,
                                  d=dds,
                                  deseq_normalized=TRUE,
                                  samplesize=num_examples,
                                  annot_source=annot_for_genome_based,
                                  al=alpha) {
  # d stands for DESeq Dataset
  
  # avoid deseq datatypes
  current_DESeq_result <- subset(as.data.frame(DEresult), !is.na(padj) & padj < al )
  transcr_ids <- row.names(current_DESeq_result)
  
  # check if sample size is bigger than amount of genes/transcripts
  if (length(transcr_ids)<samplesize) {
    samplesize <- length(transcr_ids)
  }
  
  transcr_ids %>% sample(samplesize) %>% print_pretty_counts(d=d, deseq_normalized=deseq_normalized, annot_source=annot_source)
}