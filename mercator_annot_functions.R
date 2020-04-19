
#load("mercatorCatsColoursDf_rus.Robject") # loads mercatorCatsColoursDf with russian translation

mercatorTopCats <- function(annotated_df, report_empty_categ=FALSE, categories_df=mercatorCatsColoursDf, 
                            ret_raw_freq=F, ret_tmp_df=F){
  # make frequency table for an annotated dataframe
  # returns dataframe
  require(dplyr)
  
  tmp_df <- annotated_df %>% separate(Mercator_Name, into = c("mercMainCat", NA),
                                      sep='[.]', remove = F, extra = "drop", fill = "right") %>% 
    separate(Mercator_Bincode, into = c("merc1digit", NA),
             sep='[.]', remove = F, extra = "drop", fill = "right") # -> tmp_df
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


print_pretty_counts <- function(gene_IDs_vector, d=dds) {
  
  require(ggbeeswarm)
  
  for (i in gene_IDs_vector) {
    geneCounts <- plotCounts(d, gene = i, intgroup = c("group"),
                             returnData = TRUE)
    plot <- ggplot(geneCounts, aes(x = group, y = count, color = group)) +
      
      ggtitle(i)+
      geom_boxplot()+
      theme_bw()+
      geom_beeswarm(cex = 6)+
      theme(axis.text.x = element_blank())
    
    #element_text(angle = 90, hjust = 1))
    print(plot)
  } 
}


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




