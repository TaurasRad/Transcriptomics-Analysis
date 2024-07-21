library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(tidyr)
library(openxlsx)
library(clusterProfiler)
library(patchwork)
setwd("/Users/vardaspavarde/Desktop/")


gene_annotations_file <- import('/Users/vardaspavarde/Desktop/Pseudomonas_aeruginosa_UCBPP-PA14_109.gtf')


read_and_filter_bed <- function(file_path, score_threshold, percentage_threshold) {
  # Column names for BED files
  bed_column_names <- c(
    "chrom", "chromStart", "chromEnd", "name", "score",
    "strand", "thickStart", "thickEnd", "itemRgb",
    "blockCount", "percentage", "blockStarts"
  )
  
  # Read and process the .bed file
  bed_df <- read.table(file_path, header = FALSE, sep = "\t", fill = TRUE, check.names = FALSE)
  colnames(bed_df) <- bed_column_names[1:min(ncol(bed_df), length(bed_column_names))]

  na_or_empty_names <- which(is.na(names(bed_df)) | names(bed_df) == "")
  if(length(na_or_empty_names) > 0) {
    names(bed_df)[na_or_empty_names] <- paste0("X", na_or_empty_names)
  }

  filtered_df <- bed_df %>%
    filter(!is.na(score) & score > score_threshold,
           !is.na(percentage) & percentage > percentage_threshold)
  
  return(filtered_df)
}

# File paths
file_paths <- c("/Users/vardaspavarde/Desktop/PA14_plus_DM3VIR/day0.bed",
                "/Users/vardaspavarde/Desktop/PA14_plus_DM3VIR/day01.bed",
                "/Users/vardaspavarde/Desktop/PA14_plus_DM3VIR/day02.bed",
                "/Users/vardaspavarde/Desktop/PA14_plus_DM3VIR/day03.bed",
                "/Users/vardaspavarde/Desktop/PA14_plus_DM3VIR/day04.bed",
                "/Users/vardaspavarde/Desktop/PA14_plus_DM3VIR/day05.bed")
score_threshold <- 900
percentage_threshold <- 50

filtered_beds <- lapply(file_paths, read_and_filter_bed, score_threshold, percentage_threshold)

summarize_and_compare <- function(filtered_bed, original_bed_path) {
  original_bed <- read.table(original_bed_path, header = FALSE, sep = "\t")
  
  # Summary statistics
  summary_stats <- filtered_bed %>%
    summarise(
      count_filtered = n(),
      mean_score = mean(score, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE),
      mean_percentage = mean(percentage, na.rm = TRUE),
      median_percentage = median(percentage, na.rm = TRUE)
    )
  
  # Comparative Analysis
  total_count <- nrow(original_bed)
  filtered_count <- nrow(filtered_bed)
  proportion <- filtered_count / total_count
  
  # Combine summary and comparative stats
  result_table <- cbind(summary_stats, total_count = total_count, 
                        filtered_count = filtered_count, proportion = proportion)
  return(result_table)
}


all_results <- lapply(seq_along(filtered_beds), function(i) {
  summarize_and_compare(filtered_beds[[i]], file_paths[i])
})

final_summary_table <- bind_rows(all_results)

#SUPER CIA VEIKIA PUIKIAI !!!!!!!!!!!!!!



# Filter each filtered list once more to get rid of redundant information
filtered_beds_raw <- filtered_beds
filtered_beds <- lapply(filtered_beds, function(df) {
  select(df, all_of(c("chromStart", "chromEnd", "score", "percentage")))
})
#THIS IS MUCH MORE CLEAR AND FILTERED 
merged_bed <- Reduce(function(x, y) merge(x, y, by = c("chromStart", "chromEnd"), all = TRUE), filtered_beds)

### Ziurim ar veikia- veikia

#THIS IS MANUAL
names_vector <- c(
  "chromStart", "chromEnd", 
  "score.x", "percentage.x", 
  "score.y", "percentage.y", 
  "score.x.x", "percentage.x.x", 
  "score.y.y", "percentage.y.y", 
  "score.x.x.x", "percentage.x.x.x", 
  "score.y.y.y", "percentage.y.y.y"
)
names(merged_bed) <- names_vector
identify_methylated_days <- function(row) {
  days_methylated <- c()
  if (!is.na(row$score.x)) days_methylated <- c(days_methylated, "0") #Suziureti eiles tvarka cia
  if (!is.na(row$score.y)) days_methylated <- c(days_methylated, "1")
  if (!is.na(row$score.x.x)) days_methylated <- c(days_methylated, "2")
  if (!is.na(row$score.y.y)) days_methylated <- c(days_methylated, "3")
  if (!is.na(row$score.x.x.x)) days_methylated <- c(days_methylated, "4")
  if (!is.na(row$score.y.y.y)) days_methylated <- c(days_methylated, "5")####SUZIURET ar tie failai prisikirti
  #DUOMENIS NES ZIAURIAI KEISTAI RODO CIA
  return(paste(days_methylated, collapse = ","))
}

# apply every row
merged_bed$days_methylated <- sapply(seq_len(nrow(merged_bed)), function(i) identify_methylated_days(merged_bed[i, ]))

days_table <- merged_bed %>%
  select(chromStart, days_methylated) %>%
  distinct(chromStart, days_methylated)

#CIA JAU DATA ANALIZE IR GRAFIKAI
#print(days_table)
#this is a very terrible way to show how methylation changes distributes
  # pie_chart <- ggplot(day_counts, aes(x = "", y = n, fill = days_methylated)) +
  #   geom_bar(stat = "identity", width = 1, color = "white") +
  #   coord_polar("y", start = 0) +
  #   theme_void() +
  #   theme(legend.position = "bottom") +  
  #   labs(title = "Days Methylated Distribution", fill = "Days Methylated") +
  #   guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))  
  # print(pie_chart)

day_counts <- count(days_table, days_methylated)
column_chart <- ggplot(day_counts, aes(x = days_methylated, y = n, fill = days_methylated)) +
  geom_col() + 
  theme_minimal() + 
  labs(title = "Days Methylated Distribution", 
       x = "Days Methylated", 
       y = "Count", 
       fill = "Days Methylated") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14), # Bigger x-axis labels
        axis.text.y = element_text(size = 14),    # Bigger y-axis labels
        axis.title.x = element_text(size = 16),   # Bigger x-axis title
        axis.title.y = element_text(size = 16),   # Bigger y-axis title
        plot.title = element_text(size = 20),     # Bigger plot title
        legend.title = element_text(size = 14),   # Bigger legend title
        legend.text = element_text(size = 12))    # Bigger legend text

pdf(file = "column_chart_new.pdf", width = 8, height = 8)
plot(column_chart)
dev.off()
print(column_chart)
######### THIS WORKS GOOOOOD

#I would guess that only the first day was when the spacer is held - all other days are RM system regulated 


#CIA BUS MAPPINAMI GENAI IS VISU DIENU
#######GERESNE VERSIJA :
# Function to find overlaps and count methylated points
  process_bed_data <- function(bed_data, gene_annotations_file) {
    bed_gr <- GRanges(
      seqnames = bed_data$chrom,
      ranges = IRanges(start = bed_data$chromStart, end = bed_data$chromEnd),
      strand = bed_data$strand
    )
    
    # One chromosome 
    seqlevels(gene_annotations_file) <- "1"
    seqlevels(bed_gr) <- "1"
    
    overlaps <- findOverlaps(bed_gr, gene_annotations_file)
    # Count overlaps 
    overlap_counts <- table(subjectHits(overlaps))
    
    gene_ids <- mcols(gene_annotations_file)$gene_id[as.integer(names(overlap_counts))]
    gene_names <- mcols(gene_annotations_file)$name[as.integer(names(overlap_counts))]
    locus_tags <- mcols(gene_annotations_file)$locus_tag[as.integer(names(overlap_counts))]
    transcript_ids <- mcols(gene_annotations_file)$transcript_id[as.integer(names(overlap_counts))]
    
    # Create a data frame with the new columns
    methylated_points_per_gene <- data.frame(
      gene_id = gene_ids, 
      gene_name = gene_names, 
      locus_tag = locus_tags,
      transcript_id = transcript_ids,
      count = as.integer(overlap_counts)
    )
    
    return(methylated_points_per_gene)
  }
  
  # bed data process
  gene_count_list <- lapply(seq_along(filtered_beds_raw), function(i) {
    bed_data <- filtered_beds_raw[[i]]
    methylated_points_per_gene <- process_bed_data(bed_data, gene_annotations_file)
    assign(paste0("methylated_points_per_gene_", i), methylated_points_per_gene, envir = .GlobalEnv)
    return(methylated_points_per_gene)
  })
  
topmeth <- lapply(gene_count_list, function(x) head(x[order(x[,3], decreasing = F),]))
#   

# Sitas PGD1656278 daugiausiai pasikartoja 
################Va cia eksperimentuojam su paperiais
# library(rentrez)
# 
# # Function to search PubMed for articles referencing a gene
# get_papers_for_gene <- function(locus_tag) {
#   search_results <- entrez_search(db="pubmed", term=paste(locus_tag, "AND Pseudomonas aeruginosa PA14[Title/Abstract]"))
#   pmids <- search_results$ids
#   return(pmids)
# }
# 
# # Assuming 'gene_count_list' is your list of data frames
# # These are only for the papers for top in the list gene's
# # pasiziureta ar naudoti gene id ar kazka kita
# papers_list <- lapply(topmeth, function(df) {
#   df$papers <- sapply(df$locus_tag, get_papers_for_gene)
#   return(df)
# })



#Paprasciau bandyt cia isdropint tuscius laukelius
names(gene_count_list) <- paste0("day_", seq_along(gene_count_list) - 1, "_count")

gene_count_list <- Filter(function(df) nrow(df) > 0, gene_count_list)
#### TIK CIA RISKY SUZIURET AR VISKAS GERAI NES LIKS TRYS DIENOS
###### IKICia LOADINu


for (i in seq_along(gene_count_list)) {
  gene_count_list[[i]]$day <- names(gene_count_list)[i]
}

gene_count_list_wide <- lapply(gene_count_list, function(df) {
  pivot_wider(df, names_from = day, values_from = count, values_fill = list(count = 0))
})

merged_gene_counts <- gene_count_list_wide[[1]]
for(i in 2:length(gene_count_list_wide)) {
  merged_gene_counts <- full_join(merged_gene_counts, gene_count_list_wide[[i]], 
                                  by = c("gene_id", "gene_name", "locus_tag", "transcript_id"))
}

merged_gene_counts <- merged_gene_counts %>%
  mutate_all(~replace_na(., 0))

# IKI CIA VISKAS SUPER IR OK

# # Eksportavau duomenis
# library(writexl)
# write_xlsx(merged_gene_counts, path = "~/Desktop/5day_gene_counts.xlsx")
# 


# Nuo cia eksperimentas: 
#man reikia nefiltruotu duomenu

# file_path <- "/Users/vardaspavarde/Desktop/PA14_plus_DM3VIR/day0.bed"
# 
# bed_column_names <- c(
#   "chrom", "chromStart", "chromEnd", "name", "score",
#   "strand", "thickStart", "thickEnd", "itemRgb",
#   "blockCount", "percentage", "blockStarts",'X1','X2','X3','X4'
# )
# 
# # Read and process the .bed file
# bedas <- read.table(file_path, header = FALSE, sep = "\t", fill = TRUE, check.names = FALSE)
# colnames(bedas) <- bed_column_names[1:min(ncol(bedas), length(bed_column_names))]

######################
# 
# library(dunn.test)
# 
# max_end <- max(bedas$chromEnd)
# 
# # Generate a complete sequence of positions
# complete_seq <- data.frame(chromStart = full_range[1]:(full_range[2]-1),
#                            chromEnd = (full_range[1]+1):full_range[2])
# 
# template_df <- data.frame(chrom = 'contig_1',chromStart = min_start:(max_end - 1),chromEnd = (min_start + 1):max_end, name = 'Fake', score = 0, strand = '+', thickStart = complete_seq$chromStart, thickEnd = complete_seq$chromEnd, itemRgb = '0,0,0', blockCount = 0,  percentage = 0.00,  blockStarts = 0, X1 = 0, X2 = 0, X3 = 0)
# # Assuming 'methylation_data' is prepared with columns: 'chromosome', 'start', 'end', 'methylation_level'
# 
# 
# bedas2 <- data.frame(chromStart=bedas$chromStart, chromEnd=bedas$chromEnd,name=bedas$name)
# template2 <- data.frame(chromStart=template_df$chromStart, chromEnd=template_df$chromEnd,name=template_df$name)
# 
# 
# chimera <- left_join(template2,bedas2, by = c('chromStart','chromEnd'))
# chimera <- chimera %>% select(-name.x)
# chimera$name.y[is.na(chimera$name.y)] <- "Fake"
# # cia jau calculiavimas, trying at 1000
# num_of_intervals <- round(nrow(chimera)/1000,0)
# 
# chimera$chromEnd <- NULL
# 
# chimera$interval <- (chimera$chromStart - 1) %/% 1000
# 
# methylation_summary <- chimera %>%
#   group_by(interval) %>%
#   summarise(
#     methylated_count = sum(name.y == "5mC"),
#     total_count = n(),
#     percentage_methylated = (methylated_count / total_count) * 100
#   ) %>%
#   ungroup()
# 
# 
# methylation_summary$interval_range <- paste(methylation_summary$interval * 1000, (methylation_summary$interval + 1) * 1000 - 1, sep = "-")
# 
# #####
# final_meth <- methylation_summary %>%
#   select(interval_range, percentage_methylated)
# # siaip gaunam tobulai ta lauka kur reikia wow net
# final_output_sorted <- final_output %>%
#   arrange(desc(percentage_methylated))
# 




# this is bad 
gr <- gene_annotations_file
gene_ids <- mcols(gr)$gene_id
starts <- start(gr)
ends <- end(gr)
chromosomes <- seqnames(gr)

#patogi forma
gene_ranges <- data.frame(
  gene_id = gene_ids,
  start = starts,
  end = ends
)
genome_len <- 6527555
# 
# gene - category
# nr of methylated bases        - nr of diffrentialy expressed genes
# q - is the number of methylated - 1.    cick_gene$day_0_count 
# m- nr of methylated loci in the entire genome.   sum(cick_gene$day_0_count) 
# n - unmethylated points in the entire genome.    genome_len - sum(cick_gene$day_0_count)
# k - is the gene size - nucleotides.       
# m and n  are always the same


#Intervals
merged_gene_counts_intervals <- full_join(gene_ranges,merged_gene_counts, by = "gene_id")
merged_gene_counts_intervals$transcript_id <- NULL

cick_gene<- na.omit(merged_gene_counts_intervals)

######

for (day_name in names(gene_count_list)) {
  day_count_col <- day_name
  p_value_col <- gsub("count", "p_value", day_name)
  
  total_day_count <- sum(cick_gene[[day_count_col]], na.rm = TRUE)
  
  cick_gene <- cick_gene %>%
    rowwise() %>%
    mutate(
      !!p_value_col := phyper(
        get(day_count_col) - 1,
        total_day_count, 
        genome_len - total_day_count,  
        get(day_count_col), 
        lower.tail = FALSE
      )
    ) %>%
    ungroup()  
}


p_value_cols <- grep("p_value", names(cick_gene), value = TRUE)
cick_gene <- cick_gene %>%
  mutate(across(all_of(p_value_cols), ~p.adjust(.x, method = 'BH')))


cick_gene_significant <- cick_gene %>%
  mutate(across(starts_with("p_value_day_"), ~if_else(.x > 0.05, "non-significant", as.character(.x))))
##### valio valiovalio


# so leave the full scientific p values:
# and callulate the p adjust when it doesnt have the non-significant genes
#
#  library(writexl)
#  write_xlsx(cick_gene_significant, path = "~/Desktop/significant_genes_850_85.xlsx")
# # # 
####### VISKAS IKI CIA VISKAS VEIKIA - CIA SUBPLOTAM BUS VIETa


# top methylation grafikai
topmeth <- lapply(gene_count_list, function(x) head(x[order(x[,3], decreasing = T),]))
topmeth <- lapply(topmeth, function(df) {
  select(df, all_of(c("locus_tag", "gene_name", "count")))
})

theme_thesis <- function() {
  theme_minimal(base_size = 14) + 
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(size = 0.5, color = "grey80"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
}

# Loop over each item in the list 'topmeth'
plot_list <- list()
for (i in seq_along(topmeth)) {
  day_label <- names(topmeth)[i]  # Get the name of the current list element
  plot_list[[day_label]] <- ggplot(topmeth[[i]], aes(x = gene_name, y = count, fill = as.factor(i))) +
    geom_col(width = 0.8) +  # Increase the width of the columns
    geom_text(aes(label = count), vjust = -0.3, color = "black", size = 3.5) +  # Add text on top of the columns
    theme_minimal() + 
    labs(title = paste("Top 5 Methylated Genes -", day_label),  # Use the day_label for title
         x = "Gene Function", 
         y = "Count") +
    theme_thesis()
}
#############GRAFIKAM DEMO
# for (i in seq_along(plot_list)) {
#   plot_name <- names(plot_list)[i]
#   file_name <- paste("plot_", plot_name, ".pdf", sep = "")
#   ggsave(file_name, plot = plot_list[[i]], device = "pdf", width = 11, height = 8.5)
# }


# TOLIAU GRAFIKAI LENTELIu
# tabletes5 <- methylated_points_per_gene_6 %>%
#   select(gene_id,gene_name ,count) %>%
#   arrange(desc(count)) %>%  # Arrange by column3 in descending order
#   slice_head(n = 10) 
# 
# # Saving the data frame df as a CSV file
# write.csv(methylated_points_per_gene_6, "sesta.csv", row.names = FALSE)
# 





###Surisimas su TF

# 
# cick_gene <- cick_gene %>%
#   left_join(aeruginosa %>% select(Locus.Tag, Gene.Name), by = c("locus_tag" = "Locus.Tag"))
# 
# methylated_points_per_gene_3 <- methylated_points_per_gene_3 %>%
#   left_join(aeruginosa %>% select(Locus.Tag, Gene.Name), by = c("locus_tag" = "Locus.Tag"))
# 
# methylated_points_per_gene_4 <- methylated_points_per_gene_4 %>%
#   left_join(aeruginosa %>% select(Locus.Tag, Gene.Name), by = c("locus_tag" = "Locus.Tag"))
# 
# methylated_points_per_gene_5 <- methylated_points_per_gene_5 %>%
#   left_join(aeruginosa %>% select(Locus.Tag, Gene.Name), by = c("locus_tag" = "Locus.Tag"))
# 
# methylated_points_per_gene_6 <- methylated_points_per_gene_6 %>%
#   left_join(aeruginosa %>% select(Locus.Tag, Gene.Name), by = c("locus_tag" = "Locus.Tag"))
# 
# 
# # Cia netobulai bet pazurekime ar veikia
# 
# methylated_points_per_gene_3 <- methylated_points_per_gene_3 %>%
#   left_join(tf_data%>% select(mx, label,type, group), by = c("Gene.Name" = "label"))
# 
# methylated_points_per_gene_4 <- methylated_points_per_gene_4 %>%
#   left_join(tf_data%>% select(mx, label,type, group), by = c("Gene.Name" = "label"))
# 
# methylated_points_per_gene_5 <- methylated_points_per_gene_5 %>%
#   left_join(tf_data%>% select(mx, label,type, group), by = c("Gene.Name" = "label"))
# 
# methylated_points_per_gene_6 <- methylated_points_per_gene_6 %>%
#   left_join(tf_data%>% select(mx, label,type, group), by = c("Gene.Name" = "label"))
# 
# ###### Prafiltruoti be NA
# 
# m_gene_path_3 <- na.omit(methylated_points_per_gene_3)
# m_gene_path_3$mx <- NULL
# m_gene_path_3<- m_gene_path_3[m_gene_path_3$group != 'undefined', ]
# 
# m_gene_path_4 <- na.omit(methylated_points_per_gene_4)
# m_gene_path_4$mx <- NULL
# m_gene_path_4<- m_gene_path_4[m_gene_path_4$group != 'undefined', ]
# 
# m_gene_path_5 <- na.omit(methylated_points_per_gene_5)
# m_gene_path_5$mx <- NULL
# m_gene_path_5<- m_gene_path_5[m_gene_path_5$group != 'undefined', ]
# 
# 
# m_gene_path_6 <- na.omit(methylated_points_per_gene_6)
# m_gene_path_6$mx <- NULL
# m_gene_path_6<- m_gene_path_6[m_gene_path_6$group != 'undefined', ]
# ######
# # write.csv(m_gene_path_3, "2_day_associated_genes.csv", row.names = FALSE)
# # write.csv(m_gene_path_4, "3_day_associated_genes.csv", row.names = FALSE)
# # write.csv(m_gene_path_5, "4_day_associated_genes.csv", row.names = FALSE)
# # write.csv(m_gene_path_6, "5_associated_genes.csv", row.names = FALSE)
# 
# ###GRAAPHS


topmeth <- lapply(names(topmeth), function(name) {
  df <- topmeth[[name]]
  df$table_name <- gsub("_count", "", name)  # Remove '_count' from name
  return(df)
})

# Combine the list into a single data frame
combined_topmeth <- bind_rows(topmeth)

# Create the plot
nr_genai <- ggplot(combined_topmeth, aes(x = gene_name, y = count)) +
  geom_col(fill = "cornflowerblue", color = "black") +
  geom_text(aes(label = count), vjust = -0.5, size = 4) +
  labs(
    title = "Top 10 Most Methylated Genes",
    x = "Gene Name",
    y = "Number of Methylated Loci"
  ) +
  theme_thesis() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, NA)) +  # Ensure y-axis is fixed to integers
  facet_wrap(~ table_name, scales = "free_x")  # Use free scales for x-axis

# Display the plot
print(nr_genai)


# topmeth <- lapply(names(topmeth), function(name) {
#   df <- topmeth[[name]]
#   df$table_name <- gsub("_count", "", name)  # Remove '_count' from name
#   return(df)
# })

# Combine the list into a single data frame
Day2_top <- topmeth[[1]]

#DAY 2 and Day 5 plots
nr_genai <- ggplot(combined_topmeth, aes(x = gene_name, y = count)) +
  geom_col(fill = "cornflowerblue", color = "black") +
  geom_text(aes(label = count), vjust = -0.5, size = 4) +
  labs(
    title = "Top 10 Most Methylated Genes",
    x = "Gene Name",
    y = "Number of Methylated Loci"
  ) +
  theme_thesis() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, NA)) +  # Ensure y-axis is fixed to integers
  facet_wrap(~ table_name, scales = "free_x")  # Use free scales for x-axis

# Display the plot
print(nr_genai)