library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(tidyr)
library(openxlsx)
library(clusterProfiler)
library(patchwork)
library(writexl)
setwd('/Users/vardaspavarde/Documents/GitHub/Transcriptomics-Analysis')


gene_annotations_file <- import('data/Pseudomonas_aeruginosa_UCBPP-PA14_109.gtf')

read_and_filter_bed <- function(file_path, score_threshold, percentage_threshold) {
  # Column names for BED files
  bed_column_names <- c(
    "chrom", "chromStart", "chromEnd", "name", "score",
    "strand", "thickStart", "thickEnd", "itemRgb",
    "blockCount", "percentage", "blockStarts"
  )
  # Extract the day information from the file name (assuming the day is a number in the file name)
  day <- gsub("\\D", "", basename(file_path))
  # Read and process the .bed file
  bed_df <- read.table(file_path, header = FALSE, sep = "\t", fill = TRUE, check.names = FALSE)
  colnames(bed_df) <- bed_column_names[1:min(ncol(bed_df), length(bed_column_names))]
  # Handle columns with NA or empty names
  na_or_empty_names <- which(is.na(names(bed_df)) | names(bed_df) == "")
  if(length(na_or_empty_names) > 0) {
    names(bed_df)[na_or_empty_names] <- paste0("X", na_or_empty_names)
  }
  # Filter the data frame
  filtered_df <- bed_df %>%
    filter(!is.na(score) & score > score_threshold,
           !is.na(percentage) & percentage > percentage_threshold)
  # Add the day column
  filtered_df$day <- day
  return(filtered_df)
}
# File paths to bed files (They must be loaded in the same order as the experiment took place)
file_paths <- c("data/day0.bed",
                "data/day01.bed",
                "data/day02.bed",
                "data/day03.bed",
                "data/day04.bed",
                "data/day05.bed")
#Score and percentage threshold for filtering 
score_threshold <- 500
percentage_threshold <- 50
#resulting beds
filtered_beds <- lapply(file_paths, read_and_filter_bed, score_threshold, percentage_threshold)
# function for getting summary of stats for the filtered datasets
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
# get stats for each day of the experiment
all_results <- lapply(seq_along(filtered_beds), function(i) {
  summarize_and_compare(filtered_beds[[i]], file_paths[i])
})
# Final stat table
final_summary_table <- bind_rows(all_results)
# Filter each filtered list once more to get rid of redundant information
merged_bed <- bind_rows(filtered_beds)
#methylation loci 
days_table <- merged_bed %>%
  select(chromStart, day) %>%
  distinct(chromStart, day)
#Methylation data plots
day_counts <- count(days_table, day)
column_chart <- ggplot(day_counts, aes(x = day, y = n, fill = day)) +
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
# path to plot folder
output_directory <- "/Users/vardaspavarde/Documents/GitHub/Transcriptomics-Analysis/plots"
output_file <- file.path(output_directory, "Methylation_distribution.pdf")
pdf(file = output_file, width = 8, height = 8)
plot(column_chart)
#######GERESNE VERSIJA :
# Function to find overlaps and count methylated points
  process_bed_data <- function(bed_data, gene_annotations_file) {
    bed_gr <- GRanges(
      seqnames = bed_data$chrom,
      ranges = IRanges(start = bed_data$chromStart, end = bed_data$chromEnd),
      strand = bed_data$strand
    )
    # One chromosome - hardcoded needs adjustment with multi_chromosome having organisms
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
  # Counts the amount of genes methylated
  gene_count_list <- lapply(seq_along(filtered_beds_raw), function(i) {
    bed_data <- filtered_beds_raw[[i]]
    methylated_points_per_gene <- process_bed_data(bed_data, gene_annotations_file)
    assign(paste0("methylated_points_per_gene_", i), methylated_points_per_gene, envir = .GlobalEnv)
    return(methylated_points_per_gene)
  })
#top methylated genes for each day 
topmeth <- lapply(gene_count_list, function(x) head(x[order(x[,3], decreasing = F),]))
#Gene count managemen
names(gene_count_list) <- paste0("day_", seq_along(gene_count_list) - 1, "_count")
gene_count_list <- Filter(function(df) nrow(df) > 0, gene_count_list)
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
# Data export to processed data folder
write_xlsx(merged_gene_counts, path = "processed_data/all_day_gene_counts.xlsx")



# CIA sustojau




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