# This is a script to filter the genome collection obtained from SPADES based on the genome completeness and 
# contamination score after applying CheckM software to each assembly

#### Set directory and packages ####
setwd("Path_to_your_directory")
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)

#### Reading CheckM output file ####

genome_quality <- fread("TCS_genome_quality.csv")

# We assign quality scores as specified in Table 3 of Parks et al. 2015 (DOI: 10.1101/gr.186072.114). 
# Given that the genomes were obtained from isolate populations, their completeness levels are all very high (Near Completeness),
# and thus we did not use Completeness as a quality filter

genome_quality <- genome_quality %>% mutate("Contamination level"  = factor(case_when(Contamination == 0 ~ "No detectable contamination",
                                                                     Contamination > 0 & Contamination <= 1 ~ "Very low",
                                                                     Contamination > 1 & Contamination <= 5 ~ "Low",
                                                                     Contamination > 5 & Contamination <= 10 ~ "Medium",
                                                                     Contamination > 10 & Contamination <= 15 ~ "High",
                                                                     Contamination > 15 ~ "Very high")))

genome_quality$`Contamination level` <- factor(genome_quality$`Contamination level`,levels =
                                                 c("No detectable contamination","Very low","Low","Medium","High","Very high"))

# We obtain the sequence ID of each isolate from its Kbase identifier 
genome_quality <- genome_quality %>% 
  separate(`Bin Name`, into = c("seqID", "kbase_file"), sep = "_output__")

# We plot completeness against contamination level for the unfiltered genome database
genome_quality %>% ggplot(aes(x = Completeness, y = Contamination)) +
  geom_hline(yintercept = 5, linetype = "dashed", linewidth = 0.5, color = "grey30")+
  #geom_vline(xintercept = 95, linetype = "dashed", linewidth = 0.5, color = "grey30")+
  geom_point(aes(color = `Contamination level`),size = 4, alpha = 0.6)+
  scale_color_viridis_d(option = "plasma")+
  theme_classic()+
  theme_classic(base_size = 14,base_family = "Futura Bk BT")+
  theme(
    panel.background = element_rect(color = "black",fill='grey70'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'),
    legend.position = "bottom")

#### Get taxonomy assignment ####
# We load the taxonomy output file from the Kbase pipeline 
taxonomy_file_kbase <- fread("taxonomy_file_kbase.csv")

# Then we merge the taxonomy file with the genome quality database

# Kbase_taxonomy_quality <- left_join(genome_quality,taxonomy_file_kbase, by = "seqID")
Kbase_taxonomy_quality <- full_join(taxonomy_file_kbase,genome_quality, by = "seqID")

# The whole genome database taxonomic assignment has some inconsistencies due to changes in clade nomenclature. 
# We manually replace obsolete names
Kbase_taxonomy_quality$Genus <- gsub("Pseudomonas_A", "Stutzerimonas", Kbase_taxonomy_quality$Genus)
Kbase_taxonomy_quality$Phylum <- gsub("Proteobacteria", "Pseudomonadota", Kbase_taxonomy_quality$Phylum)

# We then filter out the genomes with high contamination levels
Kbase_taxonomy_quality_filtered <- Kbase_taxonomy_quality %>% filter(`Contamination level` %in% c("No detectable contamination","Very low","Low"))


Kbase_taxonomy_quality_filtered %>% ggplot(aes(x = Family,fill = Genus))+
  geom_bar(color = "black")+
  theme_classic(base_size = 14,base_family = "Futura Bk BT")+
  scale_fill_viridis_d(option = "plasma")+
  theme(
    panel.background = element_rect(color = "black",fill='grey70'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'),
    legend.position = "bottom")

#### Getting the unfiltered taxonomy assignment keeping unique IDs ####
# We obtain a processed taxonomy file without removing contaminated isolates, 
# but keeping only one entry by isolate (some genomes were sequenced more than once)
Kbase_taxonomy_quality_unfiltered <- Kbase_taxonomy_quality %>%
  # Create a base numeric column to group by the numeric part of seqID
  mutate(base_seqID = str_extract(seqID, "^\\d+")) %>%
  # Group by the base numeric part
  group_by(base_seqID) %>%
  # Filter rows to retain the one with the minimum Contamination value
  filter(Contamination == min(Contamination)) %>%
  # In case of ties, keep only the first occurrence
  slice_head(n = 1) %>%
  # Remove "reseq" from seqID when present
  mutate(seqID = str_remove(seqID, "reseq$")) %>%
  ungroup() %>%
  dplyr::select(-base_seqID)

#### Getting a filtered taxonomy assignment keeping unique IDs ####
# We obtain a processed taxonomy file removing contaminated isolates, 
# and keeping only one entry by isolate (some genomes were sequenced more than once)
Kbase_taxonomy_quality_filtered <- Kbase_taxonomy_quality_filtered %>%
  # Create a base numeric column to group by the numeric part of seqID
  mutate(base_seqID = str_extract(seqID, "^\\d+")) %>%
  # Group by the base numeric part
  group_by(base_seqID) %>%
  # Filter rows to retain the one with the minimum Contamination value
  filter(Contamination == min(Contamination)) %>%
  # In case of ties, keep only the first occurrence
  slice_head(n = 1) %>%
  # Remove "reseq" from seqID when present
  mutate(seqID = str_remove(seqID, "reseq$")) %>%
  ungroup() %>%
  dplyr::select(-base_seqID)


write.csv(Kbase_taxonomy_quality_unfiltered,"Kbase_taxonomy_quality_unfiltered.csv",row.names = FALSE)

write.csv(Kbase_taxonomy_quality_filtered,"Kbase_taxonomy_quality_filtered.csv",row.names = FALSE)

#### Making the assembly quality data obtained from Quast compatible with the taxonomic assignment ####

# Here we obtain a list of isolates from the taxonomic assignment file
Kbase_taxonomy_quality_unfiltered_original_id <- Kbase_taxonomy_quality %>%
  mutate(base_seqID = str_extract(seqID, "^\\d+")) %>%
  group_by(base_seqID) %>%
  filter(Contamination == min(Contamination)) %>%
  slice_head(n = 1)%>%
  mutate(seqID = str_replace(seqID, "reseq","_reseq"))%>%
  mutate(seqID = str_replace(seqID, "335_reseq","335")) # This specific isolate is otherwise not found because of a difference in nomenclature

# We then load the original quast assembly quality file
quast_stats_unfiltered <- fread("quast_stats_unfiltered.csv")

# And we filter it based on the list of isolates 
quast_stats_filtered <- quast_stats_unfiltered %>% 
  filter(curve_id %in% Kbase_taxonomy_quality_unfiltered_original_id$seqID) %>%
  dplyr::select(c("curve_id","GC%","Total length")) %>%
  mutate(curve_id = str_replace(curve_id, "275_reseq","275")) # This specific isolate is otherwise not found because of a difference in nomenclature. The genomic info of both sequencing jobs was identical

# We finally write the output file that contains basic genome information
write.csv(quast_stats_filtered, "quast_stats_filtered.csv", row.names = FALSE)
