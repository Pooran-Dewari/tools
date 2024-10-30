## genomepy

#### Search and install the genome

```
#search
genomepy search gigas

#install
genomepy install --annotation cgigas_uk_roslin_v1
```
#### Extract longest CDS per gene in bed format using a gff3 file
```
library(tidyverse)
# Aim: extract longest cds in bed format using gff3 file for gimme TFBS analysis

# read gff3 file
gff <- read_tsv("Crassostrea_gigas.cgigas_uk_roslin_v1.59.chr.gff3", comment = "#", col_names = F) 

# removing first (weird) row
gff <- gff[-1,] 

# extract gene id info into a new column
gff_v1 <- gff %>% 
  filter(X3 == "CDS") %>% 
  mutate(CDSLEN = X5 - X4) %>% 
  separate(X9, into = c("a", "gene"), sep= ":") %>% 
  separate(gene, into = c("gene_id", "version"), sep="\\.")


# select longest CDS
cds_longest <- gff_v1 %>%
  group_by(gene_id) %>%
  slice_max(CDSLEN, n = 1) %>%
  ungroup()

# transform to bed format
cgigas_cds_longest_bed <- cds_longest %>% 
  select(X1, X4, X5, gene_id, X7) %>% 
  mutate(score = 0, .after=gene_id)

# extract promoter, i.e., TSS to 200 bp upstream 
cgigas_longest_cds_prom <- cgigas_cds_longest_bed %>% 
  mutate(TSS = X4, UP_200=X4-200) %>% 
  select(CHR=X1, UP_200, TSS, gene_id, score, STRAND=X7)

# write bed file
write_tsv(cgigas_cds_longest_bed, "cgigas_cds_longest_bed.bed", col_names = F)
write_tsv(cgigas_longest_cds_prom, "cgigas_longest_cds_prom.bed", col_names = F)
```

#### get genes for each cluster in bed format for gimme motifs de novo TFBS analysis
```
# get genes upregulated in each cluster as df

all_markers_up <- read_tsv("all_markers_up_for_gimme.tsv")

#all_up_bed <- left_join(all_markers_up, cgigas_cds_longest_bed, by = c("gene" = "gene_id"))

# Assuming all_up_bed_list is a named list with clusters
all_up_bed_list <- split(all_markers_up$gene, all_markers_up$cluster)

# Convert each list in all_up_bed_list into a data frame
all_up_bed_dfs <- lapply(names(all_up_bed_list), function(cluster_name) {
  # Create a data frame with a "gene" column and a "cluster" column for each list element
  data.frame(gene = all_up_bed_list[[cluster_name]], cluster = cluster_name, stringsAsFactors = FALSE)
})

# Name each data frame "cluster_0", "cluster_1", ..., "cluster_19"
names(all_up_bed_dfs) <- paste0("cluster_", names(all_up_bed_list))

# Optionally, assign each data frame to the global environment
list2env(all_up_bed_dfs, envir = .GlobalEnv)



## now get up regulated genes for each cluster in bed format for gimme

## example for one cluster
cluster_0_prom_bed <- inner_join(cluster_0, cgigas_longest_cds_prom, by= c("gene" ="gene_id")) %>% 
  distinct() %>% 
  select(-cluster, -gene)

write_tsv(cluster_0_prom_bed, "cluster_0_prom.bed", col_names = F)


## loop to get for all clusters

# Define the number of clusters
num_clusters <- 20

# Loop through each cluster from 0 to 19
for (i in 0:(num_clusters - 1)) {
  # Dynamically get the cluster data frame (e.g., cluster_0, cluster_1, ...)
  cluster_name <- paste0("cluster_", i)
  cluster_df <- get(cluster_name)  # Retrieve the data frame by name
  
  # Perform the inner join and filtering
  cluster_prom_bed <- inner_join(cluster_df, cgigas_longest_cds_prom, by = c("gene" = "gene_id")) %>% 
    distinct() %>% 
    select(-cluster, -gene)
  
  # Write the output to a .bed file
  output_file <- paste0(cluster_name, "_prom.bed")
  write_tsv(cluster_prom_bed, output_file, col_names = FALSE)
}
```

#### use the file above for gimme motifs de novo tfbs analysis
```
gimme motifs cluster_1_prom.bed cluster_1_prom --denovo -g ../Documents/parse/Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa
```
