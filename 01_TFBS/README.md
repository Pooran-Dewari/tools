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

# write bed file
write_tsv(cgigas_cds_longest_bed, "cgigas_cds_longest_bed.bed", col_names = F)
```
