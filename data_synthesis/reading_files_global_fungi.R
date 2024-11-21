library(tidyverse)
library(readxl)

otu.table.long <- readRDS("GlobalFungi//GF5_FUNGAL_OTUS_maxeval_e-40_OTUTAB_long.rds")

meta_data <- read_excel(
  "GlobalFungi//GF5_FUNGAL_OTUS_maxeval_e-40_SAMPLES_METADATA.xlsx")

gz_file <- gzfile("GlobalFungi//GF5_FUNGAL_OTUS_maxeval_e-40//GF5_FUNGAL_OTUS_maxeval_e-40//GF5_FUNGAL_OTUS_maxeval_e-40_BLAST_UNITE9_BEST.txt.gz", "rt")
data <- read.table(gz_file, header = F)
close(gz_file)

taxonomy <- 
  data.frame(
    phylum =
      str_extract(data$V3, "p__[^;]+"),
    
    order = 
      str_extract(data$V3, "o__[^;]+"),
    
    family =
      str_extract(data$V3, "f__[^;]+"),
    
    genera =
      str_extract(data$V3, "g__[^;]+"),
    
    species =
      str_extract(data$V3, "s__[^;]+"),
    
    OTU =
      str_extract(data$V1, "OTU[0-9]+")              
  )


all_data_1 <- 

left_join(otu.table.long, meta_data[,c("PermanentID",
                                       "latitude", "longitude", "sample_type","add_date",
                                       "paper_id", "continent",
                                       "sequencing_platform", "target_gene")] %>% 
            rename(SampleID = PermanentID))

all_data_2 <- left_join(all_data_1, taxonomy, by = "OTU")

all_data_2$Database <- "Global_Fungi"

saveRDS(all_data_2, "data_synthesis//long_format//Global_Fungi_long_Carlos.RDS")

