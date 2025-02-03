library(tidyverse)
library(readxl)

carlos_computer <- "C://Users//caaguila//Dropbox//ESIL_data//"

otu.table.long <- readRDS(paste(carlos_computer,
                                "GlobalFungi//GF5_FUNGAL_OTUS_maxeval_e-40_OTUTAB_long.rds",
                                sep = ""))

meta_data <- read_excel(paste(carlos_computer,
  "GlobalFungi//GF5_FUNGAL_OTUS_maxeval_e-40_SAMPLES_METADATA.xlsx",
  sep = ""))

gz_file <- gzfile(paste(carlos_computer,
  "GlobalFungi//GF5_FUNGAL_OTUS_maxeval_e-40//GF5_FUNGAL_OTUS_maxeval_e-40//GF5_FUNGAL_OTUS_maxeval_e-40_BLAST_UNITE9_BEST.txt.gz",
  sep = ""), "rt")
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
                                       "sequencing_platform", "target_gene", "primers")] %>% 
            rename(SampleID = PermanentID))

all_data_2 <- left_join(all_data_1, taxonomy, by = "OTU")

all_data_2$Database <- "Global_Fungi"

saveRDS(all_data_2,
        "C://Users//caaguila//Documents//fungal_dispersal//data_synthesis//Global_Fungi_long_format.RDS")

