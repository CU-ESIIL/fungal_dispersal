# Global Spore Sampling Project Data Prep for ESILL Group Analysis

# Set working directory and load required packages
setwd("C:/Users/egan_/Documents/R work/ESIIL_Fungal_Dispersal/Global_Spore_Sampling_Project/")
install.packages("kgc")
install.packages("seqinr")
install.packages("tidyverse")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("maps")
library(ggplot2)
library(maps)
library(reshape2)
library(kgc)
library(seqinr)
library(tidyverse)


# Load open data available from publication (https://zenodo.org/api/records/11125610/files-archive)
load("C:/Users/egan_/Documents/R work/ESIIL_Fungal_Dispersal/Global_Spore_Sampling_Project/sample_data.RData")
load("C:/Users/egan_/Documents/R work/ESIIL_Fungal_Dispersal/Global_Spore_Sampling_Project/site_data.RData")


# Pre-process data using S1_generate_data_tables R script available from download
datafolder="C:/Users/egan_/Documents/R work/ESIIL_Fungal_Dispersal/Global_Spore_Sampling_Project/bioinformatics_pipeline_outputs/bioinformatics pipeline outputs"
taxonomy = readRDS(file=paste0(datafolder,"/otu_taxonomy.rds"))
for(j in 4:10) taxonomy[,j] = as.factor(taxonomy[,j])

otu.table = t(readRDS(file=paste0(datafolder,"/otu_table.rds")))
otu.sequences = read.fasta(file=paste0(datafolder,"/otu.gz"),as.string = T,forceDNAtolower=F)

ma1 = match(rownames(otu.table),meta$sample.id)
ma1 = ma1[!is.na(ma1)]
meta = droplevels(meta[ma1,])

ma2 = match(meta$sample.id,rownames(otu.table))
otu.table = otu.table[ma2,]

b=read.table(paste0(datafolder,"/read_counts.tsv"),header=T)
ma3 = match(meta$sample.id,b$sample)
b=b[ma3,]

meta$numnonspikes = b$fungi_nread
meta$numspikes = b$nochim2_nread-b$nospike_nread
w = ((meta$numnonspikes+1)/(meta$numspikes+1))/meta$duration

meta$dna_amount = w * meta$spike_dilution*45/2/0.78 #See Ovaskainen et al. FEE (2020)
meta$dna_amount = log(meta$dna_amount)

lat = rep(NA,nrow(meta))
lon = rep(NA,nrow(meta))
temp.mean = rep(NA,nrow(meta))
CZ = rep(NA,nrow(meta))
for(i in 1:nrow(site.data)){
  lat[meta$site==site.data$Site.code[i]] = site.data$LAT[i]
  lon[meta$site==site.data$Site.code[i]] = site.data$LONG[i]
  temp.mean[meta$site==site.data$Site.code[i]] = site.data$temp.mean.site[i]
}
meta$lat = lat
meta$lon = lon
meta$temp.mean = temp.mean

my.order = order(meta$site,meta$date)

meta = meta[my.order,]
otu.table = otu.table[my.order,]
rownames(meta) = meta$sample.id

if(!all(rownames(taxonomy)==names(otu.sequences))) print("problem with otu names")
taxonomy$sequence = as.character(otu.sequences)
taxonomy$ref_seq_id = NULL

all(colnames(otu.table)==rownames(taxonomy))

save(meta, otu.table, taxonomy,
     file = "C:/Users/egan_/Documents/R work/ESIIL_Fungal_Dispersal/Global_Spore_Sampling_Project/processed_data/allData.Rdata")


# Begin Processing Data Tables to be in Communal ESILL Format

write.csv2(meta,file = "processed_data/metadata.csv",row.names = FALSE)
write.csv2(taxonomy,file = "processed_data/taxonomy.csv")
write.csv2(otu.table,file = "processed_data/otu.table.csv")
spore_otu <- t(otu.table)

# Combine spore_otu and taxonomy by row names
otu_tax.combined <- merge(spore_otu, taxonomy, by = "row.names", all = TRUE)

# Reshape spore_otu to long format using reshape2::melt() 
otu_tax.long <- melt(otu_tax.combined, varnames = c("Row.names", "kingdom","phylum","class",
                                                    "order","family","genus","species","sequence"),
                     variable.name = "sample.id", 
                     value.name = "Abundance")


# The next step is computationally heavy. Be prepared!
gc()             # Run garbage collection to free memory


# Add in metadata using tidyverse::left_join
global_spore.combined <- left_join(otu_tax.long, meta, by = "sample.id")

# Remove unneeded columns from new dataframe
global_spore.esiil <- global_spore.combined[, !names(global_spore.combined) %in% c("kingdom", "sequence","seqrun","site",
                                                                                      "yday","duration",
                                                                   "water","insect","unst.tweezers","spike_dilution",
                                                                   "numnonspikes","numspikes","dna_amount","temp.mean")]


# Update to have new column names to match other databases and add any columns needed
 "OTU"                 "SampleID"            "Abundance"           "latitude"            "longitude"          
 "sample_type"         "add_date"            "paper_id"            "sequencing_platform" "target_gene"        
 "phylum"              "order"               "family"              "genera"              "species"            
 "Database"     

# Update column names
global_spore.esiil <- global_spore.esiil %>% rename(OTU = Row.names, 
                                                          SampleID = sample.id, latitude = lat,
                                                          longitude = lon, add_date = date)

# Add a new column names and add in values (if possible)
global_spore.esiil$sample_type <- "air"
global_spore.esiil$paper_id <- "Ovaskainen et al. 2024"
global_spore.esiil$sequencing_platform <- "MiSeq"
global_spore.esiil$target_gene <- "ITS2"
global_spore.esiil$Database <- "Global Spore"


# Plot map of Europe
world_map <- map_data("world")

# Count the occurrences of each species and find the most abundant species from the dataframe
species_count <- table(global_spore.esiil$species)
most_abundant_species <- names(species_count[which.max(species_count)])
cat("The most abundant species is:", most_abundant_species)

# Filter the dataframe for the species of interest, e.g., "specific_species_name"
# Replace "specific_species_name" with the actual species name you're interested in
filtered_species <- global_spore.esiil[global_spore.esiil$species == "Abortiporus_biennis_283905", ]

# Plot the distribution map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray") +
  geom_point(data = filtered_species, aes(x = longitude, y = latitude, shape = sample_type, color = sample_type), size = 3) +
  coord_cartesian(xlim = c(-25, 45), ylim = c(35, 70)) +  # Europe boundaries
  theme_minimal() +
  labs(
    title = paste("Geographic Distribution of", "Abortiporus biennis", "in Europe"), 
    x = "Longitude", 
    y = "Latitude"
  )
