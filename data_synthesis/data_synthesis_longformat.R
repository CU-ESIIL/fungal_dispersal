### Here is to merge the datasets in the long format#####

rm(list = ls()); gc()

# Global fungi (as provided by Petr Baldrian on September 2024)

Global_fungi <- readRDS("C://Users//caaguila//Dropbox//ESIL_data//fungal_dispersal//data_synthesis//Global_Fungi_long_format.RDS")
Global_fungi$species <- gsub("s__", "", Global_fungi$species)

Global_fungi_s <- Global_fungi[c(1:1000),]
Global_fungi_s$species <- gsub("s__", "", Global_fungi_s$species)


spore_data <- readRDS("data_synthesis//carlos_spore_data_ver12Nov21.RDS")
spore_data$names_to_use <- gsub(" ", "_", spore_data$names_to_use)

library(gggibbous)

global_species <- unique(Global_fungi$species)
spore_species <- unique(spore_data$names_to_use)

any(unique(spore_data$names_to_use)%in%unique(Global_fungi$species))

length(which(global_species%in%spore_species))/length(global_species)

ggplot(data.frame(x = 1, y = 0, ratio = c(0.25, 0.75), right = c(TRUE, FALSE)),
       aes(x = x, y = y)) +
  geom_moon(aes(ratio = ratio, fill = right, right = right),
            size = 40,  key_glyph = draw_key_full_moon) +
  theme_void()

simple_global <- distinct(Global_fungi[,c("phylum", "order", "family","genera", "species")])

table(simple_global$phylum[which(simple_global$species%in%spore_species)])
table(simple_global$phylum)

asco_percentage <- 3173/11900
basidio_percentage <- 2197/9247



load('miseq/output/workspace/geoextent_brm_host.RData')



## geo distance function
geo_distance <-
  function(x){  # at least 3 samples
    x <- x %>% dplyr::select(longitude, latitude) %>% distinct
    if(nrow(x) <= 2) res <- NA
    
    if(nrow(x) > 2) res <- max(distVincentyEllipsoid(x[, c('longitude', 'latitude')]))  # for maximum extent
    # if(!is.na(res)) if(res <= 0.1) res <- NA
    return(res)
  }


fungal_list <- 
  split(Global_fungi, interaction(Global_fungi$species,
                                  Global_fungi$target_gene,
                                  sep='__', drop=TRUE))

ditancias <- 
sapply(fungal_list, geo_distance)

saveRDS(ditancias, "distancias.RDS")

distancias <- readRDS("distancias.RDS")

distancias <- data.frame(species_gene = names(distancias),
                        geo_distances =  distancias, row.names = NULL)

####

geo_area <- function(x){
  x <- x %>% dplyr::select(longitude, latitude) %>% distinct
  if(nrow(x) <= 2) res <- NA
  if(nrow(x) > 2) res <- tryCatch(
    {
      getDynamicAlphaHull(x, coordHeaders=c('longitude','latitude'), 
                          clipToCoast='terrestrial', initialAlpha=3)  # Try to compute the alpha hull
    },
    error = function(e) {
      message("Error encountered: ", e$message)  # Print error message
      return(NULL)  # Return NULL if an error occurs
    }
  ) # for alpha hull
  return(res)

  }

ranges <- sapply(fungal_list, geo_area)
geo.alpha <- ranges

geo.area <- sapply(geo.alpha[sapply(geo.alpha,function(x)class(x[[1]]) == 'SpatialPolygons')],
                   function(x)gArea(spTransform(x[[1]], CRS('+proj=moll +ellps=WGS84'))))



y <- fungal_list$Agaricales_sp__ITS1 %>% select(longitude, latitude) %>% distinct
y <- fungal_list[[1]] %>% select(longitude, latitude) %>% distinct

x <- crotalus[which(crotalus$genSp == 'Crotalus_atrox'),]
x <- x[sample(1:nrow(x), 50),]
x <- x[,c("decimallongitude", "decimallatitude")]

range <- getDynamicAlphaHull(x, coordHeaders=c('decimallongitude','decimallatitude'), 
                             clipToCoast = 'terrestrial')

range_y <- getDynamicAlphaHull(y, coordHeaders=c('longitude','latitude'), 
                             clipToCoast = 'terrestrial', initialAlpha = 3)


geo_distance <-
  function(x){  # at least 3 samples
    x <- x %>% dplyr::select(longitude, latitude) %>% distinct
    if(nrow(x) <= 2) res <- NA
    
    if(nrow(x) > 2) res <- max(distVincentyEllipsoid(x[, c('longitude', 'latitude')]))  # for maximum extent
    # if(!is.na(res)) if(res <= 0.1) res <- NA
    return(res)
  }




# Global Spore Sampling project (as donwloaded from the Scientific Data publication)



# Global Mycologycal Initiative (Leho Tedersoo)