# Data Processing Documentation

## Overview
Brief description of the data processing objectives and scope. Reminder to adhere to data ownership and usage guidelines.

## ACCESS proposal
Dispersal is a key mechanism driving the geographic distributions of species on Earth, but dispersal theory and methods are based primarily on macroorganisms with microbial dispersal paradigms emerging only recently. In fungi, numerous traits related to dispersal (e.g. spore traits, fruiting body traits, dispersal syndromes) are likely linked to fungal biogeographic patterns, but these hypotheses remain largely untested. We aim to harmonize fungal dispersal trait data with DNA sequence-based taxon occurrence data to test trait-based predictions regarding the dispersal capabilities of fungi across spatial scales. We will also assess the potential for fungal dispersal to buffer against range shifts predicted with global climate change. This work will contribute to our understanding of global fungal biodiversity and ecosystem function, as well as aid in predicting plant and human fungal disease outbreaks. Finally, we will integrate fungal dispersal models with global climate change predictions to assess the potential for fungal range shifts in a changing world. 

To accomplish these tasks, most data wrangling and analyses will be performed using the language R. Much of the eDNA data will be in raw sequence form, which will require advanced computing capabilities to enable an efficient bioinformatics pipeline; merging trait-based spore dispersal modeling with climate envelope models will also require substantial computing power. Due to the novelty of this work we are unsure about our specific requirements, but we are excited to learn more about our GPU and storage requirements through benchmarking.

## Data Sources
List and describe data sources used, including links to cloud-optimized sources. Highlight permissions and compliance with data ownership guidelines.

## CyVerse Discovery Environment
Instructions for setting up and using the CyVerse Discovery Environment for data processing. Tips for cloud-based data access and processing.

## Data Processing Steps

### Using GDAL VSI
Guidance on using GDAL VSI (Virtual System Interface) for data access and processing. Example commands or scripts:
```bash
gdal_translate /vsicurl/http://example.com/data.tif output.tif
```

## Cloud-Optimized Data
Advantages of using cloud-optimized data formats and processing data without downloading. Instructions for such processes.

## Data Storage

Information on storing processed data, with guidelines for choosing between the repository and CyVerse Data Store.

## Best Practices

Recommendations for efficient and responsible data processing in the cloud. Tips to ensure data integrity and reproducibility.

## Challenges and Troubleshooting

Common challenges in data processing and potential solutions. Resources for troubleshooting in the CyVerse Discovery Environment.

## Conclusions

Summary of the data processing phase and its outcomes. Reflect on the methods used.

## References

Citations of tools, data sources, and other references used in the data processing phase.
