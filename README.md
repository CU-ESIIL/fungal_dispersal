[![DOI](https://zenodo.org/badge/780491969.svg)](https://zenodo.org/doi/10.5281/zenodo.11167974)

# Advancing fungal dispersal ecology through traits and data harmonization

Welcome to the **fungal dispersal** repository, an integral part of the Environmental Data Science Innovation and Inclusion Lab (ESIIL). This repository is the central hub for our working group, encompassing our project overview, proposals, team member information, codebase, and more.

## Our Project
Dispersal is a key mechanism driving the geographic distributions of species on Earth, but dispersal theory and methods are based primarily on macroorganisms with microbial dispersal paradigms emerging only recently. In fungi, numerous traits related to dispersal (e.g. spore traits, fruiting body traits, dispersal syndromes) are likely linked to fungal biogeographic patterns, but these hypotheses remain largely untested. We aim to harmonize fungal dispersal trait data with DNA sequence-based taxon occurrence data to test trait-based predictions regarding the dispersal capabilities of fungi across spatial scales. We will also assess the potential for fungal dispersal to buffer against range shifts predicted with global climate change. This work will contribute to our understanding of global fungal biodiversity and ecosystem function, as well as aid in predicting plant and human fungal disease outbreaks. Finally, we will integrate fungal dispersal models with global climate change predictions to assess the potential for fungal range shifts in a changing world.

## Documentation
- Access detailed documentation on our [GitHub Pages site](https://your-gh-pages-url/).
- Find comprehensive guides, tutorials, and additional resources.

## Project Proposal
[Fungal Dispersal Working Group Proposal_11-1-23.pdf](https://github.com/CU-ESIIL/fungal_dispersal/files/14826675/Fungal.Dispersal.Working.Group.Proposal_11-1-23.pdf)

## Group Members
- PI: Dr. Bala Chaudhary, Dartmouth College.
- Co-PI: Dr. Cameron Egan, University of Southern California.
- Co-PI: Dr. Kabir Peay, Stanford University.
- Co-technical lead: Dr. Carlos A. Aguilar-Trigueros, University of Jyväskylä.
- Co-technical lead: Sarah Cuprewich, Dartmouth College.
- Dr. Michelle Afkhami, University of Miami.
- Kristin Barbour, UC Irvine.
- Dr. Priscila Chaveri, Bowie State University.
- Dr. Veera Norros, Finnish Environment Institute.
- Dr. Anne Pringle, University of Wisconsin Madison.
- Dr. Adriana Romero-Olivares, New Mexico State University.
- Dr. Agnese Seminara, University of Genoa.
- Dr. Ryan Stephens, Eastern Tennessee State University.
- Lauren Ward, Stanford University.
- Dr. Kira Lynn, University of Pretoria.
- Dr. Robert Ramos, ESIIL.

## Data assemble

We targeted two main sources of data: 

1. Occurrence of fungi based on DNA based approaches (metabarcoding)

we compiled data from three databases: global fungi, GSMC and GSSP. Based on this data we found very skewed sampling intensity at global scale
![Alt text](Figures/sampling_intensity.png)

Looking at each data location from each dataset:
![Alt text](Figures/GF_data.png)

![Alt text](Figures/GSMC_data.png)

![Alt text](Figures/GSSP_data.png)


2. Dispersal traits across described species of fungi

The biggest dataset correspond to spore traits. Combining species occureences with spore size data, we obtained roughly 20-30% coverage

![Alt text](Figures/spore_data_available.png)


Correlating spores size to absolute spore range (measured as the maximum distance at which species are recorded), we found that large spore species in the phylum Ascomycota tend to have larger ranges (tests will follow)

![Alt text](Figures/distances_spore_volume.png)


## Repository Structure
- **Analysis Code**: Scripts for data analysis, statistical modeling, etc.
- **Data Processing**: Scripts for cleaning, merging, and managing datasets.
- **Visualization**: Code for creating figures, charts, and interactive visualizations.

## Meeting Notes and Agendas
- Regular updates to keep all group members informed and engaged with the project's progress and direction.

## Contributing to This Repository
- Contributions from all group members are welcome.
- Please adhere to these guidelines:
  - Ensure commits have clear and concise messages.
  - Document major changes in the meeting notes.
  - Review and merge changes through pull requests for oversight.

## Getting Help
- If you encounter any issues or have questions, please refer to the [ESIIL Support Page](https://esiil-support-page-url/) or contact the repository maintainers directly. https://cu-esiil.github.io/Postdoc_OASIS/resources/cyverse_hacks/

## Customize Your Repository
- **Edit This Readme**: Update with information specific to your project.
- **Update Group Member Bios**: Add detailed information about each group member's expertise and role.
- **Organize Your Code**: Use logical structure and clear naming conventions.
- **Document Your Data**: Include a data directory with README files for datasets.
- **Outline Your Methods**: Create a METHODS.md file for methodologies and tools.
- **Set Up Project Management**: Use 'Issues' and 'Projects' for task tracking.
- **Add a License**: Include an appropriate open-source license.
- **Create Contribution Guidelines**: Establish a CONTRIBUTING.md file.
- **Review and Merge Workflow**: Document your process for reviewing and merging changes.
- **Establish Communication Channels**: Set up channels like Slack or Discord for discussions.

Remember, the goal is to make your repository clear, accessible, and useful for all current and future members of your working group. Happy researching!
