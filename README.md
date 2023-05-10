# Krogan Lab R Shiny Tools

## Tools Available:
- [Proteomics Summary Analyses](#Proteomics-Summary-Analyses)
- [PPI Data Scoring](#PPI-Data-Scoring)

---

**What can I find here?**  
User-friendly Shiny apps, source code, and documentation for Krogan Lab mass spec proteomics analysis pipelines.  
**How to access these apps?**  
Links are provided on this github page for all finished apps. The apps on this page are hosted on a UCSF server and can only be accessed while logged in to a UCSF internet network. Non UCSF personnel will have to download the source code and run locally in RStudio or host the code on a shiny server ([shinyapps.io](https://www.shinyapps.io/?_gl=1*yakayf*_ga*NzE3MDY1MjQ0LjE2ODMyMzYxMjc.*_ga_8QJS108GF1*MTY4MzMzNjAxNC4yLjEuMTY4MzMzNzU4My4wLjAuMA..*_ga_2C0WZ1JHG0*MTY4MzMzNjAxNC4yLjEuMTY4MzMzNzU4My4wLjAuMA..) is free and easy to use).  
**What is Shiny?**  
[Shiny](https://shiny.rstudio.com/) is an R package for building web apps to make data processing, exploration, and visualization in R accessible by a GUI for non-R-users. Shiny apps can be run locally from the source code in [RStudio](https://posit.co/download/rstudio-desktop/) or hosted on a server and accessed over the web.  
**What's the plan?**  
The general plan for this resource is to add apps for the most routine data preparation and analysis steps first. Creating a Shiny app with adequate flexibility, even from existing analyses in R, can be very time consuming, so we're prioritizing work on the most common and standardized steps. Eventually, the whole analysis pipeline from processing peptide intensities to making figures from specialized analyses like kinase activity prediction, network modularization, and others might be implemented as easy-to-use, modular Shiny apps. The first things on the list are proteomics data preparation with MSstats, proteomics quality control analyses, gene ontology enrichment analysis, and kinase 'activity' enrichment analysis.  


### Proteomics Summary Analyses
Link to app server: http://higgs2.ucsf.edu:3838/bpolacco/shinyKrogan/ProteomicsSummaryAnalyses/ 

This app is for creating volcano plots, bar charts, PCA plots, and heatmaps from proteomics data (abundance or PTMs) that has already been run through MSstats. Labels, colors, thresholds, and common labelling/formatting schemes for the figures are all configurable. Files are downloadable individually as .PDFs.  
- Requires the "GroupComparison" and "ProteinLevelData" .csv files from MSstats (a .gz zipped .csv is also acceptable).

*If you have a working understanding of R, you can generate MSstats results files yourself using the [Proteomics_spec2MSStats_results.rmd](https://github.com/kroganlab/shinyKrogan/blob/main/ProteomicsMSsPreparation/Proteomics_spec2MSStats_results.rmd) file in shinyKrogan/ProteomicsMSsPreparation or using [artMS](https://github.com/biodavidjm/artMS)---or ask someone from the computational corner for help!*

### PPI Data Scoring
Link to app server: http://higgs2.ucsf.edu:3838/yzhou5/R_shiny_for_PPI/PPI_scoring

This app is for generating artMS quality control figures for AP-MS data, and calculating SAINT, MIST, and CompPASS scores. Figures and the updated data file can be downloaded as zipped .PDF and .txt files, respectively. *To remove a run from the analysis, delete the row with it's rawfile and condition names from the keys.txt file.*
- Requires evidence.txt, keys.txt, and the .fasta file from MaxQuant. 
- For examples of properly formatted evidence, keys, and fasta files, see [sample_data](https://github.com/kroganlab/shinyKrogan/blob/main/PPI_scoring/sample_data) in shinyKrogan/PPI_scoring. 
- For more information on how to use the PPI data scoring app, see the [slide presentation](https://github.com/kroganlab/shinyKrogan/blob/main/PPI_scoring/PPI%20scoring.pdf) in shinyKrogan/Documentation.
