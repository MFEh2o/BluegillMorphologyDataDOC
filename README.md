# Bluegill Shape Data
Relating lake DOC to Geometric Morphometric Analyses of bluegill captured from study lakes.

### Abstract
This repo contains code and data related to Chelsea Bishop's project on Bluegill morphometrics. The data, including a version of the MFE database (where the original landmark data resides) can be found in the data/ folder. Code is in the code/ folder. Archived scripts, which didn't end up being necessary to the complete minimal reproducible pipeline (e.g. earlier drafts of Chelsea's analysis) are thrown into the archived/ folder. The [renv](https://rstudio.github.io/renv/articles/renv.html)/ folder preserves package versions used in this analysis, including older versions of packages that may not be current on CRAN. The resources/ folder has various relevant literature documents and manuscript drafts. The figures/ folder contains some drafts

Here is a diagram illustrating the reproducible data and code pipeline.

![2ChelseaDataSchematic (1)](https://user-images.githubusercontent.com/37053323/113429309-4ec01f00-93a6-11eb-99d3-bc1b5de15085.png)

### Directory structure
- archived/ contains older scripts and data files that aren't part of the final reproducible workflow. Experimentation etc.
- code/ contains final R code files, represented in the schema above.
- data/ contains final data files, including minimal necessary inputs and outputs. 
     - unclassified/ contains a few data files carried over from the archived data folder, to use as comparison for the recreated data files. Where there are discrepancies between these files and the final files in outputs/, they are noted in the recreate*.R scripts.
- figures/ contains saved images that are either final or near-final versions of the figures needed for the manuscript. These are generated in ReviewApril2020_KGedit.R. When further editing is needed, these can be loaded into Photoshop or similar and modified. Because some figures contain multiple parts and I wanted to leave room for arranging those parts as needed, I've created a separate folder for each figure, containing one or multiple image files.
- renv/ contains data about package versions used in preparing this project, to facilitate re-running the analyses at a later date even if packages have been updated.
- resources/ contains various manuscripts and literature that was helpful in preparing these analyses.


```
.
├── BluegillMorphologyDataDOC.Rproj
├── README.md
├── archived
│   ├── code
│   └── data
├── code
│   ├── DOC_binning.R
│   ├── ReviewApril2020_KGedit.R
│   ├── createIdentifiersUpdated.R
│   ├── dbUtil.R
│   ├── hucIDforLakes.R
│   ├── map.R
│   ├── recreateEyewidths.R
│   ├── recreateGillRakers.R
│   ├── recreatePecFinAngles.R
│   └── recreatePecFins.R
├── data
│   ├── MFEdb_20210305.db
│   ├── inputs
│   ├── outputs
│   └── unclassified
├── figures
│   ├── eyeWidths
│   ├── fishShapes_pc1_pc2
│   ├── gillRakers
│   └── pecFins
├── renv
│   ├── activate.R
│   ├── library
│   └── settings.dcf
├── renv.lock
└── resources
    ├── BISHOP_THESIS2020_FINALSUBMISSION.pdf
    ├── Craigetal2017_Ecology&Evolution.pdf
    ├── Kaeuffer_et_al_2012.pdf
    ├── LandmarkDiagrams.pdf
    └── Manuscript Draft_FiguresInText.docx
```

### Contributions
Repository created by Chelsea Bishop and Kaija Gahm in 2020-21. 

Data were collected by Chelsea Bishop with the help of Alex Ross, Colin Dassow, Matthew Farragher, Henry Chung and undergrads of the University of Notre Dame in 2018. Data and scripts were reorganized and reviewed by Kaija Gahm in consultation with Chelsea Bishop, and Kaija prepared the data reproducibility pipeline.

Geometric Morphometrics Analyses were conducted with the help of Madlen Stange. 
Fish morphometric data are contributed to the [MFE database](https://figshare.com/articles/MFE_database_Data_from_ecosystem_ecology_research_by_Jones_Solomon_and_collaborators_on_the_ecology_and_biogeochemistry_of_lakes_and_lake_organisms_in_the_Upper_Midwest_USA/7438598). 

