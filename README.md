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
     - This analysis draws from MFEdb_20210423.db, which is version 4.5.2 of the database. This is the same version that's stored on [the MFE Figshare repository](https://caryinstitute.figshare.com/articles/dataset/MFE_database_Data_from_ecosystem_ecology_research_by_Jones_Solomon_and_collaborators_on_the_ecology_and_biogeochemistry_of_lakes_and_lake_organisms_in_the_Upper_Midwest_USA/7438598) as Figshare v4 (https://doi.org/10.25390/caryinstitute.7438598.v4).
     - unclassified/ contains a few data files carried over from the archived data folder, to use as comparison for the recreated data files. Where there are discrepancies between these files and the final files in outputs/, they are noted in the recreate*.R scripts.
- figures/ contains saved images that are either final or near-final versions of the figures needed for the manuscript. These are generated in ReviewApril2020_KGedit.R. When further editing is needed, these can be loaded into Photoshop or similar and modified. Because some figures contain multiple parts and I wanted to leave room for arranging those parts as needed, I've created a separate folder for each figure, containing one or multiple image files.
- renv/ contains data about package versions used in preparing this project, to facilitate re-running the analyses at a later date even if packages have been updated.
- resources/ contains various manuscripts and literature that was helpful in preparing these analyses.


```
.
├── BluegillMorphologyDataDOC.Rproj
├── README.md
├── archived/
│   ├── code/
│   └── data/
├── code/
│   ├── DOC_binning.R
|   ├── calculateDOC.R
|   ├── landmarksToMeasurements.R
│   ├── ReviewApril2020_KGedit.R
│   ├── createIdentifiersUpdated.R
│   ├── dbUtil.R
│   ├── hucIDforLakes.R
│   ├── recreateEyewidths.R
│   ├── recreateGillRakers.R
│   ├── recreatePecFinAngles.R
│   └── recreatePecFins.R
├── data/
│   ├── MFEdb_20210423.db
│   ├── inputs/
│   ├── outputs/
│   └── unclassified/
├── figures/
│   ├── eyeWidths/
│   ├── fishShapes_pc1_pc2/
│   ├── gillRakers/
│   └── pecFins/
├── renv/
│   ├── activate.R
│   ├── library/
│   └── settings.dcf
├── renv.lock
└── resources/
    ├── BISHOP_THESIS2020_FINALSUBMISSION.pdf
    ├── Craigetal2017_Ecology&Evolution.pdf
    ├── Kaeuffer_et_al_2012.pdf
    ├── LandmarkDiagrams.pdf
    └── Manuscript Draft_FiguresInText.docx
```
### Data sources
- **MFEdb_20210402.db** Many of the recreate*.R scripts draw from data in the FISH_MORPHOMETRICS table of the MFE database. The database version that is used in this analysis is version 4.5.4, from 2021-04-23. It can be downloaded from [the MFE Figshare repository](https://caryinstitute.figshare.com/articles/dataset/MFE_database_Data_from_ecosystem_ecology_research_by_Jones_Solomon_and_collaborators_on_the_ecology_and_biogeochemistry_of_lakes_and_lake_organisms_in_the_Upper_Midwest_USA/7438598) as Figshare v4 (https://doi.org/10.25390/caryinstitute.7438598.v4).
- **Bishop_NotoriousBLG.xlsx** manually-created file with raw landmark measurements, mostly copied/exported from the [tpsDIG program](http://www.sbmorphometrics.org/soft-dataacq.html).
- **fishBodyPhotos_fileNames.txt** created by running `ls()` on the folder "fishBodyPhotos", stored in the MFE Box drive under "MFE/notoriousBLGPhotoArchive".
- **FULL_2018_TPS_FILE_UPDATED_09-25-19.TPS** created by tpsDIG
- **Full_body_links.txt** created manually by Chelsea Bishop, based on knowledge of how the body landmarks interact. This file serves as a precursor to geometric morphometric analyses.
- **lakeInfo.csv** mostly pulled from LAKES in the database, but with DOC values originating elsewhere. #XXX need to flesh this out and create a full workflow!
- **PecFins2018.TPS** created by tpsDIG
- **proposedFishScapeLakeList_20180327.csv** # XXX flesh out the DOC workflow here.
- **region4/** and **region7/** shapefiles are from the USGS NHDPlus Version 2 dataset. Source links: [region 4 (Great Lakes)](https://nhdplus.com/NHDPlus/NHDPlusV2_04.php), [region 7 (Upper Mississippi)](https://nhdplus.com/NHDPlus/NHDPlusV2_07.php), [general NHDPlus page](https://nhdplus.com/NHDPlus/NHDPlusV2_data.php)

### Contributions
Repository created by Chelsea Bishop and Kaija Gahm in 2020-21. 

Data were collected by Chelsea Bishop with the help of Alex Ross, Colin Dassow, Matthew Farragher, Henry Chung and undergrads of the University of Notre Dame in 2018. Data and scripts were reorganized and reviewed by Kaija Gahm in consultation with Chelsea Bishop, and Kaija prepared the data reproducibility pipeline.

Geometric Morphometrics Analyses were conducted with the help of Madlen Stange. 
Fish morphometric data are contributed to the [MFE database](https://figshare.com/articles/MFE_database_Data_from_ecosystem_ecology_research_by_Jones_Solomon_and_collaborators_on_the_ecology_and_biogeochemistry_of_lakes_and_lake_organisms_in_the_Upper_Midwest_USA/7438598). 

