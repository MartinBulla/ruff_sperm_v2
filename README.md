## Supporting information for "Sperm traits of the three genetic morphs in the ruff sandpiper "

![Satellite](./Illustrations/morphs.png)

### **Overview**

Contains data collection and preregistration protocols, recordings and pictures of sperm to estimate motility and measure morphology of sperm in three ruff morphs, meta-data on individual males, and computer code to generate outputs of the analysis. 

### **Folders and files**
[Supplement](Supplement.pdf): Supplementary methods, figures and tables

[Data](Data/): raw data and manipulated data generated with R-scripts starting on "DAT_"
- [sperm_recs](Data/sperm_recs/): Recordings of sperm motility
- [sperm_pics](Data/sperm_pics/): Pictures of individual sperm cells ([original](/Data/sperm_pics/original/)), their inverted version ([original](/Data/sperm_pics/inverted/)) of which 920 were chosed for measurements in SpermSizer ([original](/Data/sperm_pics/measured/))
- [sperm_morpho](Data/sperm_morpho/): Folders (named as date and sperm batch) containing 
    - *measure* folder with inverted randomised pictures measured in a given batch; if further manipulation of the picture was needed, an alternative picture is included
    - *Results* folder created by Sperm Sizer software, having date and time of the measurement in the name, and containing pictures of each measured part with indicated measurement lines and line lengths (in pixels), and respective excel files with meaasurements, all corresponding to the picture names in *measure*

- [sim](Data/sim/): mcmc chains

[R](R/)-scripts used in the analysis; those starting with 
- "DAT" prepare data for analysis
- "Fig or Table" create figures or tables
- "MET_sample-sizes.R" provides sample sizes for the Methods and Results section
- "EXP_relatedness_effects.R" explores mixing of the mcmcm chains in the analyses that control for relatedness
- "tools.R" tools necessary for most of the R-scripts and sources within those
- "Out_velocity-videos.Rmd" script to create html docuement showing the ruff and zebra finch motility
- "Preregistration_v3.Rmd" script to create html document of a priory methods of data collection and analyses
- "protocol_sperm.Rmd" script to create html document of sperm sampling, sample preparation, and picture taking protocol
- "_runRmarkdown.R generates htmls from Rmd scripts

[Outputs](Outputs/): all outputs used in the manuscript

[Protocols](Protocols/): Protocols used in this study, including Rmd source code and example videos for respective html documents

[Illustrations](Illustrations): ruff morph pictures used in the graphs, drawn by Yifan Pei, and two example videos of ruff and zebrafinch sperm motility (cut from the original videos from [sperm_recs](Data/sperm_recs/))

[LICENSE](LICENSE): terms of reuse - applicable only after this work is published as a preprint or in a scientific journal, until then the data are not available for reuse.