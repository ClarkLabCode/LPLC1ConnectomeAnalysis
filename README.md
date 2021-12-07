# LPLC1ConnectomeAnalysis
This repository contains three sets of codes to perform analysis of the hemibrain dataset presented in Tanaka & Clark (2021) bioRxiv as follows:

(1) Identification of putative T5 neurons (**identify_T5**)

(2) Clustering of lobula columnar inputs into LPLC1 (**lobula_clustering**)

(3) Identification of downstream pathways of LPLC1 (**downstream_analysis**)

### How to run the scripts
The scripts here interface with the hemibrain dataset using [neuprint-python](https://github.com/connectome-neuprint/neuprint-python).
To access the dataset, you need to sign up on [neuPrint+](https://neuprint.janelia.org/) (Google account required). Find your Auth Token from the Account menu on neuPrint+, and save it under this directory as a plane text file named **authtoken** (without extension). The scripts read this **authtoken** file to connect to the hemibrain database.
Other than neuprint-python, the scripts require common packages like numpy, scipy, pandas, matplotlib, sklearn etc.
