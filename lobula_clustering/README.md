# Lobula clustering
The scripts in this folder performs the clustering analysis of fragmented terminals in lobula presynaptic to LPLC1 based on their connectivity as well as morphology, to categorize them into putative cell types.

To replicate the results in Tanaka & Clark (2021) bioRxiv, run scripts in the following order:

1. **GetSmallInputsConnectivity.py** This script identifies small lobula intrinsic neurons without cell type labels which are presynaptic to a specified cell type (e.g. LPLC1), and saves their connectivity to labeled cell types as a csv.

2. **GetAndSaveSynapse.py** This script will download and save the XYZ coordinates of all the presynapses of LT1 neurons. LT1 neurons are wide-field monostratified neuron innervating the layer 2 of lobula, which we will use as a landmark relative to which we judge layer affiliation of small terminal synapses.

3. **CalcMorphoStats_SDandLayers.py** This script performs PCA on the LT1 synapse data downloaded in **GetAndSaveSynapse.py**, and then fit a quadratic shell model on them. We then rotate the postsynapses of each fragmented terminal we identified in **GetSmallInputsConnectivity.py** into this LT1-based PC space, and compute the depth of their synapses relative to the LT1 shell model. We then calculate the spatial spread of their synapse locations as their standard deviation along the three PC axis, and also calculate the histogram of their synapse depth. We append these two morphological metrics (spread and depth) to the previously saved connectivity matrix.

4. **ClusterSmallInputsByConnectivityAndMorphology.py** This script performs hierarchical clustering on the fragmented terminals based on their connectivity and morphology.

The results folder contains a pre-computed results after step 3.

### Quantification of T2/T3 connectivity on LPLC1
The script **totalInputFractionT2T3.py** read the annotated list of putative T2/T3 cells (putativeT2T3_annotated.csv) and count how much of the LPLC1 inputs they account for.
