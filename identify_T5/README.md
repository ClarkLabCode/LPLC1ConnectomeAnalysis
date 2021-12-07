# T5 identification
The scripts in this folder identify candidate T5 neurons from the hemibrain dataset based on their synapse locations and their connectivity to the monostratified LPTCs (or lack thereof). The visually annotated lists of putative T5 cells by subtype are provided in the **results** folder as a part of the repository, which is read by **identifyT5c.py** as well as **examineLPLCT5connectivity.py**.

## identifyT5abd.py
This script identifies candidate T5 neurons as cells that are intrinsic to lobula and lobula plate, and have more inputs in lobula and more outputs in lobula plate. Among these putative T5 cells, we select T5a as ones connected to HS and CH, T5b as ones connected to H2, and T5d as ones connected to VS (See for example [Boergens et al. 2018](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0207828) for the anatomy of these cells). We performed visual annotation on the results of this script, which is provided under the **results** folder.

## identifyT5c.py
This script identify T5c neurons as a lobula/lobula plate intrinsic neuron that lacks connectivity to HS, CH, H2, or VS.

## examineLPTCT5connectivity.py
This script examines connectivity from the identified T5 neurons to a specified downstream cell type (e.g. LPLC1) by subtypes.
