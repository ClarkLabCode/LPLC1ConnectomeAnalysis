# List input and output neuron types of a specified celltype
# Basically the same idea as "subburst" visualization on neuPrint website,
# but combined for all cells in a given celltype

## load packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

## Define parameters
respath = './downstream_analysis/results/' # save results here
cellType = 'LPLC1' # change this to your favourite cell
# cellType = 'PLP219'
# cellType = 'PVLP112'
# cellType = 'PVLP113'
# cellType = 'DNp03'
# cellType = 'DNp06'

## Connect to the server
# NOTE: For this script to run, you must put "authtoken" file under the top directory
f = open("authtoken","r")
tokenstr = f.read()
c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token=tokenstr)
c.fetch_version()

## Part 1: Get all instances of the cell type
# prepare query
q = """\
    MATCH (a:Neuron)
    WHERE a.type='%s'
    RETURN a.bodyId as bodyId
    """ % cellType
cellDF = c.fetch_custom(q)
nCells = len(cellDF)
print('found ',nCells,cellType)

## Part 2: Go through all instances one by one and get in/out
for ii in range(nCells):
    # Get outputs
    q = """\
         MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron)
         WHERE a.bodyId = %s
         RETURN DISTINCT b.bodyId AS bodyId, b.type AS type, w.weight as w
        """ % cellDF.at[ii,"bodyId"]
    dfout = c.fetch_custom(q)
    # Get inputs
    q = """\
         MATCH (b:Neuron)-[w:ConnectsTo]->(a:Neuron)
         WHERE a.bodyId = %s
         RETURN DISTINCT b.bodyId AS bodyId, b.type AS type, w.weight as w
        """ % cellDF.at[ii,"bodyId"]
    dfin = c.fetch_custom(q)

    if ii==0:
        alldfout = dfout
        alldfin  = dfin
    else:
        alldfout = pd.concat([alldfout,dfout])
        alldfin  = pd.concat([alldfin,dfin])

## Part 3 Sort results by cell type
# Separate ones that have and don't have assigned cell types
notypedfout   = alldfout.loc[pd.isnull(alldfout['type'])]
notypedfin    = alldfin.loc[pd.isnull(alldfin['type'])]
withtypedfout = alldfout.loc[~pd.isnull(alldfout['type'])]
withtypedfin  = alldfin.loc[~pd.isnull(alldfin['type'])]

# prepare output structure
resdfout = pd.DataFrame(columns=['Id', 'Wtotal', 'Wpercell']) # here id is either bodyId or type
resdfin  = pd.DataFrame(columns=['Id', 'Wtotal', 'Wpercell'])

# For ones that don't have cell types, sort by bodyId
# Out
uniqueIdsNotypeOut = notypedfout.bodyId.unique()
for thisId in uniqueIdsNotypeOut:
    thisAllW = notypedfout.loc[notypedfout['bodyId']==thisId].w.sum()
    resdfout = resdfout.append({'Id': thisId, 'Wtotal': thisAllW, 'Wpercell': thisAllW/nCells}, ignore_index=True)
# In
uniqueIdsNotypeIn = notypedfin.bodyId.unique()
for thisId in uniqueIdsNotypeIn:
    thisAllW = notypedfin.loc[notypedfin['bodyId']==thisId].w.sum()
    resdfin = resdfin.append({'Id': thisId, 'Wtotal': thisAllW, 'Wpercell': thisAllW/nCells}, ignore_index=True)

# For ones that have cell types, sort by type
# Out
uniqueTypeOut = withtypedfout.type.unique()
for thisType in uniqueTypeOut:
    thisAllW = withtypedfout.loc[withtypedfout['type']==thisType].w.sum()
    resdfout = resdfout.append({'Id': thisType, 'Wtotal': thisAllW, 'Wpercell': thisAllW/nCells}, ignore_index=True)
# In
uniqueTypeIn = withtypedfin.type.unique()
for thisType in uniqueTypeIn:
    thisAllW = withtypedfin.loc[withtypedfin['type']==thisType].w.sum()
    resdfin = resdfin.append({'Id': thisType, 'Wtotal': thisAllW, 'Wpercell': thisAllW/nCells}, ignore_index=True)

resdfout.sort_values(by=['Wtotal'], ascending=False).to_csv(respath+cellType+'_out.csv')
resdfin.sort_values(by=['Wtotal'], ascending=False).to_csv(respath+cellType+'_in.csv')

# visualize
fig, ax = plt.subplots(2,1)

out_con = np.sort(resdfout['Wtotal'])[::-1]
in_con  = np.sort(resdfin['Wtotal'])[::-1]


# ignoring everything below 10
ax[0].plot(np.log10(out_con[out_con>=10]),'g.-')
ax[1].plot(np.log10(in_con[in_con>=10]),'g.-')
ax[0].set_title('log10 out synapses')
ax[1].set_title('log10 in synapses')
plt.show()
