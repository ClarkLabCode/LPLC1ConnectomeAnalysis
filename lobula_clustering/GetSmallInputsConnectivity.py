# This scripts identifies small, lobula intrinsic terminals presynaptic to a specified cell types
# and save their id as well as their connectivity in a csv file
# This will take a while to run

## load packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

## Define parameters
cellType = 'LPLC1' # change this to your favourite cell
savePath = './lobula_clustering/results/'

maxSynapses = 300 # Get only small ones
minMainConnection = 3; # We only care about terminals that have more than this much synapses on the target type
minSubConnection = 2;  # We ignore postsynaptic cells that have only synpases less than this much

## Connect to the server
# NOTE: For this script to run, you must put "authtoken" file under the top directory
f = open("authtoken","r")
tokenstr = f.read()
c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token=tokenstr)
c.fetch_version()

## Part 1: Get lobula intrinsic inputs with small size
# prepare query
q = """\
    MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron)
    WHERE b.type='%s' AND apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].pre=a.pre
    AND apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].post=a.post
    AND a.pre+a.post<%s AND w.weight>=%s
    RETURN DISTINCT a.bodyId as bodyId
    """ % (cellType,maxSynapses,minMainConnection)
df = c.fetch_custom(q)
nCells = len(df)
print('found ',nCells,' unique lobula intrinsic neurons connecting to '+cellType+' with less than '+str(maxSynapses)+' synapses, and at least '+str(minMainConnection)+ ' synapses to one '+cellType)


## Part 2
# prepare output structure
outdf = df # keep appending new columns to this

# Go through all the identified cells, get their outputs
for ii in range(len(df)):
    thisId = df.bodyId[ii]
    q = """\
        MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron)
        WHERE a.bodyId=%s AND w.weight>=%s
        RETURN DISTINCT b.bodyId as bodyId, b.type as type, w.weight as w
        """ % (thisId,minSubConnection)
    df2 = c.fetch_custom(q)
    outtypes = df2.type.unique()
    if ii%10==0:
        print('Working on cell #'+str(ii),flush=True)
    for thisType in outtypes:
        if thisType is not None: # ignore ones with no type
            thisCon = np.sum(df2.loc[df2["type"]==thisType].w)
            # if this postsynaptic type has not been already registered
            if thisType not in outdf.columns:
                outdf[thisType] = np.zeros(nCells)
            outdf.at[ii, thisType] = thisCon

outdf.to_csv(savePath+'lobulaTerminals_'+cellType+'_totSynBelow'+str(maxSynapses)+'_minCon'+str(minMainConnection)+'_'+str(minSubConnection)+'.csv');
