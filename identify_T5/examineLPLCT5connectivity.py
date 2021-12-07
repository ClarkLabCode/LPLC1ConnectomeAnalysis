# This script reads annotated list of T5s and examine their connectivity to specific
# postsynaptic cell types as the totla number of synapses as well as the fraction of
# T5 cells connected to the given cell type per subtype

## load packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

## Connect to the server
# NOTE: For this script to run, you must put "authtoken" file under the top directory
f = open("authtoken","r")
tokenstr = f.read()
c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token=tokenstr)
c.fetch_version()

## parameters
targetType = 'LPLC1'
T5types = ('a','b','c','d')
totalCounts = []
fractionConnected = []

## Go through T5 subtypes, quantify their connectivity to the specified postsynaptic type
for type in T5types:
    df = pd.read_csv('./identify_T5/results/candidateT5'+type+'_annotated.csv')
    thisIdList = df['bodyId'].to_list()
    q = """\
        MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron)
        WHERE a.bodyId IN %s AND b.type='%s'
        RETURN DISTINCT a.bodyId as preId, b.bodyId as postId, w.weight as weight
        """ % (thisIdList,targetType)
    df = c.fetch_custom(q)
    print('Found',np.sum(df['weight']),'synapses from T5',type,'to',targetType)
    totalCounts.append(np.sum(df['weight']))
    fractionConnected.append(len(np.unique(df['postId']))/len(thisIdList))
    print(len(np.unique(df['postId'])),'out of',len(thisIdList),'T5',type,'were connected to',targetType)
    for postId in list(set(df['postId'].to_list())):
        print(postId, np.sum(df.loc[df['postId']==postId].weight))

## Visualization
fig, ax = plt.subplots(2, 1)
ax[0].bar(range(4), totalCounts)
ax[0].set_xticks(range(4))
ax[0].set_xticklabels(T5types)
ax[0].set_ylabel('total synapse counts')

ax[1].bar(range(4), fractionConnected)
ax[1].set_xticks(range(4))
ax[1].set_xticklabels(T5types)
ax[1].set_ylabel('fraction connected')

plt.show()
