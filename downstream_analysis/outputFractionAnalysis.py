# List major postsynaptic neurons of LPLC1 and show
# (1) how much of LPLC1 output they account for
# (2) how much of their input LPLC1 account for
# We use the knowledge obtained by running GenericInOutAnalysisByType

## load packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

## Define parameters
respath = './lplc1/DownstreamAnalysis/'
pretype = 'LPLC1' # change this to your favourite cell
posttypes = ['LPLC1','PLP219','PVLP112','PVLP113','DNp03','DNp06']

## Connect to the server
# NOTE: For this script to run, you must put "authtoken" file under the top directory
f = open("authtoken","r")
tokenstr = f.read()
c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token=tokenstr)
c.fetch_version()

counts = []

## Part 1: Get all instances of the cell type
# prepare query
q = """\
    MATCH (a:Neuron)
    WHERE a.type='%s'
    RETURN DISTINCT a.bodyId as bodyId, a.downstream as pre
    """ % pretype
df = c.fetch_custom(q)
postcount = len(df)
allout = np.sum(df['pre'])
counts.append(allout)

print('found ',postcount,pretype,'(flies have 750 ommatidia on one eye)')
print('total presynapse',allout)

allcon = 0
for cell in posttypes:
    q = """\
        MATCH (b:Neuron)
        WHERE b.type='%s'
        RETURN DISTINCT b.bodyId as bodyId, b.upstream as post
        """ % cell
    celldf = c.fetch_custom(q)
    q = """\
        MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron)
        WHERE a.type='%s' AND b.type='%s'
        RETURN w.weight as con
        """ % (pretype, cell)
    condf = c.fetch_custom(q)
    thisallin = np.sum(celldf['post'])
    thiscon = np.sum(condf['con'])
    allcon += thiscon
    print('found ',thiscon,'synapses from',pretype,'to',cell)
    print('This amounts to',thiscon/allout*100,'% of total',pretype,'output')
    print(cell,'has total',thisallin,'inputs, so synapses from LPLC1 accounts for',thiscon/thisallin*100,'% of it')
    counts.append(thiscon)

print('These',len(posttypes),'cell types have total',allcon,'synapses from LPLC1')
print('Which amounts to',allcon/allout*100,'% of total',pretype,'output')

## visualize
fig, ax = plt.subplots()
piecounts = [counts[0]-np.sum(counts[1:])] # all other output
piecounts = piecounts + counts[1:]
labels = ['Other'] + posttypes
ax.pie(piecounts,labels=labels,autopct='%1.1f%%')
plt.show()
