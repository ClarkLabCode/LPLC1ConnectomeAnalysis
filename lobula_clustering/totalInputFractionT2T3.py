# Estimate how much of total inputs into LPLC1 T2 and T3 respectively account for

## load packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

## Connect to the server
# NOTE: For this script to run, you must put "authtoken" file under the top directory
f = open("authtoken","r")
tokenstr = f.read()
c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token=tokenstr)
c.fetch_version()

# Get the list of manually annotated T2 and T3s
T2T3file = './lobula_clustering/results/putativeT2T3_annotated.csv'
T2T3df = pd.read_csv(T2T3file)
T2list = T2T3df.loc[T2T3df['Putative cell type']=='T2'].bodyId.to_list();
T3list = T2T3df.loc[T2T3df['Putative cell type']=='T3'].bodyId.to_list();
otherlist = T2T3df.loc[T2T3df['Putative cell type']=='other'].bodyId.to_list();

# Get the list of all "small terminals" we analyzed
alllobulafile = './lobula_clustering/results/lobulaTerminals_LPLC1_totSynBelow300_minCon3_2WithMorphoStats.csv'
lobuladf = pd.read_csv(alllobulafile)
print('Total small lobula terminal synapses of LPLC1: ',np.sum(lobuladf['LPLC1']))
SmallTerminalCon = np.sum(lobuladf['LPLC1'])

### Get total synapse count from T2 and T3 to all LPLC1
### T2
q = """\
    MATCH (b:Neuron)-[w:ConnectsTo]->(a:Neuron)
    WHERE a.type='LPLC1' AND b.bodyId IN %s
    RETURN DISTINCT b.bodyId as bodyId, w.weight as weight
    """ % T2list
df = c.fetch_custom(q)
print('Annotated T2: ',len(T2list))
print('Total T2 connection: ',np.sum(df['weight']))
T2con = np.sum(df['weight'])

### T3
q = """\
    MATCH (b:Neuron)-[w:ConnectsTo]->(a:Neuron)
    WHERE a.type='LPLC1' AND b.bodyId IN %s
    RETURN DISTINCT b.bodyId as bodyId, w.weight as weight
    """ % T3list
df = c.fetch_custom(q)
print('Annotated T3: ',len(T3list))
print('Total T3 connection: ',np.sum(df['weight']))
T3con = np.sum(df['weight'])

### other
q = """\
    MATCH (b:Neuron)-[w:ConnectsTo]->(a:Neuron)
    WHERE a.type='LPLC1' AND b.bodyId IN %s
    RETURN DISTINCT b.bodyId as bodyId, w.weight as weight
    """ % otherlist
df = c.fetch_custom(q)
print('Annotated other: ',len(otherlist))
print('Total other connection: ',np.sum(df['weight']))


### Get total lobula postsynapses of LPLC1
q = """\
    MATCH (a:Neuron)
    WHERE a.type='LPLC1'
    RETURN DISTINCT a.bodyId as bodyId, apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].post as lopost
    """
df = c.fetch_custom(q)
print('Total lobula post synapses of LPLC1: ',np.sum(df['lopost']))
LobulaCon = np.sum(df['lopost'])

### Get total postsynapses of LPLC1
q = """\
    MATCH (a:Neuron)
    WHERE a.type='LPLC1'
    RETURN DISTINCT a.bodyId as bodyId, a.post as allpost
    """
df = c.fetch_custom(q)
print('Total post synapses of LPLC1: ',np.sum(df['allpost']))
AllCon = np.sum(df['allpost'])

### Concatenate all the synapse count in a MECE list for a pie chart
conlist = [SmallTerminalCon-T2con-T3con, T2con, T3con]
fig, ax = plt.subplots()
ax.pie(conlist,labels=('Small','T2','T3'))
plt.show()
