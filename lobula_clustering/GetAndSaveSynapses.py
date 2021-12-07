# This script gets and saves synapses locations of a specified cell type

## load packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

## Define parameters
cellType = 'LT1' # change this to your favourite cell
savePath = './lobula_clustering/results/synapse_locations/'+cellType

# make the directory (if this is the first time running this)
try:
    os.makedirs(savePath)
except OSError as error:
    print(error)

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
df = c.fetch_custom(q)
nCells = len(df)
print('found ',nCells,cellType,'(flies have 750 ommatidia on one eye)')

## Part 2: Get synapses
meanLoc = [] # mean synapse loation
allLocs = {} # store all synapse locations by cell
for ii in range(nCells): # we do this one cell at a time because otherwise we get server timeout
    thisId = df.at[ii,'bodyId'] # id of the cell we look at here
    # prepare query
    q = """\
        MATCH (a:Neuron)-[:Contains]->(:SynapseSet)-[:Contains]->(s:Synapse)
        WHERE a.bodyId=%s AND s.`LO(R)` AND s.type='post'
        RETURN DISTINCT s.location.x as x, s.location.y as y, s.location.z as z
        """ % thisId
    df2 = c.fetch_custom(q) # get the result
    # print(df2.head()) # this was for debug
    if not df2.empty: # I need this to avoid zero division error
        if ii==0:
            alldf = df2
        else:
            alldf = pd.concat([df2,alldf])

alldf.to_csv(savePath+'/allPostSynapses.csv')

## Part 3: Plot Synapse
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.scatter(alldf['x'],alldf['y'],alldf['z']) # mean synapse location
ax1.set_xlabel('x (in 8nm px)')
ax1.set_ylabel('y (in 8nm px)')
ax1.set_zlabel('z (in 8nm px)')
plt.show()
