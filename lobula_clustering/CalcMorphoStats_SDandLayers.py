# Based on the downloaded LT1 synapses, this script appends morphological metrics
# to the previously downloaded connectivity matrix of fragmented terminals
# this can take a while to run

## load packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA

## Define parameters
savePath = './lobula_clustering/results/'
csvName = 'lobulaTerminals_LPLC1_totSynBelow300_minCon3_2' # the name of the output of GetSmallInputsConnectivity
fitSynapsePath = './lobula_clustering/results/synapse_locations/LT1/allPostSynapses.csv'
synapseType = 'post'

## Connect to the server
# NOTE: For this script to run, you must put "authtoken" file under the top directory
f = open("authtoken","r")
tokenstr = f.read()
c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token=tokenstr)
c.fetch_version()

## Part 1: prepare LT1 synapse based shell model of lobula layers
df = pd.read_csv(fitSynapsePath)

## transform data into numpy format
X = df["x"].to_numpy().flatten()
Y = df["y"].to_numpy().flatten()
Z = df["z"].to_numpy().flatten()
XYZ = np.array([X,Y,Z]).T

# Run PCA
pca = PCA(n_components=3)
PCs = pca.fit_transform(XYZ)

# flatten data
PC1 = PCs[:,0].flatten()
PC2 = PCs[:,1].flatten()
PC3 = PCs[:,2].flatten()

# define and fit quadric model to PCs
A = np.array([PC1*0+1, PC1, PC2, PC1**2, PC2**2, PC1*PC2]).T
coeff, r, rank, s = np.linalg.lstsq(A,PC3)

# calculate R squared
rsq = 1 - r/np.sum((PC3-np.mean(PC3))**2)
print('R2 of the quadratic fit: ',rsq)

## Visualize the quadratic fit
meshX, meshY = np.meshgrid(np.arange(-90,90,5)*1000/8, np.arange(-40,80,5)*1000/8) # in 8nm pixels
meshX = meshX.flatten()
meshY = meshY.flatten()
B = np.array([meshX*0+1, meshX, meshY, meshX**2, meshY**2, meshX*meshY]).T
meshZ = np.dot(B,coeff)

# only show a fraction of synapses
inds = np.random.uniform(size=len(X)) < 0.1
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_trisurf(meshX*8/1000,meshY*8/1000,meshZ*8/1000,
                cmap = 'cool', alpha=0.5)
ax.scatter(PC1[inds]*8/1000,PC2[inds]*8/1000,PC3[inds]*8/1000,
           c='black', s=1) # mean synapse location
ax.set_xlabel('PC1 (micron)');
ax.set_ylabel('PC2 (micron)');
ax.set_zlabel('PC3 (micron)');
plt.show()


## Part 2: load synapse info of each small input, and calculate summary stats
# Read csv of small inputs
df = pd.read_csv(savePath+csvName+'.csv')
nCells = len(df)

# prepare columns for morpho stats
df['std1'] = np.zeros(nCells)
df['std2'] = np.zeros(nCells)
df['std3'] = np.zeros(nCells)

# bins for histogram
bins = [-10,-5,0,5,10,15,20,25,30,35,40,45]
for ii in range(len(bins)-1):
    df['L'+str(ii)] =  np.zeros(nCells)

# Go through all the cells
for ii in range(nCells):
    thisId = df.at[ii,'bodyId']
    if ii%20==1:
        print('Processing cell #',ii)
    # Get synapses
    q = """\
        MATCH (a:Neuron)-[:Contains]->(:SynapseSet)-[:Contains]->(s:Synapse)
        WHERE a.bodyId=%s AND s.`LO(R)` AND s.type='%s'
        RETURN DISTINCT s.location.x as x, s.location.y as y, s.location.z as z
        """ % (thisId,synapseType)
    sdf = c.fetch_custom(q)

    # transform synapse location data into a format compatible to the model
    X = sdf["x"].to_numpy().flatten()
    Y = sdf["y"].to_numpy().flatten()
    Z = sdf["z"].to_numpy().flatten()
    XYZ = np.array([X,Y,Z]).T

    # Rotate synapse
    PCs = pca.transform(XYZ)
    PC1 = PCs[:,0].flatten()
    PC2 = PCs[:,1].flatten()
    PC3 = PCs[:,2].flatten()

    # calculate expected synapse z depth from the shell model
    A = np.array([PC1*0+1, PC1, PC2, PC1**2, PC2**2, PC1*PC2]).T
    depthError = PC3 - np.dot(A,coeff) # positive = deeper
    depthError = depthError*8/1000 # in micron

    # calculate diffuseness in the three PC directions
    df.at[ii, 'std1'] = np.std(PC1*8/1000)
    df.at[ii, 'std2'] = np.std(PC2*8/1000)
    df.at[ii, 'std3'] = np.std(PC3*8/1000)

    # sort synapses by depth at 5 micron bin
    synHist = np.histogram(depthError,bins=bins)
    print(synHist)
    for jj in range(len(bins)-1):
        df.at[ii, 'L'+str(jj)] = synHist[0][jj]

# save results
df.to_csv(savePath+csvName+'WithMorphoStatsLayers.csv')
