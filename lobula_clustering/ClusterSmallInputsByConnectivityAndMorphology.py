# This scripts runs hiearchical clustering on the saved connectivity and morphology data

## load packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as sch

# Note: agglomerateive clustering is deterministics so you can run this without worrying about noise

## Define parameters
dataPath = './lobula_clustering/results/'
fileName = 'lobulaTerminals_LPLC1_totSynBelow300_minCon3_2WithMorphoStatsLayers.csv'

nComp = 4          # for PCA (this is only for visualization purpose here)
nClst = 15         # number of clusters
nReportTarget = 10
nVisualizeTarget = 30
# preppare colormap
cmap = cm.get_cmap('gist_ncar')

# 1. Read data
df = pd.read_csv(dataPath+fileName)
# find where the actual connectivity data starts in the df
cols = df.columns
dataStart = cols.get_loc("bodyId")+1
postTypes = cols[dataStart:-14]
mat = df.iloc[:,dataStart:].to_numpy()

## 2: Run PCA before clustering (just for visualization)
pca = PCA(n_components=nComp)
PCs = pca.fit_transform(mat)# show variance explained

## 3: Cluter cell
linkage = sch.linkage(mat, method='ward')
dendrogram = sch.dendrogram(linkage, truncate_mode='lastp',p=nClst)
label = sch.fcluster(linkage, nClst, criterion='maxclust')
for ii in np.unique(label): # sanity check that things are ordered correctly
    print('Cluster',ii,np.sum(label==ii),'cells')
colors = (label+1)/(max(label)+1)*2*np.pi

# plot all to all scatter plot
fig2 = plt.figure()
for ii in range(nComp):
    ax = fig2.add_subplot(nComp,nComp,ii*nComp+ii+1)
    ax.hist(PCs[:,ii])
    ax.set_ylabel('PC #'+str(ii+1))
    for jj in range(ii+1,nComp):
        ax = fig2.add_subplot(nComp,nComp,ii*nComp+jj+1)
        if ii==0:
            ax.set_xlabel('PC #'+str(jj+1) )
            ax.xaxis.set_label_position('top')
        for kk in np.unique(label):
            ax.scatter(PCs[label==kk,ii],PCs[label==kk,jj],s=0.5,c=np.array([cmap(kk/(nClst+1))]))

# add a legend panel
ax = fig2.add_subplot(nComp,nComp,nComp*(nComp-1)+1)
for kk in np.unique(label):
    sc = ax.scatter(kk,1,c=np.array([cmap(kk/(nClst+1))]))
ax.set_yticks([])
ax.set_xticks(range(1,nClst+1))
ax.set_title('Clusters')

## 4: Sort original connectivity matrix with cluster label and Visualize
for ii in range(1,nClst+1):
    thisClst = mat[label==ii,:]
    if ii==1:
        sortmat = thisClst
    else:
        sortmat = np.concatenate((sortmat,thisClst),axis=0)

# only visualize important targets for simplicity
conmat = sortmat[:,:-14]
totConPerTarg = np.sum(conmat,axis=0)
importantTargInd = np.argsort(-totConPerTarg)[:nVisualizeTarget]

fig3, ax = plt.subplots()
im = ax.imshow(conmat[:,importantTargInd].T,aspect='auto')

# add border lines
for ii in range(1,nClst):
    plt.plot([np.sum(label<=ii)-0.5,np.sum(label<=ii)-0.5],[-0.5,len(postTypes)-0.5],'w--')
ax.set_ylim(nVisualizeTarget-0.5,-0.5)

# annotate
ax.set_yticks(np.arange(nVisualizeTarget))
ax.set_yticklabels(postTypes[importantTargInd])
ax.tick_params(labelright=True)
cb0 = fig3.colorbar(im,ax=ax)
cb0.set_label('synapses')

# Visualize morphological statistics
fig4 = plt.figure()
ax1 = fig4.add_subplot(2,1,1)
im1 = ax1.imshow(sortmat[:,-14:-11].T,aspect='auto')
for ii in range(1,nClst):
    plt.plot([np.sum(label<=ii)-0.5,np.sum(label<=ii)-0.5],[-0.5,2.5],'w--')

ax2 = fig4.add_subplot(2,1,2)
im2 = ax2.imshow(sortmat[:,-11:].T,aspect='auto')
for ii in range(1,nClst):
    plt.plot([np.sum(label<=ii)-0.5,np.sum(label<=ii)-0.5],[-0.5,10.5],'w--')

ax1.set_ylim(2.5,-0.5)
ax1.set_yticks(np.arange(3))
ax1.set_yticklabels(['std1','std2','std3'])
cb1 = fig4.colorbar(im1,ax=ax1)
cb1.set_label('microns')

ax2.set_ylim(0.5,-0.5)
ax2.set_yticks(np.arange(11))
ax2.set_yticklabels(['L'+str(ii) for ii in range(11)])
cb2 = fig4.colorbar(im2,ax=ax2)
cb2.set_label('synapses')

## 5. Report top target for each cluster
for ii in range(1,nClst+1):
    thisClst = mat[label==ii,:-14]
    totCon = np.sum(thisClst,axis=0)
    mainTargIdx = np.argsort(-totCon)[:nReportTarget]
    print('Main targets of cluster #'+str(ii)+' is: ')
    print(postTypes[mainTargIdx])
    thisMorphoStats = mat[label==ii,-14:]
    meanMS = np.mean(thisMorphoStats,axis=0)
    print('Std1: '+str(meanMS[0]))
    print('Std2: '+str(meanMS[1]))
    print('Std3: '+str(meanMS[2]))
    for ii in range(11):
        print('L'+str(ii)+': '+str(meanMS[3+ii]))

# save sorted result
df['label'] = label
df = df.sort_values(by='label')
df.to_csv(dataPath+fileName[:-4]+'HC'+str(nClst)+'.csv')

plt.show()
