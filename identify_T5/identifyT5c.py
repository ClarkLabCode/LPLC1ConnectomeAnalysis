# This script identifies candidate T5c based on their synapse locations and the lack
# of their connectivity to LPTCs (there is no identified monostratified LPTC in the
# layer 3 of lobula plate right now). We also use the statistics of candidate T5a,
# b, and d subtypes to get the sense of the size of T5c neurons.
# The list of putative T5 neurons after visual annotation is provided as a supplementary
# file of the manuscript.

## load packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

## Define parameters
nonTargetList = ['HSE','HSN','HSS','VCH','DCH','H2','VS'] # T5c cannot be connected to these cells

## Connect to the server
# NOTE: For this script to run, you must put "authtoken" file under the top directory
f = open("authtoken","r")
tokenstr = f.read()
c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token=tokenstr)
c.fetch_version()

## Read the results of T5a,b,d identification as a reference
referenceTypes = ('a','b','d')

list_pre_Lo = []
list_pre_Lp = []
list_post_Lo = []
list_post_Lp =[]

for type in referenceTypes:
    df = pd.read_csv('./identify_T5/results/candidateT5'+type+'_annotated.csv')
    thisIdList = df['bodyId'].to_list()
    q = """\
        MATCH (a:Neuron)
        WHERE a.bodyId IN %s
        RETURN DISTINCT a.bodyId as bodyId, apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].pre as LOpre, apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].post as LOpost,
        apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].pre as LOPpre, apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].post as LOPpost
        """ % thisIdList
    df = c.fetch_custom(q)
    list_pre_Lo.append(df['LOpre'].to_list())
    list_pre_Lp.append(df['LOPpre'].to_list())
    list_post_Lo.append(df['LOpost'].to_list())
    list_post_Lp.append(df['LOPpost'].to_list())

## Get the maximum numbers of synapses T5a, b, d cells have in lobula/lobula plate
# Our candidate T5c cells should fit withint these maxima
max_lo_pre = np.max([N for sublist in list_pre_Lo for N in sublist])
max_lp_pre = np.max([N for sublist in list_pre_Lp for N in sublist])
max_lo_post = np.max([N for sublist in list_post_Lo for N in sublist])
max_lp_post = np.max([N for sublist in list_post_Lp for N in sublist])
print('max number of (presynapse in LO, presynapse in LP, postsynapse in LO, postsynapse in LP) for T5a,b,d = ',(max_lo_pre,max_lp_pre,max_lo_post,max_lp_post))

## get LO/LOP intrinsic cell with less synapses than these maxima
# also they need to have more inputs in lobula and more outputs in lobula plate
q = """\
    MATCH (a:Neuron)
    WHERE apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].pre+apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].pre=a.pre
    AND apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].post+apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].post=a.post
    AND apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].pre>apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].pre
    AND apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].post>apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].post
    AND apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].pre<%s AND apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].pre < %s
    AND apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].post < %s AND apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].post < %s
    RETURN DISTINCT a.bodyId as bodyId
    """ % (max_lo_pre,max_lp_pre,max_lo_post,max_lp_post)
df = c.fetch_custom(q)
candidateT5cList = df['bodyId'].to_list()
print(len(candidateT5cList))

## convert no-target cell types into a list of bodyIds
nonoList = []
for thisType in nonTargetList:
    q = """\
        MATCH (a:Neuron)
        WHERE a.type='%s'
        RETURN DISTINCT a.bodyId as bodyId
        """ % thisType
    df = c.fetch_custom(q)
    nonoList.extend(df['bodyId'].to_list())

## Find ones that have unwanted cross connectivity
q = """\
    MATCH (a:Neuron)-[:ConnectsTo]->(b:Neuron)
    WHERE a.bodyId IN %s and b.bodyId IN %s
    RETURN DISTINCT a.bodyId as bodyId
    """ % (candidateT5cList,nonoList)
df = c.fetch_custom(q)
candidateT5cList = list(set([ids for ids in candidateT5cList if ids not in df['bodyId'].to_list()]))
print('#exclusively connected cells found',len(candidateT5cList))

# reformat the output
outdf = pd.DataFrame(columns=['bodyId'],data=candidateT5cList)
outdf.to_csv('./identify_T5/results/candidateT5c.csv')
