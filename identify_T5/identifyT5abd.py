# This script identifies candidate T5a, b, d subtypes based on their synapse locations
# as well as their connectivity to identified LPTCs in the hemibrain dataset.
# Candidate T5 cells are defined as cells that have exclusive connectivity to specific
# sets of monostratified LPTCs, are intrinsic to lobula / lobula plate, and have
# more inputs in lobula and more outputs in lobula plate.
# The csv list of putative T5 neurons after visual annotation is provided as a
# part of the result folder (which is used for T5c identification)

## load packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

## Define parameters
targetTypes = [['HSE','HSN','HSS','VCH','DCH'],['H2'],['VS']] # T5a is connected to HS/CH, T5b to H2, and T5d to VS
candidateT5list = []

## Connect to the server
# NOTE: For this script to run, you must put "authtoken" file under the top directory
f = open("authtoken","r")
tokenstr = f.read()
c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token=tokenstr)
c.fetch_version()

## Loop through the three sets of postsynaptic LPTC types
allIdList = []
for ii in range(3):
    thisList = targetTypes[ii]
    IdList = []

    # Create a list of bodyIds of all the cells that belong to the set of cell types
    for thisType in thisList:
        q = """\
            MATCH (a:Neuron)
            WHERE a.type='%s'
            RETURN DISTINCT a.bodyId as bodyId
            """ % thisType
        df = c.fetch_custom(q)
        IdList.extend(df['bodyId'].to_list())
    # This is a nested list of lists of bodyIds (used later to make sure connections are exclusive)
    allIdList.append(IdList)

    # Find all the neurons that are connected to specific LPTCs, are intrinsic to lobula / lobula plate,
    # have more inputs in lobula and more outputs in lobula plate
    q = """\
        MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron)
        WHERE b.bodyId IN %s AND apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].pre+apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].pre=a.pre
        AND apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].post+apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].post=a.post
        AND apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].pre>apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].pre
        AND apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].post>apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].post
        RETURN a.bodyId as bodyId, apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].pre as LOpre, apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].post as LOpost,
        apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].pre as LOPpre, apoc.convert.fromJsonMap(a.roiInfo)["LOP(R)"].post as LOPpost
        """ % IdList
    df = c.fetch_custom(q)
    candidateT5list.insert(ii,df['bodyId'].to_list())

# Among identified T5 cells, remove ones that have connection to more than one sets
# of LPTCs
for ii in range(3):
    thisList = candidateT5list[ii]
    nonoList = [ids for sublist in allIdList for ids in sublist if ids not in allIdList[ii]]
    q = """\
        MATCH (a:Neuron)-[:ConnectsTo]->(b:Neuron)
        WHERE a.bodyId In %s AND b.bodyId In %s
        RETURN DISTINCT a.bodyId as bodyId
        """ % (thisList, nonoList)
    df = c.fetch_custom(q)
    badList = df['bodyId'].to_list()
    candidateT5list[ii] = list(set([ids for ids in thisList if ids not in badList]))
    print('#exclusively connected cells found',len(candidateT5list[ii]))

# save the results by subtypes
type = ('a','b','d')
for ii in range(3):
    outdf = pd.DataFrame(columns=['bodyId'],data=candidateT5list[ii])
    outdf.to_csv('./identify_T5/results/candidateT5'+type[ii]+'.csv')
