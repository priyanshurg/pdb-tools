#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 22:06:09 2023

@author: pgupta
"""


# =============================================================================
# Q2
from Bio.PDB.internal_coords import * 
from Bio.PDB.ic_rebuild import write_PDB, IC_duplicate, structure_rebuild_test 
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.internal_coords import * 
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import PDBList
import numpy as np
import math
pdb_list = PDBList()
ppb = PPBuilder()
pdb_id = input("Enter RCSB code (ex. 1bty for Trypsin, 1u7g for AmtB): \n")
pdb_filename = pdb_list.retrieve_pdb_file(pdb_id, pdir="pdb_files/", file_format="pdb")

userinput = int(input("Choose 1 to fix bin width, or Choos2 to fix number of bins (Ex. 12bins or 30deg bin width): \n"))
if userinput == 1:
    binwidth = int(input("Enter bin width desrired: \n"))
    nbins = 360/binwidth
elif userinput==2:
    nbins = int(input("Enter no. of bins desired: \n"))
           
phi_bins = np.linspace(-180,180,nbins+1)
psi_bins =np.linspace(-180,180,nbins+1)
def get_phi_psi_angles(chain):
    # Initialize lists to store phi and psi angles
    phi_angles = []
    psi_angles = []    
    IC_Chain.MaxPeptideBond = 4.0
    peptide = ppb.build_peptides(chain)
    for pp in peptide:
        angle_list = pp.get_phi_psi_list()
        for i,j in enumerate(angle_list):
            if i == 0 or i == len(angle_list):
                phi_angles.append(angle_list[i][0]) 
                psi_angles.append(angle_list[i][1]*180/math.pi)
            elif i == len(angle_list)-1:
                phi_angles.append(angle_list[i][0]*180/math.pi) 
                psi_angles.append(angle_list[i][1])
            else:
                phi_angles.append(angle_list[i][0]*180/math.pi) 
                psi_angles.append(angle_list[i][1]*180/math.pi)

    return phi_angles, psi_angles

parser=PDBParser()
structure = parser.get_structure('STS', 'pdb_files/pdb'+pdb_id+".ent")
a=list(structure.get_chains())[0]
phi_list, psi_list = get_phi_psi_angles(a)
phi_list2 = []
psi_list2 = []
for i in range(0,len(phi_list)):
    if [phi_list != np.array(None)][0][i]:
        phi_list2.append(phi_list[i])
    if [psi_list != np.array(None)][0][i]:
        psi_list2.append(psi_list[i])

phi_store = [[] for _ in range(nbins)]
psi_store = [[] for _ in range(nbins)]

for phi, psi in zip(phi_list,psi_list):
    for i in range(0,nbins):
        if phi:
            if phi >= phi_bins[i] and phi <= phi_bins[i+1]:
                phi_store[i].append(phi)
                continue
        if psi:
            if psi >= psi_bins[i] and psi <= phi_bins[i+1]:
                psi_store[i].append(psi)
                continue

phi_store_normalized = []                
length_of_peptide = max(len(phi_list),len(psi_list))
ram_phi0 = []
ram_psi0 = []
for i in range(0,nbins):
    ram_phi0.append(len(phi_store[i])/length_of_peptide)
    ram_psi0.append(len(psi_store[i])/length_of_peptide)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.set_xticks(phi_bins)

ax.grid(True)
ax.hist(phi_list2,bins =phi_bins,density=True,label="phi-angles")
ax.hist(psi_list2,bins =psi_bins, density=True,label="psi-angles")

h, xedges, yedges = np.histogram2d(phi_list2, psi_list2, bins=(phi_bins, psi_bins), density=True)
fig.legend()
xbin = np.linspace(-180, 180, nbins)  # phi axis for contour plot
ybin = np.linspace(-180, 180, nbins)  # psi axis

hval = h.T  # transpose the count matrix since contour flips the column & row

fig2, ax2 = plt.subplots()
ax2.grid(False)
contour= ax2.contourf(xbin, ybin, hval)
cbar = plt.colorbar(contour) # show the scale bar

plt.show()
   




 
# =============================================================================
