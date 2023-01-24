import numpy as np
import pandas
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm as cm
import sys
import math
import os
#np.set_printoptions(threshold=np.nan)


def vector(p1, p2):
    return [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]

def vabs(a):
    return math.sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2))
def checkIfNative(xyz_CAi, xyz_CAj):
    v = vector(xyz_CAi, xyz_CAj)
    r = vabs(v)
    if r<12.0: return True
    else: return False

def isNative(r):
        if r<12.0: return True
        else: return False

def get_ca_s_atoms(pdb_part):
    ca_atoms = []
    for line in pdb_part:
        if line.split()[0] != 'END':
           if line.split()[2] == 'CA' or line.split()[2] == 'S':
           #print line
              x=float(line[30:38])
              y=float(line[38:46])
              z=float(line[46:54])
              atom = [x,y,z]
              ca_atoms.append(atom)
   # print coord
    return ca_atoms

def compute_frustmap(pmffile,n_res):
    width = n_res
    frustmap = np.zeros((width,width))
    with open(pmffile,"r") as fopen:
         for line in fopen.readlines():
             if len(line) > 60:
                 if line[0:10]=="draw_links":
                   index1=int(line.split(",")[0].split()[2])
                   index2=int(line.split(",")[1].split()[1])
                   if (line.split()[-1][-1]) == "n":
                    frustmap[index1][index2] = 1
                    frustmap[index2][index1] = 1
                   if (line.split()[-1][-1]) == "d":
                    frustmap[index1][index2] = -1
                    frustmap[index2][index1] = -1
         print (frustmap)
    return frustmap

def merge_matrix(contactmap1,contactmap2,n_atoms):
    matrix = np.zeros((n_atoms,n_atoms))
    for i in range(n_atoms):
        for j in range(i+1,n_atoms):
            if contactmap1[i][j] == 1 and contactmap2[i][j] != 1: 
               matrix[i][j] = 1
            elif contactmap1[i][j] == -1 and contactmap2[i][j] != -1:
               matrix[i][j] = -1
            matrix[j][i] = matrix[j][i]
    return matrix

def draw_contactmap(matrix):
    (n_atomsx,n_atomsy) = np.shape(matrix)
    if n_atomsx == n_atomsy:
       n_atoms = n_atomsx
    else:
       print  ("matrix is not square")
       sys.exit()
    plt.figure(figsize=(10,10))
    X=np.arange(1,n_atoms+1);
    cmap = cm.get_cmap('RdYlGn')
    Y=np.arange(1,n_atoms+1);
    #print matrix
    #cmax = np.max(matrix)*1.2
    for i in range(n_atomsx):
       for j in range(i+5,n_atomsx):
            if matrix[i][j] == 1:
               plt.plot([i],[j],'gs')
               plt.plot([j],[i],'gs')
            if matrix[i][j] == -1:
               plt.plot([i],[j],'rs')
               plt.plot([j],[i],'rs')
    #for i in range(int(n_atoms)):
    #  for j in range(int(n_atoms)):
      #   print matrix
      #   print i,j,matrix[i][j]
    #     if matrix[int(i)][int(j)] == 0.0:
    #         matrix[i][j]=float("inf")
    #ax = plt.pcolor(X,Y,matrix,cmap=cmap)
    #,edgecolors='k')
#colorbar;
    #plt.colorbar()
    #cmax = max(matrix)
    #plt.clim(-1,1)
    #plt.xlabel(xname,fontsize=30)
    #plt.ylabel(yname,fontsize=30)
    #plt.title(title,fontsize=30)
    major_ticks = np.arange(1, n_atoms+1,8)
    minor_ticks = np.arange(1, n_atoms+1,8)
    #ax = plt.axes()
    #ax.set_xticks(major_ticks)
    #ax.set_xticks(minor_ticks, minor=True)
    #ax.set_yticks(major_ticks)
    #ax.set_yticks(minor_ticks, minor=True)
#plt.grid()
    ax=plt.gca()
    #ax.set_yticks([72,124,368,441,491,736], minor=False)
    #ax.set_yticks([72,124,368,441,491,736], minor=True)
    #ax.yaxis.grid(True, which='major')
    #ax.yaxis.grid(True, which='minor')
    #ax.set_xticks([72,124,368,441,491,736], minor=False)
    #ax.set_xticks([72,124,368,441,491,736], minor=True)
    #ax.xaxis.grid(True, which='major')
    #ax.xaxis.grid(True, which='minor')
    #plt.axis([0,n_atoms+1, 0, n_atoms+1])
    plt.xticks(rotation="vertical",fontsize=20)
    plt.yticks(fontsize=20)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig('frust.png')
    plt.show()

def main():
    #if len(sys.argv) != 6:
    #   print "*py txtfile pdbfile xname yname title"
    #   sys.exit()
    print "*py txtfile pmlfile1 pmlfile2 #residue "
    print "two protein should be same length"
    print "green: minimally frustrated in pmlfile1 but not in pmlfile2; red: highly frustrated in pmlfile1 but not in pmlfile2"
    pmlfile1 = sys.argv[1]
    pmlfile2 = sys.argv[2]
    n_res = int(sys.argv[3])
    #xname = sys.argv[3]
    #yname = sys.argv[4]
    #title = sys.argv[5]
    #contactmap2 = np.loadtxt(txtfile2)
    #with open(pdbfile,'r') as fopen:
    #     pdb_part = fopen.readlines()
    #ca_atoms = get_ca_s_atoms(pdb_part)
    frustmap1 =  compute_frustmap(pmlfile1,n_res)
    frustmap2 =  compute_frustmap(pmlfile2,n_res)
    matrix = merge_matrix(frustmap1,frustmap2,n_res) 
    draw_contactmap(matrix)
    
if __name__ == '__main__':
    main()

