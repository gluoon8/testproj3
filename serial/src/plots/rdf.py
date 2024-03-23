# RDF for Molecular Dynamics Simulation - Serial Version
#
# Manel Serrano 
#
#
"""This script is used to plot the radial distribution function (RDF) of a molecular dynamics simulation. 
The RDF is a histogram that shows the probability of finding a particle at a distance r from a reference particle. 
The RDF is calculated as follows:
    - We calculate the distance between all pairs of particles in the simulation box.
    - We bin the distances in a histogram.
    - We normalize the histogram by the number of particles and the volume of the box.
    - We plot the RDF.
"""
# Import libraries
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy as sp

# Set current directory as default path
path = "/home/manel/Project-III/serial"
os.chdir(path)


# Load data from files
positions = np.loadtxt('pos_ini.dat', dtype=float)
print(positions.shape)

L = 5.6312394

def distance(a, b):
    dx = abs(a[0] - b[0])
    x = min(dx, abs(L - dx))
     
    dy = abs(a[1] - b[1])
    y = min(dy, abs(L - dy))
     
    dz = abs(a[2] - b[2])
    z = min(dz, abs(L - dz))
 
    return np.sqrt(x**2 + y**2 + z**2)


n = len(positions)

def rdf(positions, L, dr):
    n = len(positions)
    hist = np.zeros(int(L/(2*dr)))
    for i in range(n):
        for j in range(i+1, n):
            r = distance(positions[i], positions[j])
            hist[int(r/(2*dr))] += 2
    return hist

dr = 0.1
rdf = rdf(positions, L, dr)
r = np.arange(0, L/2, dr)[:len(rdf)]  # Adjust the range of np.arange() to match the length of rdf

print(r, rdf)

g = rdf/(4*np.pi*r**2*dr*n)
plt.plot(r, g)
plt.show()






