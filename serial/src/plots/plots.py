import numpy as np
import os

import matplotlib.pyplot as plt

# Set current directory as default path
path = "/Users/manel/EIA-G3-1/src/plots"

os.chdir(path)

data = np.loadtxt('thermodynamics.dat', skiprows=3)

# Replace inf values with NaN
data[data == np.inf] = np.nan
data[data == -np.inf] = np.nan

# Plot the data with matplotlib

#
#       PLOT ENERGIES vs STEP
#

plt.figure(figsize=(10, 5))
plt.plot(data[:, 0], data[:, 1], label='Potential Energy')
plt.plot(data[:, 0], data[:, 2], label='Kinetic Energy')
plt.plot(data[:, 0], data[:, 3], label='Total Energy')
plt.xlabel('Step')
plt.ylabel('Energy')
plt.legend()
plt.title('Energy vs Step')
plt.show()

#
#      PLOT TEMPERATURE vs STEP
#

plt.figure(figsize=(10, 5))
plt.xlabel('Step')
plt.ylabel('Temperature')
plt.legend()
plt.title('Temperature vs Step')
plt.show()
