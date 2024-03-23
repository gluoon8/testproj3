import numpy as np
import os
import matplotlib.pyplot as plt

# Set current directory as default path
path = "/home/manel/Project-III/serial"
os.chdir(path)


# Load data from files
energy = np.loadtxt('energy_verlet.dat', skiprows=4, dtype=float)
temperature = np.loadtxt('Temperatures_verlet.dat', skiprows=1, dtype=float)


#
#       PLOT ENERGIES vs time
#

plt.figure(figsize=(10, 5))
plt.plot(energy[:, 0], energy[:, 1], label='Potential Energy')
plt.plot(energy[:, 0], energy[:, 2], label='Kinetic Energy')
plt.plot(energy[:, 0], energy[:, 3], label='Total Energy')
plt.xlabel('Step')
plt.ylabel('Energy')
plt.legend()
plt.title('Energy vs Step')
plt.show()


#
#      PLOT MOMENTUM vs time
#

plt.figure()
plt.plot(energy[:, 0], energy[:, 4], label='Momentum')
plt.ylim(np.mean(energy[:, 4])-1, np.mean(energy[:, 4])+1)
plt.xlabel('Time')
plt.ylabel('Momentum')
plt.legend()
plt.title('Momentum vs Time')
plt.savefig('Momentum.png')
plt.show()


#
#      PLOT TEMPERATURE vs time
#

plt.figure()
plt.plot(temperature[:, 0], temperature[:, 1], label='Temperature')
plt.xlabel('Time')
plt.ylabel('Temperature')
plt.legend()
plt.title('Temperature vs Time')
plt.savefig('Temperature.png')
plt.show()


#
#      PLOT PRESSURE vs time
#

plt.figure()
plt.plot(temperature[:, 0], temperature[:, 2], label='Pressure')
plt.xlabel('Time')
plt.ylabel('Pressure')
plt.legend()
plt.title('Pressure vs Time')
plt.savefig('Pressure.png')
plt.show()
