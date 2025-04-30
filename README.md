# OTIS
The Orbital Timescale for Sinking Python Library

# Description

This library calculates the orbital sinking timescales for BH infalling into cosmic structures according to Chandrasekhar dynamical friction formula. The code is designed to integrate the equations of motion of a BH of arbitrary mass and initial positon in a NFW halo eventually provided with a central stellar bulge. 
Here is the a scatchy flowchart of the library:
![image](https://github.com/user-attachments/assets/03f474e7-9268-45ca-8b10-23ddd07944c5)
The code workflow starts with the initialization of the DM halo parameters, including the concentration parameter and virial quantities, eventually followed by the initialization of the stellar Hernquist bulge. The density, mass, and velocity dispersion profiles of both the halo and the bulge can be hence retrieved. The BH is then initialized with its mass and phase-space coordinates. The integral of the velocity distribution function is precomputed on a phase-space grid, enabling fast interpolation during the numerical integration of the equations of motion. User-configurable options allow customization of the interpolation technique, integration method, time steps, and impact parameters.

 

