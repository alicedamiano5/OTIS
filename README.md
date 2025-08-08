# OTIS
The Orbital Timescale for Sinking Python Library

# Description

This library calculates the orbital sinking timescales for BH infalling into cosmic structures according to Chandrasekhar dynamical friction formula. The code is designed to integrate the equations of motion of a BH of arbitrary mass and initial positon in a NFW halo eventually provided with a central stellar bulge. 
Here is the a scatchy flowchart of the library:
![image](https://github.com/user-attachments/assets/03f474e7-9268-45ca-8b10-23ddd07944c5)
The code workflow starts with the initialization of the DM halo parameters, including the concentration parameter and virial quantities, eventually followed by the initialization of the stellar Hernquist bulge. The density, mass, and velocity dispersion profiles of both the halo and the bulge can be hence retrieved. The BH is then initialized with its mass and phase-space coordinates. The integral of the velocity distribution function is precomputed on a phase-space grid, enabling fast interpolation during the numerical integration of the equations of motion. User-configurable options allow customization of the interpolation technique, integration method, time steps, and impact parameters.

 
# How to install OTIS

Dependencies: SciPy, NumPy

OTIS is easy to install using pip:

```
pip3 install https://github.com/alicedamiano5/OTIS.git
```

Then just import otis in you python code: 
```
import otis

```
hand have fun!


# Example code

In the tests directory you will find an example to see OTIS in action !


# Integrator accuracy and default options

The accuracy of the integration depends on the specific Scipy integrator and grid size. To reduce the integrator inaccuracies, OTIS default maximum timestep size is 0.01 and the default integrator is the Scipy Radau. Integrators ad timestep sizes critically depend on the choice of the problem, when moving from the default option, always double-check your results using a different integrators. 
The following figure analyses the time to solution as a function of the specific Scipy integrator adopted, with bar color-coded with the number of acceleration evaluation. Each group of bars corresponds to a different maximum allowed timestep, in Gyr: default (left), 0.05 (middle), and 0.01 (right). The number on top of each bar indicates the number of steps taken. The corresponding distance to the halo center is displayed in the bottom figure. For the specific problem analysed, the high accuracy integrator RK45 (suitable fro non stiff problems) leads to the wrong solution if the maximum timestep size is not constrained. On the other hand, the Radau integrator shows stability across the the timestep sizes. 

[rr_timetosolution_threeplots.pdf](https://github.com/user-attachments/files/21686507/rr_timetosolution_threeplots.pdf)

[rr_Orbits.pdf](https://github.com/user-attachments/files/21686557/rr_Orbits.pdf)



