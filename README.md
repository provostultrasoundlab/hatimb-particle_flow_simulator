#       ⚪⚪⚪ Particle Flow Simulator ⚪⚪⚪
This simulator takes a .swc tree graph as input ans outputs  3D particles positions as trajectories and steady state frames.

## [Input] 
.swc file containing a tree graph with source and target positions with vessel radii

## [Output] 
### Particle Trajectories
Individual 3D particles trajectories in a cell variable stored in a .mat file
### Steady State Frames
Temporal frames in a steady state flow (constant amount of particles in time) saved in a .mat file. 

Here are the columns signification: 
1. Particle ID from the trajectories file
1. Particle position index in it's trajectory from trajectory file
1. Particle X spatial position (um)
1. Particle Y spatial position (um)
1. Particle Z spatial position (um)

![Main Image](/images/bubbles_20umDiameter_5000.png)
