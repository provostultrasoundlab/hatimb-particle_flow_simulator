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

![FIG](/images/Pulsatility_GIF.gif)

![FIG](/images/DEMO_20umMB_150000_total_5000_SSF.gif)
### Parameters used for previous figure:
1. samp_freq = 1000;   
1. n_bubbles = 150000;   
1. n_bubbles_steady_state = 5000;   
1. t_steady_state = 1;   
1. bubble_size = 20;     
1. pulsatility = 1; 
1. bypass_N_vs_d_stats = 1;

## Citing the Particle Flow Simulator
If you use this simulator, please cite:

Preprint:
- [An Anatomically and Hemodynamically Realistic Simulation Framework for 3D Ultrasound Localization Microscopy, 2021](https://www.biorxiv.org/content/10.1101/2021.10.08.463259v1.full.pdf)

@article {Belgharbi2021.10.08.463259,
	author = {Belgharbi, Hatim and Por{\'e}e, Jonathan and Damseh, Rafat and Perrot, Vincent and Milecki, L{\'e}o and Delafontaine-Martel, Patrick and Lesage, Fr{\'e}d{\'e}ric and Provost, Jean},
	title = {An Anatomically and Hemodynamically Realistic Simulation Framework for 3D Ultrasound Localization Microscopy},
	elocation-id = {2021.10.08.463259},
	year = {2021},
	doi = {10.1101/2021.10.08.463259},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2021/10/09/2021.10.08.463259},
	eprint = {https://www.biorxiv.org/content/early/2021/10/09/2021.10.08.463259.full.pdf},
	journal = {bioRxiv}
}

Peer-reviewed:
- (Not yet available)

## Related data
If you use the provided data, please cite:
- [3D Simulated Microbubble Flow in a Mouse Vascular Network](https://www.frdr-dfdr.ca/repo/dataset/ae5706dc-22fb-41c2-be22-082d893fcb9a)

@misc{Belgharbi:2022,
  author = {Belgharbi, Hatim and Porée, Jonathan and Damseh, Rafat and Perrot, Vincent and Milecki, Léo and Delafontaine-Martel, Patrick and Lesage, Frederic and Provost, Jean},
  title = {3D Simulated Microbubble Flow in a Mouse Vascular Network},
  year = {2022},
  howpublished= {10.20383/102.0494}
} 

## Licence
Free software: MIT license