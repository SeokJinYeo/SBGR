# SBGR

# Spot-Based Global Registration (SBGR)

## Overview

SBGR is a novel strategy for achieving precise image stitching at the single-molecule level in tissue-scale spatial transcriptomics. This method shifts the focus from raw images to identified molecular spots, enabling high-resolution image alignment with exceptional accuracy.

## Features

- **Spot-Based Data:** SBGR utilizes identified molecular spots for image alignment, offering a robust evaluation of translation estimates and stitching performance.
  
- **Sub-pixel Accuracy:** Achieves sub-pixel accuracy (83 ± 36 nm) outperforming existing image-based stitching methods.
  
- **Duplicate Spot Removal:** Incorporates a mechanism to surgically remove duplicate spots in overlapping regions, maximizing information recovery.



# : SBGR is Fast and accurate Matlab based 'susceptibility-induced B0 inhomogeneity' calculation programs
ppm2Hz is based on generalized susceptibility voxel convolution (gSVC) method which is rapid and artifact-free[1].
For application in various cases, the method was extended to arbitrary orientations and spatially varying applied fields[2].
Additionally another static magnetic field perturbation calculation method, k-space-discretized (KD) is coded in Matlab (ppm2Hz_KD). If you want more details (theory, applications) regarding the method, see the references below.

# Usage examples


# References
[1] Seokjin Yeo, Alex W Schrader, Juyeon Lee, Marisa Asadian, Hee-Sun Han, Spot-Based Global Registration for Sub-pixel Resolution Stitching of Single-Molecule Resolution Images for Tissue-Scale Spatial Transcriptomics.
