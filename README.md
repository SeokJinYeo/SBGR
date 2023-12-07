# Spot-Based Global Registration (SBGR)

## Overview

SBGR is a novel strategy for achieving precise image stitching at the single-molecule level in tissue-scale **spatial transcriptomics**. This method shifts the focus from raw images to identified molecular spots, enabling high-resolution image alignment with exceptional accuracy.

![SBGR workflow](https://github.com/SeokJinYeo/SBGR/blob/main/Wrokflow.png)
## Features

- **Spot-Based Data:** SBGR utilizes identified molecular spots for image alignment, offering a robust evaluation of translation estimates and stitching performance.
  
- **Sub-pixel Accuracy:** Achieves sub-pixel accuracy outperforming existing image-based stitching methods.
  
- **Duplicate Spot Removal:** Incorporates a mechanism to surgically remove duplicate spots in overlapping regions, maximizing information recovery.


# Usage examples
Put SBGR.py on your repository.

All packages in the SBGR.py are pre-requirements.

Prepare decoded spots file and microscope position file.

Make sure that files have same format with barcode_warped.csv, positions.csv in exmample data

```
import SBGR
stitched_spots,stitched_positions = SBGR.SBGR(pix_size, save_folder, spot_loc, position_loc, edge_pix, z_step)
```


# References
[1] Seokjin Yeo, Alex W Schrader, Juyeon Lee, Marisa Asadian, Hee-Sun Han, "Spot-Based Global Registration for Sub-pixel Resolution Stitching of Single-Molecule Resolution Images for Tissue-Scale Spatial Transcriptomics."
