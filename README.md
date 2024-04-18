# Spot-Based Global Registration (SBGR)

## Overview

SBGR is a novel strategy for achieving precise image stitching at the single-molecule level in tissue-scale **spatial transcriptomics**. This method shifts the focus from raw images to identified molecular spots, enabling high-resolution image alignment with exceptional accuracy.

![SBGR workflow](https://github.com/SeokJinYeo/SBGR/blob/main/Wrokflow.png)
## Features

- **Spot-Based Data:** SBGR utilizes identified molecular spots for image alignment, offering a robust evaluation of translation estimates and stitching performance.
  
- **Sub-pixel Accuracy:** Achieves sub-pixel accuracy outperforming existing image-based stitching methods.
  
- **Duplicate Spot Removal:** Incorporates a mechanism to surgically remove duplicate spots in overlapping regions, maximizing information recovery.


## Input Parameters

- **`pixel_size`**: The physical size of one pixel in the image, specified in micrometers (µm). This parameter is essential for converting pixel measurements to real-world units.

- **`output_loc`**: The file path to the directory where output files from the SBGR function will be stored.

- **`spots_loc`**: The path to the CSV file containing information about detected spots (`barcodes_20230426_raw_warped.csv`).

- **`positions_loc`**: The path to the CSV file containing the locations of each FOV before stitching (`positions.csv`).

- **`z_step`**: The step size along the z-axis for images taken at different depths, specified in micrometers (µm).

## Input Files

### `Spot file (ex. barcodes_20230426_raw_warped.csv)`

Contains detailed information about detected spots from spatial transcriptomics data, with the following columns:

- **`barcode_id`**: Unique identifier for each spot's barcode (gene identity).
- **`fov`**: Index of the FOV where the spot was detected.
- **`global_x`, `global_y`, `global_z`**: Global coordinates of the spot, mapping it onto the entire sample or experiment's reference frame.

### `Position file (ex. positions.csv)`

A two-column CSV file with each row representing the x and y coordinates (in µm) of the top left corner of a FOV. This file is crucial for determining the absolute positions of the FOVs within the global space of the sample.

## Output files

The SBGR function will output below data sets:

- **`stitched_spots.csv`**: Contains the processed spots data, including their new global positions after stitching.

- **`stitched_positions.csv`**: Contains the adjusted positional information of the FOVs.

- **`offset.csv`**: Contains the d_err and the number of total spots and duplicate spots etc.. for each overlapping region.

- **`dup_spots.csv`**: Contains the identities of every identified duplicate spots.

## Usage Notes

Please ensure the input CSV files are correctly formatted with the expected headers and data types. The output directory must be a valid path with the necessary write permissions.

## Usage examples
Put **SBGR.py** on your repository. Edit the inputs.
```
from SBGR import SBGR
if __name__ == '__main__':  
    pixel_size = 0.103 #um
    fov_size = 2048 #pixel
    spots_loc =  '/barcodes_20230426_raw_warped.csv'
    positions_loc = '/positions.csv'
    output_loc = r'/JIn\SBGR'
    SBGR(spots_loc,positions_loc,pix_size,fov_size, outout_loc)
```


# References
[1] Seokjin Yeo, Alex W Schrader, Juyeon Lee, Marisa Asadian, Hee-Sun Han, "Spot-Based Global Registration for Sub-pixel Resolution Stitching of Single-Molecule Resolution Images for Tissue-Scale Spatial Transcriptomics."
