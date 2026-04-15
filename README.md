# prolong-spt

MATLAB tools for single-particle detection and tracking in long-duration live-cell single-molecule movies.
In practice, the pipeline scans an input folder for image files supported by the Bio-Formats toolbox, detects candidate particles, and performs subpixel localization with Gaussian refinement using routines adapted from Picasso. Particle trajectories are then generated using u-track.
## Requirements

### MATLAB
- MATLAB **R2018b or newer** recommended
### MATLAB toolboxes
- Image Processing Toolbox
- Parallel Computing Toolbox
## Installation
1. Clone or download this repository.
2. Add this repository to the MATLAB path.
Example:

```matlab
addpath(genpath('path/to/prolong-spt'));
addpath(genpath('path/to/u-track'));
addpath(genpath('path/to/bfmatlab'));
savepath;
```

## Usage

### `prolong-spt`
#### Function signature

```matlab
prolong-spt(inputFolder, boxSize, minNetGradient, voxelSize, frameInterval)
```

#### Required arguments

- `inputFolder`  
  Path to the folder containing one or more movie files.

- `boxSize`  
  Detection box size in pixels. It must be an odd integer. If an even value is supplied, the function normalizes it to the next odd integer.

- `minNetGradient`  
  Minimum net-gradient threshold used to filter detected particle candidates after the initial detection step.

- `voxelSize`  
  Pixel size in **micrometers per pixel**.

- `frameInterval`  
  Time between frames in **seconds**.
#### Minimal example

```matlab
inputFolder = 'D:\movies';
boxSize = 7;
minNetGradient = 2000;
voxelSize = 0.1;     % um/pixel
frameInterval = 0.05;  % s/frame

spt_headless_mode(inputFolder, boxSize, minNetGradient, voxelSize, frameInterval);
```

### `spt_dataset_export`

Builds a dataset-level summary from a pipeline `Results` folder. The function scans each immediate subfolder for `processingInfo.mat`, aggregates track statistics across datasets, exports summary tables, and generates histogram and trajectory-map outputs.
#### Function signature

```matlab
spt_dataset_export(Results_folderpath)
```
#### Arguments
- `Results_folderpath`  
  Path to a `Results` folder. Each immediate subfolder is treated as one dataset entry and should contain `processingInfo.mat`.
#### Minimal example
```matlab
Results_folderpath = 'D:\movies\Results';
spt_dataset_export(Results_folderpath);
```

#### Exported statistics
The export summarizes per-dataset tracking statistics, including:
- diffusion coefficient
- travel distance
- travel time
- major axis
- minor axis
- track count
#### Exported files
In the main `Results` folder, the function writes:
- `datasetInfo.mat` — MATLAB structure containing raw and mean dataset statistics
- `datasetInfo.csv` — CSV summary file with one row per dataset subfolder
- `D_hist.jpg` and `D_hist.pdf` — diffusion coefficient histogram
- `Major_axis.jpg` and `Major_axis.pdf`
- `Minor_axis.jpg` and `Minor_axis.pdf`
- `Travel_distance.jpg` and `Travel_distance.pdf`
- `Travel_time.jpg` and `Travel_time.pdf`
- `traj_map.png` — combined trajectory-map grid

In each dataset subfolder, the function may also write:
- `traj_map.fig` — per-movie trajectory map

#### Notes
- Histogram export requires `histogram_dataset_log10` to be available on the MATLAB path.
- Trajectory-map export requires `traj_map_plot` to be available on the MATLAB path.
- If plotting helpers are missing, `datasetInfo.mat` and `datasetInfo.csv` are still exported.
- The combined `traj_map.png` includes up to the first 50 subfolders in folder order.
### Output structure
- For each batch run, the function creates a results directory inside the input folder.
- Inside the results directory, the pipeline writes:

```text
Results/
├── NNB_log.txt
├── movie_1/
|   └── processingInfo.mat
└── movie_2/
    └── processingInfo.mat
```

#### `NNB_log.txt`

A plain-text batch log containing:

- run start information
- input folder and parameter summary
- per-file progress messages
- warnings and error messages
- final batch summary

#### `processingInfo.mat`

A per-movie MATLAB file that stores run metadata and analysis outputs. Depending on the stage reached successfully, it may contain:

- input settings used for the run
- source image path and image metadata
- particle coordinates and subpixel localization results
- number of detected particles
- particle tracking outputs
- error information if the movie failed during analysis

### Additional outputs from `spt_dataset_export`

After running `spt_dataset_export(Results_folderpath)`, the `Results` folder may also contain:

```text
Results/
├── datasetInfo.mat
├── datasetInfo.csv
├── D_hist.jpg
├── D_hist.pdf
├── Major_axis.jpg
├── Major_axis.pdf
├── Minor_axis.jpg
├── Minor_axis.pdf
├── Travel_distance.jpg
├── Travel_distance.pdf
├── Travel_time.jpg
├── Travel_time.pdf
├── traj_map.png
├── movie_1/
│   ├── processingInfo.mat
│   └── traj_map.fig
└── movie_2/
    ├── processingInfo.mat
    └── traj_map.fig
```

### `datasetInfo.mat`

A MATLAB file containing aggregated dataset-level outputs, including raw per-track values and mean statistics across subfolders.

### `datasetInfo.csv`

A CSV summary table containing one row per dataset subfolder. Exported columns include:

- `D (µm²/s)`
- `travelDistances (µm)`
- `travelTime (s)`
- `majorAxis (µm)`
- `minorAxis (µm)`
- `trackNo`
- `immobile fraction` *(only when `immobile_threshold` is provided)*

### Histogram and trajectory-map files

The export function can also generate dataset-level histogram files for diffusion coefficient, major axis, minor axis, travel distance, and travel time, along with per-subfolder trajectory maps and a combined trajectory-map overview image.
## Troubleshooting

### `Missing required functions on the MATLAB path`

One or more required project or dependency functions cannot be found. Re-check your MATLAB path and confirm that `u-track`, `bfmatlab`, and this repository's function folders are available.

### `Bio-Formats Java classes are still unavailable after setup`

MATLAB could not initialize Bio-Formats correctly. Confirm that:

- MATLAB is running with the JVM enabled
- `bioformats_package.jar` exists inside the `bfmatlab` folder
- the `bfmatlab` folder is either near the repository root or added to the MATLAB path

### `Unsupported4DStack`

The input file contains both multiple Z slices and multiple time points. Export or pre-process the movie into a supported 2D time series before running the pipeline.

### `no centers found` or `no centers after filtering`

The current detection settings may be too strict for the dataset. Re-evaluate:

- `minNetGradient`
- acquisition signal-to-noise ratio
- labeling density
- movie quality and background level

## License
- This repository is distributed under the **MIT License**. See `LICENSE` for the complete license text.
- Third-party dependencies, including **u-track** and **Bio-Formats**, remain subject to their respective licenses and are not relicensed by this repository. **Picasso** is distributed under the **MIT License**, and **u-track** is distributed under the **GPL-3.0 License**.
