# mdCATH Dataset Processing Project for RMSF Prediction

## Overview

This project provides a comprehensive data processing pipeline for the mdCATH protein dynamics dataset. It extracts, transforms, and organizes the data into formats optimized for training machine learning models that predict Root Mean Square Fluctuation (RMSF) from protein structure information.

The resulting datasets enable the development of ML architectures that can accurately predict protein dynamics from structural features.

## Features

- Extract RMSF data for each simulation and replica
- Process PDB data from H5 files with proper cleaning
- Calculate core/exterior classification for protein residues
- Generate comprehensive ML features including DSSP secondary structure
- Voxelize protein structures using aposteriori
- Create visualizations to analyze and validate the data

## Installation

### Prerequisites

- Python 3.7 or higher
- aposteriori (for voxelization)
- MSMS (optional, for improved core/exterior classification)

### Setup

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/mdcath-processing.git
   cd mdcath-processing
   ```

2. Create a virtual environment:
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install the required packages:
   ```
   pip install -r requirements.txt
   ```

4. Install aposteriori:
   ```
   pip install aposteriori
   ```

5. Install the package in development mode:
   ```
   pip install -e .
   ```

## Usage

### Basic Usage

Process all domains in the default location:

```
python main.py
```

### Configuration

You can customize the processing by creating your own configuration file:

```
python main.py -c path/to/your/config.yaml
```

See `src/mdcath/config/default_config.yaml` for available configuration options.

### Output Structure

The processed data will be organized in the following directory structure:

```
outputs/
├── RMSF/
│   ├── replicas/
│   │   ├── replica_0/
│   │   │   ├── 320/
│   │   │   ├── ...
│   │   ├── ...
│   └── replica_average/
│       ├── 320/
│       ├── ...
│       └── average/
├── pdbs/
├── frames/
├── voxelized/
├── ML_features/
└── visualizations/
```

## Processing Steps

1. **Data Loading**: Extract data from H5 files
2. **RMSF Processing**: Calculate and average RMSF values
3. **PDB Processing**: Clean and format PDB files
4. **Core/Exterior Classification**: Identify core and exterior residues
5. **ML Feature Generation**: Create features for machine learning
6. **Voxelization**: Create 3D grid representations for CNN models
7. **Visualization**: Generate plots and visualizations

## License

MIT License - See LICENSE file for details.

## Acknowledgments

- mdCATH dataset created by [Institution]
- MSMS software by [Creator]
- aposteriori voxelization tool by [Creator]
# mdcath-processor2
