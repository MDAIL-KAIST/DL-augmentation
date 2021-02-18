# DL augmentation package

DL augmentation package is implemented for neural network-assisted atomic electron tomography.
This package has two parts; data generator codes (written by matlab) & deep learning codes (written by python, pytorch package).

Devoloper: Juhyeok Lee

For the detail, please see the paper: J. Lee, C. Jeong and Y. Yang, “Single-atom level determination of 3-dimensional surface atomic structure via neural network-assisted atomic electron tomography”, Nature Communications (2021).


## Data generator for DL augmentation (matlab)
DL data generator consists of generating input/target data.

### Requirements
- matlab (>= 2017a)
- GENFIRE (reconstruction algorithm)

To use this data generator, please download the GENFIRE package.
(https://www.physics.ucla.edu/research/imaging/dataSoftware.html)


### How to start
```
cd ./DL_data_generator
```
1. Choose DLDGen_main.m for crystal structure or DLDGen_main_amorphous.m for amorphous structure.
2. Set a folder location, a file name, and some parameters for generating atomic model.
3. Run DLDGen_main.m or DLDGen_main_amorphous.m


## DL augmentation (python, pytorch package)
DL augmentation 

### Requirements
- pytorch (>=1.1) (we tested 1.1, 1.2, 1.3)
- numpy
- scipy

### Reconmended environment
- Linux
- CUDA (we tested 10.1)
- GPU memory (>= 4 GB for 144<sup>3</sup> voxel size, >= 8 GB for 256^3 voxel size) (if used)
- conda environment


### How to start

```
conda install test
```
