# DL augmentation

DL augmentation is the package for neural network-assisted atomic electron tomography.
DL augmentation package has two parts; data generator codes (written by matlab) & deep learning codes (written by python, pytorch package).

Devoloper: Juhyeok Lee

For the detail, please see the paper: J. Lee, C. Jeong and Y. Yang, “Single-atom level determination of 3-dimensional surface atomic structure via neural network-assisted atomic electron tomography”, Nature Communications (2021).


## Data generator for DLaugmentation (matlab)
DL data generator consists of generating input/target data.

### ** Requirements **
- matlab (>= 2017a)
- GENFIRE (reconstruction algorithm)

To use this data generator, please download the GENFIRE package.
(https://www.physics.ucla.edu/research/imaging/dataSoftware.html)


### ** Test **
```
cd ./DL_data_generator
```
1. DLDGen_main.m for crystal structure. DLDGen_main_amorphous for amorphous structure.
2. Set the folder location, file name, and some parameters for generating atomic model.
3. run ./DL



## DL augmentation (python, pytorch package)


```
conda install test
```
