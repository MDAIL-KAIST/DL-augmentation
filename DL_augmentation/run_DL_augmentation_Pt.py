#!/usr/bin/env python
# coding: utf-8

# Author: J.Lee, KAIST (Korea), 2020.
# 
# Y.Yang, Multi-Dimensional Atomic Imaging Lab, KAIST
# 
# DL augmentation code
# 
# If you use the DL augmentation or related materials, we would appreciate citation to the following paper:
# J. Lee, C. Jeong and Y. Yang, “Single-atom level determination of 3-dimensional surface atomic structure via neural network-assisted atomic electron tomography”, Nature Communications (2021)


import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import numpy as np
import torch.nn.functional as F
import torch.optim as optim
import scipy.io
import time
from torch.utils import data

import DL_aug as DLa


### checking that cuda is available or not ###
USE_CUDA=torch.cuda.is_available()
DEVICE=torch.device("cuda" if USE_CUDA else "cpu")
print("CUDA: {}".format(USE_CUDA))
###

###########################################################
### input parameter ###
INPUT_PATH='./Pt_inputdata'

INPUT_FILE_NAME='Pt_input_1_real_intepolation_zero_padding'
INPUT_INSIDE_NAME='ESTvol'
data_size=144;

OUTPUT_PATH='./Pt_inputdata'
OUTPUT_FILE_NAME='./Pt_output_1'
OUTPUT_INSIDE_NAME='output'
###
###########################################################

### generate DL-augmenation model(Unet) ###
aut = DLa.UnetGenerator_3d(in_dim=1,out_dim=1,num_filter=12).to(DEVICE)
print("model contructing: OK!")
###



### loading a previous saved model parameter ###
PATH = './DL_aug_save_file/DL_aug_Pt_FCC_Bf5' # FCC Bf5
#PATH = './DL_aug_save_file/DL_aug_Pt_FCC_Bf3.2' # FCC Bf3.2
#PATH = './DL_aug_save_file/DL_aug_Pt_amorphous_Bf5' # amorphous Bf5
aut.load_state_dict(torch.load(PATH, map_location={'cuda:0': 'cpu'}))
#aut.load_state_dict(torch.load(PATH))
aut.to(DEVICE)

print("loading save file: OK!")
###


###
input_data = scipy.io.loadmat('{}/{}.mat'.format(INPUT_PATH,INPUT_FILE_NAME))[INPUT_INSIDE_NAME];
aut.eval()
with torch.no_grad():
    inputs = torch.tensor(input_data).view(-1,1,data_size,data_size,data_size).float().to(DEVICE);
    
    ### choose intensity scale factor
    inputs = inputs*13  # for Bf5  (FCC+Bf5, amorphous+Bf5)
    #inputs = inputs*8  # for Bf3.2 (FCC+Bf3.2)
    ###
    
    outputs = aut(inputs)
    outputs = outputs.data[0][0].cpu().numpy()
    scipy.io.savemat('{}/{}.mat'.format(OUTPUT_PATH,OUTPUT_FILE_NAME), {'{}'.format(OUTPUT_INSIDE_NAME):outputs}) # save
###




