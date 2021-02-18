#!/usr/bin/env python
# coding: utf-8

# Author: J.Lee, KAIST (Korea), 2020.
# Y.Yang, Multi-Dimensional Atomic Imaging Lab, KAIST
# DL augmentation
# If you use the DL augmentation or related materials, we would appreciate citation to the following paper:
# J. Lee, C. Jeong and Y. Yang, “Single-atom level determination of 3-dimensional surface atomic structure via neural network-assisted atomic electron tomography”, Nature Communications (2021)


import torch
import torch.nn as nn
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

#########################################################
### input parameter for training ###
N_of_data=10000;  # training set. number
N_of_vdata=1000;  # validation set, number
N_of_tdata=1000;  # test set, number

batch_size=1;  # batch_size
N_epoch=100;     # epoch
data_size=144;  # data_size (ex 144x144x144 volume -> 144)
###

### input & target folder path ###
INPUT_PATH='input data path';
TARGET_PATH='target data path';
###
#########################################################


### data loading setting ###
params1 = {'INPUT_FILE_NAME':'Pt_input', 'TARGET_FILE_NAME' : 'Pt_target','INPUT_INSIDE_NAME':'GF_Vol', 'TARGET_INSIDE_NAME':'target'}
params2 = {'batch_size': batch_size, 'shuffle': True, 'num_workers': 12}

N_of_start=1
train_ID_set = range(N_of_start,N_of_start+N_of_data);
validation_ID_set = range(N_of_start+N_of_data,N_of_start+N_of_data+N_of_vdata)
test_ID_set = range(N_of_start+N_of_data+N_of_vdata,N_of_start+N_of_data+N_of_vdata+N_of_tdata)


train_dataset = DLa.Dataset(train_ID_set,INPUT_PATH, TARGET_PATH, DLa.image_shift(4), **params1)
validation_dataset = DLa.Dataset(validation_ID_set,INPUT_PATH, TARGET_PATH, DLa.image_shift(4), **params1)
test_dataset = DLa.Dataset(test_ID_set,INPUT_PATH, TARGET_PATH, None, **params1)


train_generator = data.DataLoader(train_dataset, **params2)
validation_generator = data.DataLoader(validation_dataset, **params2)
test_generator = data.DataLoader(test_dataset, **params2)
###


### generate model 
aut = DLa.UnetGenerator_3d(in_dim=1,out_dim=1,num_filter=12).to(DEVICE)
print("model contructing: OK!")
###


### define loss function & optimizer ###
criterion = nn.MSELoss()
optimizer = optim.Adam(aut.parameters(), lr=0.0002, betas=(0.9, 0.999), eps=1e-08, weight_decay=0, amsgrad=False)
#optimizer = optim.SGD(aut.parameters(), lr=0.05, momentum=0.9)
###

 
### training process ###
print("Training starts")

epoch_print=5000//batch_size;
total_startTime = time.time()
for epoch in range(N_epoch):  # loop over the dataset multiple times
    aut.train()
    startTime = time.time()
    running_loss = 0.0
    for i, (inputs, target) in enumerate(train_generator):
        # input & target data
        inputs=inputs.view(-1,1,data_size,data_size,data_size).float().to(DEVICE);
        target=target.view(-1,1,data_size,data_size,data_size).float().to(DEVICE);
        
        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs = aut(inputs)
        loss = criterion(outputs, target)
        loss.backward()
        optimizer.step()

        # print train imformation
        running_loss += (loss.item())**0.5
        #print(running_loss)
        #if i % epoch_print== epoch_print-1:    # print every 1000 mini-batches
        ##print('[i: %d, %4d %%] loss: %.10f' %(i + 1, (i + 1)/N_of_data*batch_size*100, running_loss / epoch_print))
        #running_loss = 0.0

    endTime = time.time() - startTime


    # calculating loss of training set & validation set 
    aut.eval()
    with torch.no_grad():
        loss_sum_test=0
        for j, (inputs, target) in enumerate(train_generator):
            inputs=inputs.view(-1,1,data_size,data_size,data_size).float().to(DEVICE);
            target=target.view(-1,1,data_size,data_size,data_size).float().to(DEVICE);

            outputs = aut(inputs)
            loss_test = criterion(outputs, target)
            loss_sum_test += (loss_test.item())**0.5 # MSE -> RMSE
        print('[epoch: %d, %3d %%] training set loss: %.10f '  %(epoch + 1, (epoch + 1)/N_epoch*100 , loss_sum_test/(j+1)))

        loss_sum_test=0
        for j, (inputs, target) in enumerate(validation_generator):
            inputs=inputs.view(-1,1,data_size,data_size,data_size).float().to(DEVICE);
            target=target.view(-1,1,data_size,data_size,data_size).float().to(DEVICE);

            outputs = aut(inputs)
            loss_test = criterion(outputs, target)
            loss_sum_test += (loss_test.item())**0.5 # MSE -> RMSE
        print('[epoch: %d, %3d %%] validation set loss: %.10f time: %.3f '  %(epoch + 1, (epoch + 1)/N_epoch*100 , loss_sum_test/(j+1) ,endTime))

total_endTime = time.time() - total_startTime
print('Training has been finished')
print('Total time: %.3f'  %(total_endTime))
###


### calculate loss of test set ####
aut.eval()
with torch.no_grad():
    loss_sum_test=0
    for j, (inputs, target) in enumerate(test_generator):
        inputs=inputs.view(-1,1,data_size,data_size,data_size).float().to(DEVICE);
        target=target.view(-1,1,data_size,data_size,data_size).float().to(DEVICE);

        outputs = aut(inputs)
        loss_test = criterion(outputs, target)
        loss_sum_test += (loss_test.item())**0.5 # MSE -> RMSE
    print('test set loss: %.10f '  %(loss_sum_test/(j+1)))
###



### save the result ###
PATH = './DL_aug_save_file/'
#torch.save(aut.state_dict(), PATH)
#print('saving model: OK!')
###
