#!/usr/bin/env python
# coding: utf-8

import torch
import torch.nn as nn
import numpy as np
import torch.nn.functional as F
import torch.optim as optim
import scipy.io
import time
from torch.utils import data



### define image trainslation ###
class image_shift(object):

    def __init__(self, translation_range=0):
        #assert isinstance(translation_range, (float, tuple))
        self.translation_range = translation_range

    def __call__(self, input_data, target_data):
        inputs, target = input_data, target_data
        
        dx, dy, dz = np.random.randint(self.translation_range*2+1, size=3)-self.translation_range
        
        inputs = np.roll(inputs, dz, axis=0)
        inputs = np.roll(inputs, dy, axis=1)
        inputs = np.roll(inputs, dx, axis=2)
        
        target = np.roll(target, dz, axis=0)
        target = np.roll(target, dy, axis=1)
        target = np.roll(target, dx, axis=2)
        
        if dz>0:
            inputs[:,:,:dz] = 0
            target[:,:,:dz] = 0
        elif dz<0:
            inputs[:,:,dz:] = 0
            target[:,:,dz:] = 0
        if dy>0:
            inputs[:,:dy,:] = 0
            target[:,:dy,:] = 0
        elif dy<0:
            inputs[:,dy:,:] = 0
            target[:,dy:,:] = 0
        if dx>0:
            inputs[:dx,:, :] = 0
            target[:dx,:, :] = 0
        elif dx<0:
            inputs[dx:, :,:] = 0
            target[dx:, :,:] = 0
        return inputs, target
###



### For data loading ###
class Dataset(data.Dataset):
    'Characterizes a dataset for PyTorch'
    def __init__(self, list_IDs, INPUT_PATH,TARGET_PATH,Transform=None,INPUT_FILE_NAME='Pt_input',TARGET_FILE_NAME='Pt_target',INPUT_INSIDE_NAME='GF_Vol',TARGET_INSIDE_NAME='target'):
        'Initialization'
        self.list_IDs = list_IDs
        self.input_path = INPUT_PATH
        self.target_path = TARGET_PATH
        self.transform = Transform
        self.input_file_name = INPUT_FILE_NAME
        self.target_file_name = TARGET_FILE_NAME
        self.input_inside_name = INPUT_INSIDE_NAME
        self.target_inside_name = TARGET_INSIDE_NAME

    def __len__(self):
        'Denotes the total number of samples'
        return len(self.list_IDs)

    def __getitem__(self, index):
        # ID_set
        ID = self.list_IDs[index]
        'Generates one sample of data'
        # Load data and get target_data for matlab file
        input_data = scipy.io.loadmat('{}/{}_{}.mat'.format(self.input_path,self.input_file_name,ID))[self.input_inside_name];
        target_data = scipy.io.loadmat('{}/{}_{}.mat'.format(self.target_path,self.target_file_name,ID))[self.target_inside_name];

        if self.transform:
            input_data, target_data = self.transform(input_data, target_data)

        return input_data, target_data

###


### define DL_aug(3D_Unet) function ###
def Conv_1(in_dim,out_dim,act_fn,ks=3,st=2,pa=1): #conv+batch+act
    model = nn.Sequential(
        nn.Conv3d(in_dim,out_dim, kernel_size=ks, stride=st, padding=pa),
        nn.BatchNorm3d(out_dim),
        act_fn,
    )
    return model

def Conv_2(in_dim,out_dim,act_fn,ks=3,st=2,pa=1,drp=0.2): #conv+batch+dropout+act
    model = nn.Sequential(
        nn.Conv3d(in_dim,out_dim, kernel_size=ks, stride=st, padding=pa),
        nn.BatchNorm3d(out_dim),
        nn.Dropout3d(p=drp),
        act_fn,
    )
    return model

def Maxpool_1(ks=2, st=2, pa=0): #maxpooling
    pool = nn.MaxPool3d(kernel_size=ks, stride=st, padding=pa)
    return pool

def Conv_trans_1(in_dim,out_dim,act_fn,ks=3,st=2,pa=1,op=1): #trans_conv+batch+act
    model = nn.Sequential(
        nn.ConvTranspose3d(in_dim,out_dim, kernel_size=ks, stride=st, padding=pa,output_padding=op),
        nn.BatchNorm3d(out_dim),
        act_fn,
    )
    return model

def Conv_trans_2(in_dim,out_dim,act_fn,ks=3,st=2,pa=1,op=1,drp=0.2): #trans_conv+batch+dropout+act
    model = nn.Sequential(
        nn.ConvTranspose3d(in_dim,out_dim, kernel_size=ks, stride=st, padding=pa,output_padding=op),
        nn.BatchNorm3d(out_dim),
        nn.Dropout3d(p=drp),
        act_fn,
    )
    return model

def Conv_trans_3(in_dim,out_dim,act_fn,ks=3,st=2,pa=1,op=1): #trans_conv+act
    model = nn.Sequential(
        nn.ConvTranspose3d(in_dim,out_dim, kernel_size=ks, stride=st, padding=pa,output_padding=op),
        act_fn,
    )
    return model
###


### Generate DL-augmentation ###
class UnetGenerator_3d(nn.Module):

    def __init__(self,in_dim,out_dim,num_filter):
        super(UnetGenerator_3d,self).__init__()
        self.in_dim = in_dim
        self.out_dim = out_dim
        self.num_filter = num_filter
        act_fn1 = nn.LeakyReLU(0.2, inplace=True)
        act_fn2 = nn.ReLU()

        self.down_1 = Conv_2(self.in_dim,self.num_filter*2,act_fn1)
        self.pool_1 = Maxpool_1()
        self.down_2 = Conv_1(self.num_filter*2,self.num_filter*4,act_fn1)
        self.pool_2 = Maxpool_1()

        self.bridge = Conv_2(self.num_filter*4,self.num_filter*8,act_fn1,3,1,1,0.2)

        self.trans_up_11 = Conv_trans_1(self.num_filter*8,self.num_filter*8,act_fn1)
        self.up_1 = Conv_1(self.num_filter*12,self.num_filter*4,act_fn1,3,1,1)
        self.trans_up_12 = Conv_trans_1(self.num_filter*4,self.num_filter*4,act_fn1)
        self.trans_up_21 = Conv_trans_1(self.num_filter*4,self.num_filter*4,act_fn1)
        self.up_2 = Conv_2(self.num_filter*6,self.num_filter*2,act_fn1,3,1,1,0.2)
        self.out = Conv_trans_3(self.num_filter*2,self.out_dim,act_fn2)

    def forward(self,x):
        down_1 = self.down_1(x)
        pool_1 = self.pool_1(down_1)
        down_2 = self.down_2(pool_1)
        pool_2 = self.pool_2(down_2)

        bridge = self.bridge(pool_2)

        trans_up_11  = self.trans_up_11(bridge)
        concat_1 = torch.cat([trans_up_11,down_2],dim=1)
        up_1     = self.up_1(concat_1)
        trans_up_12  = self.trans_up_12(up_1)

        trans_up_21  = self.trans_up_21(trans_up_12)
        concat_2 = torch.cat([trans_up_21,down_1],dim=1)
        up_2     = self.up_2(concat_2)

        out = self.out(up_2)

        return out
    
###

