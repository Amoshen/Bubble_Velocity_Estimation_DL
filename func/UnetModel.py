# -*- coding: utf-8 -*-
"""
Created on Feb 2018
@author: fangshuyang (yfs2016@hit.edu.cn)
"""

################################################
########        DESIGN   NETWORK        ########
################################################

import torch.nn as nn
import torch
import torch.nn as nn
import torch.nn.functional as F


class unetConv2(nn.Module):
    def __init__(self, in_size, out_size, is_batchnorm):
        super(unetConv2, self).__init__()
        # Kernel size: 3*3, Stride: 1, Padding: 1
        if is_batchnorm:
            self.conv1 = nn.Sequential(nn.Conv2d(in_size, out_size, 3, 1, 1),
                                       nn.BatchNorm2d(out_size),
                                       nn.ReLU(inplace=True),)
            self.conv2 = nn.Sequential(nn.Conv2d(out_size, out_size, 3, 1, 1),
                                       nn.BatchNorm2d(out_size),
                                       nn.ReLU(inplace=True),)
        else:
            self.conv1 = nn.Sequential(nn.Conv2d(in_size, out_size, 3, 1, 1),
                                       nn.ReLU(inplace=True),)
            self.conv2 = nn.Sequential(nn.Conv2d(out_size, out_size, 3, 1, 1),
                                       nn.ReLU(inplace=True),)
    def forward(self, inputs):
        outputs = self.conv1(inputs)
        outputs = self.conv2(outputs)
        return outputs


class unetDown(nn.Module):
    def __init__(self, in_size, out_size, is_batchnorm):
        super(unetDown, self).__init__()
        self.conv = unetConv2(in_size, out_size, is_batchnorm)
        self.down = nn.MaxPool2d(2, 2, ceil_mode=True)

    def forward(self, inputs):
        outputs = self.conv(inputs)
        outputs = self.down(outputs)
        return outputs




class  UnetModel(nn.Module):
    def __init__(self, n_classes, in_channels ,is_deconv, is_batchnorm, data_dim, label_dim):
        super(UnetModel, self).__init__()
        self.is_deconv     = is_deconv
        self.in_channels   = in_channels
        self.is_batchnorm  = is_batchnorm
        
        filters = [4, 8, 16, 32, 64]
        
        self.down1   = unetDown(self.in_channels, filters[0], self.is_batchnorm)
        self.down2   = unetDown(filters[0], filters[1], self.is_batchnorm)
        self.down3   = unetDown(filters[1], filters[2], self.is_batchnorm)
        self.down4   = unetDown(filters[2], filters[3], self.is_batchnorm)
        self.center  = unetConv2(filters[3], filters[4], self.is_batchnorm)
        self.mlp     = nn.Linear(189*filters[4],int(label_dim[1]))
        
    def forward(self, inputs,label_dsp_dim):
        down1  = self.down1(inputs)
        down2  = self.down2(down1)
        down3  = self.down3(down2)
        down4  = self.down4(down3)
        center = self.center(down4)
        mlp    = self.mlp(center.reshape(center.size(0),1,-1)).reshape(center.size(0),1,label_dsp_dim[0],label_dsp_dim[1])
        
        return mlp.contiguous()
    
    # Initialization of Parameters
    def  _initialize_weights(self):
          for m in self.modules():
            if isinstance(m, nn.Conv2d):
                n = m.kernel_size[0] * m.kernel_size[1] * m.out_channels
                m.weight.data.normal_(0, sqrt(2. / n))
                if m.bias is not None:
                    m.bias.data.zero_()
            elif isinstance(m, nn.BatchNorm2d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()
            elif isinstance(m,nn.ConvTranspose2d):
                n = m.kernel_size[0] * m.kernel_size[1] * m.out_channels
                m.weight.data.normal_(0, sqrt(2. / n))
                if m.bias is not None:
                    m.bias.data.zero_()
            elif isinstance(m, nn.Linear):
                m.weight.data.normal_(0, 1)
                m.bias.data.zero_()
