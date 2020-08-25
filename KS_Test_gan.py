# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 08:53:08 2020

@author: 赵晨宇
"""

import ROOT as rt
import numpy as np
import h5py 
import sys 
import gc
import math
import argparse
rt.gROOT.SetBatch(rt.kTRUE)
from sklearn.utils import shuffle



def get_parser():
    parser = argparse.ArgumentParser(
        description='root to hdf5',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', action='store', type=str,
                        help='input root file')
    parser.add_argument('--output', action='store', type=str,
                        help='output root file')
    parser.add_argument('--tag', action='store', type=str,
                        help='tag name for plots')
    parser.add_argument('--str_particle', action='store', type=str,
                        help='e^{-}')


    return parser




if __name__ == '__main__':
     
    # load data
    file_in = '/hpcfs/bes/mlgpu/liaoyp/one_dim_CGAN/Pred/Pred_em_gan.h5'
    hf = h5py.File(file_in, 'r')
    df_dict = {}
    for key in hf.keys():
        #print('key=',key)
        df_dict[key] = hf[key][:]
    hf.close()
    
    #KS test
    h_real_dedx = rt.TH1F('s_H_real_dedx', '', 10000, -3, 3)
    h_fake_dedx = rt.TH1F('s_H_fake_dedx', '', 10000, -3, 3)
    
    df_dict['Pred_epoch99'][:,0] = df_dict['Pred_epoch99'][:,0]*2
    df_dict['Pred_epoch99'][:,1] = np.arccos(df_dict['Pred_epoch99'][:,1])*180/math.pi
        
    for i in range(df_dict['Pred_epoch99'].shape[0]):
        h_real_dedx.Fill(df_dict['Pred_epoch99'][i,2])
        h_fake_dedx.Fill(df_dict['Pred_epoch99'][i,3])
        
    Prob = h_real_dedx.KolmogorovTest(h_fake_dedx, '')
    print(Prob)
    
    if Prob > 0: print('Done!')
    print(Prob)
    

    
    
    
    
    
    
    
    