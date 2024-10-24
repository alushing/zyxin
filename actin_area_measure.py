#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import division

from glob import glob
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from skimage.color import label2rgb
import numpy as np
import os.path as op
import sh
import csv
import warnings
import pandas as pd
from pims import Bioformats
from PIL import Image
import skimage.filters as sf
import skimage.io as si
import skimage.measure as sme
import skimage.morphology as smo
import scipy.ndimage as sni
import pickle
from tifffile import imwrite

def compute_patches(folder_tif):
    
    actin_area_d = []
    
    for tif_filename in sorted(glob(op.join(folder_tif, '*.tif*'))):
        
        print('Loading ... ' + tif_filename)
        filename = (tif_filename.split('.')[0]).split('tiff')[1]
        
        name, raw_tyxc = load_TIF_series(tif_filename)
        
        print('Making masks...')
        actin_mask, actin_area = make_masks(raw_tyxc)
        
        actin_area_d.append(actin_area)
    
    
    output_txt = ('output/' + filename + '.txt')
    sh.mkdir('-p', op.dirname(op.abspath(output_txt)))
    with open(output_txt, 'w') as f:
        for d, actin_area in enumerate(actin_area_d):
            f.write('%r\n' % (actin_area))
                    
        
def load_TIF_series(filename):
    raw_tyxc = si.imread(filename)
    name = op.splitext(op.basename(filename))[0]
    return name, raw_tyxc

def make_masks(raw_tyxc, t_start=0, window_size=3):
    
    nt, ny, nx, nc = raw_tyxc.shape
    
    ###  generate actin mask
    raw_actin = raw_tyxc[t_start, :, :, 0]
    
    actin_thresh = raw_actin > sf.threshold_otsu(raw_actin)
    actin_mask = smo.remove_small_objects(actin_thresh, 4, 2)
    
    si.imsave('actin.tiff', actin_mask, check_contrast=False)
    
    actin_area = np.sum(actin_mask)
    
    return (actin_mask, actin_area)


def filter_props(folder_txt, lifetime_cutoff=0, ir_min=1, ir_max = 1000):
    
    filename = (folder_txt.split('/')[1])
    output_txt = ('output/filtered/' + filename + '.txt')    
    print(output_txt)

    lifetime_d = []
    intensity_ratio_max_d = []
    area_max_d = []
    aspect_ratio_max_d = []
    
    for file in sorted(glob(op.join(folder_txt, '*.txt'))):
        
        print('Loading: ' + file)
        
        with open(file, 'r') as f:
            lines=f.readlines()
            for x in lines:
                unclean_split = x.split(' ')
                clean_split = []
                for item in unclean_split:
                    if len(item) >= 1:
                        clean_split.append(float(item))
                
                lifetime_d.append(clean_split[0])
                intensity_ratio_max_d.append(clean_split[1])
                area_max_d.append(clean_split[2])
                aspect_ratio_max_d.append(clean_split[4])
            f.close()

    lifetime_D = []
    intensity_ratio_max_D = []
    area_max_D = []
    aspect_ratio_max_D = []
    for (lifetime, intensity_ratio_max, area_max,
         aspect_ratio_max) in zip(
             lifetime_d, intensity_ratio_max_d, area_max_d, aspect_ratio_max_d):
        if (lifetime > lifetime_cutoff) and (ir_max > intensity_ratio_max >
                                             ir_min):
            lifetime_D.append(lifetime)
            intensity_ratio_max_D.append(intensity_ratio_max)
            area_max_D.append(area_max)
            aspect_ratio_max_D.append(aspect_ratio_max)
            
    # Convert number of frames to actual lifetime.
    # Imaging interval: 1.17s
    lifetime_D = [lifetime * 1.17 for lifetime in lifetime_D]
    

    # Convert number of pixels to micron2. Pixel size: 0.267 um.
    area_max_D = [area_max * 0.0256 for area_max in area_max_D]
    
    sh.mkdir('-p', op.dirname(op.abspath(output_txt)))
    with open(output_txt, 'w') as f:
        for D, lifetime in enumerate(lifetime_D):
            f.write('%r  %.4f  %.4f  %.4f\n' % (
                lifetime, intensity_ratio_max_D[D], area_max_D[D],
                aspect_ratio_max_D[D]))


            
if __name__ == '__main__':
    for folder_tif in sorted(glob('tiff/*')):
        compute_patches(folder_tif)
    
    #for folder_txt in sorted(glob('output-f*/*')):
        #filter_props(folder_txt)

