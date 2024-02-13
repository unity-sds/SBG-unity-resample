#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus
"""

import json
import shutil
import sys
import os
import hytools as ht
from hytools.io.envi import WriteENVI
import numpy as np
from scipy.interpolate import interp1d
from skimage.util import view_as_blocks
from PIL import Image

def main():
    ''' Perform a two-step spectral resampling to 10nm. First wavelengths are aggregated
    to approximateley 10nm, then aggregated spectra a interpolated to exactly 10nm using a
    piecewise interpolator.
    '''

    run_config_json = sys.argv[1]

    with open(run_config_json, 'r') as in_file:
        run_config =json.load(in_file)

    os.mkdir('output')

    print ("Resampling reflectance")

    rfl_base_name = os.path.basename(run_config['inputs']['reflectance_dataset'])
    sister,sensor,level,product,datetime,in_crid = rfl_base_name.split('_')

    crid = run_config['inputs']['crid']

    rfl_file = f'input/{rfl_base_name}/{rfl_base_name}.bin'
    rfl_met = rfl_file.replace('.bin','.met.json')

    out_rfl_file =  f'output/SISTER_{sensor}_L2A_RSRFL_{datetime}_{crid}.bin'
    out_rfl_met = out_rfl_file.replace('.bin','.met.json')

    resample(rfl_file,out_rfl_file)

    generate_metadata(rfl_met,out_rfl_met,
                      {'product': 'RSRFL',
                      'processing_level': 'L2A',
                      'description' : '10nm resampled reflectance'})
    generate_quicklook(out_rfl_file)

    print ("Resampling uncertainty")

    unc_base_name = os.path.basename(run_config['inputs']['uncertainty_dataset'])
    sister,sensor,level,product,datetime,in_crid,subproduct = unc_base_name.split('_')

    unc_file = f'input/{unc_base_name}/{unc_base_name}.bin'
    unc_met = unc_file.replace('.bin','.met.json')

    out_unc_file =  f'output/SISTER_{sensor}_L2A_RSRFL_{datetime}_{crid}_UNC.bin'
    out_unc_met = out_unc_file.replace('.bin','.met.json')

    resample(unc_file,out_unc_file)

    generate_metadata(unc_met,out_unc_met,
                      {'product': 'RSUNC',
                      'processing_level': 'L2A',
                      'description' : '10nm resampled uncertainty'})

    shutil.copyfile(run_config_json,
                    out_rfl_file.replace('.bin','.runconfig.json'))

    shutil.copyfile('run.log',
                    out_rfl_file.replace('.bin','.log'))


def generate_metadata(in_file,out_file,metadata):

    with open(in_file, 'r') as in_obj:
        in_met =json.load(in_obj)

    for key,value in metadata.items():
        in_met[key] = value

    with open(out_file, 'w') as out_obj:
        json.dump(in_met,out_obj,indent=3)

def gaussian(x,mu,fwhm):

    c = fwhm/(2* np.sqrt(2*np.log(2)))
    return np.exp(-1*((x-mu)**2/(2*c**2)))

def resample(in_file,out_file):

    image = ht.HyTools()
    image.read_file(in_file,'envi')

    if image.wavelengths.max()< 1100:
        new_waves = np.arange(400,991,10)
    else:
        new_waves = np.arange(400,2501,10)

    bins = int(np.round(10/np.diff(image.wavelengths).mean()))
    agg_waves  = np.nanmean(view_as_blocks(image.wavelengths[:(image.bands//bins) * bins],
                                           (bins,)),axis=1)

    if bins ==1 :
        agg_fwhm  = image.fwhm

    else:
        hi_res_waves = np.arange(300,2600)
        fwhm_array = np.zeros((image.bands,2600-300))

        for i,(wave,fwhm) in enumerate(zip(image.wavelengths,image.fwhm)):
            fwhm_array[i] = gaussian(hi_res_waves,wave,fwhm)

        sum_fwhm  = np.nansum(view_as_blocks(fwhm_array[:(image.bands//bins) * bins],
                                               (bins,2600-300)),axis=(1,2))
        agg_fwhm = []
        for i in range(len(agg_waves)):
            arg_max = np.argmax(sum_fwhm[i])
            half_max =  sum_fwhm[i].max()/2
            diff = np.abs(sum_fwhm[i]-half_max)
            end = arg_max + np.argmin(diff[arg_max:])
            start = np.argmin(diff[:arg_max])
            agg_fwhm.append(hi_res_waves[end] - hi_res_waves[start])

    #True resampled FWHM is difficult to determine, using a simple nearest neighbor approximation
    resampled_fwhm = interp1d(agg_waves,agg_fwhm,fill_value = 'extrapolate', kind = 'nearest')(new_waves)

    print(f"Aggregating every {bins} bands")

    out_header = image.get_header()
    out_header['bands'] = len(new_waves)
    out_header['wavelength'] = new_waves.tolist()
    out_header['fwhm'] = resampled_fwhm
    out_header['default bands'] = []

    if  "UNC" in in_file:
        out_header['description'] ='10 nm resampled reflectance uncertainty'
    else:
        out_header['description'] ='10 nm resampled reflectance'

    writer = WriteENVI(out_file,out_header)
    iterator =image.iterate(by = 'line')

    while not iterator.complete:
        line = iterator.read_next()[:,:(image.bands//bins) * bins]
        line  = np.nanmean(view_as_blocks(line,(1,bins,)),axis=(2,3))
        interpolator = interp1d(agg_waves,line,fill_value = 'extrapolate', kind = 'cubic')
        line = interpolator(new_waves)
        writer.write_line(line,iterator.current_line)

def generate_quicklook(input_file):

    img = ht.HyTools()
    img.read_file(input_file)
    image_file = input_file.replace('.bin','.png')

    if 'DESIS' in img.base_name:
        band3 = img.get_wave(560)
        band2 = img.get_wave(850)
        band1 = img.get_wave(660)
    else:
        band3 = img.get_wave(560)
        band2 = img.get_wave(850)
        band1 = img.get_wave(1660)

    rgb=  np.stack([band1,band2,band3])
    rgb[rgb == img.no_data] = np.nan

    rgb = np.moveaxis(rgb,0,-1).astype(float)
    bottom = np.nanpercentile(rgb,5,axis = (0,1))
    top = np.nanpercentile(rgb,95,axis = (0,1))
    rgb = np.clip(rgb,bottom,top)
    rgb = (rgb-np.nanmin(rgb,axis=(0,1)))/(np.nanmax(rgb,axis= (0,1))-np.nanmin(rgb,axis= (0,1)))
    rgb = (rgb*255).astype(np.uint8)

    im = Image.fromarray(rgb)
    im.save(image_file)

if __name__ == "__main__":
    main()
