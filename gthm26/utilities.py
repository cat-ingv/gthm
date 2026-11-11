import numpy as np
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import pandas as pd
import ast
import os
import json
import glob
import matplotlib.colors as mcolors
import re

import sys
import time
import argparse
import pickle


def get_idx_poi(poi_id, pois_file):
    """
    recover idx lon lat from poi id

    poi_id: the poi name (string)
    idx_poi: poi index (integer) 
    lon_poi: poi longitude (float) 
    lat_poi: poi latitude (float) 
    
    """

    #### TODO add depth 08.05.26
    poi_id_file = pd.read_csv(pois_file,
                              sep = r"\s+",
                              header = None,
                              names = ['lon', 'lat', 'id', 'depth'])

    id_lst = poi_id_file['id'].values
    lon_poi_all = poi_id_file['lon'].values
    lat_poi_all = poi_id_file['lat'].values
    depth_poi_all=poi_id_file['depth'].values
    idx_poi = np.where(id_lst==poi_id)[0][0] # check if there a more elegant way to do it
    lon_poi = lon_poi_all[idx_poi]
    lat_poi = lat_poi_all[idx_poi]
    dep_poi= depth_poi_all[idx_poi]
    idx_poi = idx_poi + 1 # correction to be coherent with the counter in the csv file

    return idx_poi, lon_poi, lat_poi, dep_poi


def greens_functions_db(poi_idx, gf_path, subduction, n_timesteps):
    '''
    '''

    gf_subfaults_path = os.path.join(gf_path, subduction, 'GF2')
    gf_subfaults = [f for f in os.listdir(gf_subfaults_path)]# if not f.startswith('slip')]
    n_subfaults = len(gf_subfaults)

    print(f'Subduction subfaults folder is: {gf_subfaults_path}')
    print(f'N. of subfaults for {subduction} is: {n_subfaults}')

    # pattern = r"GF2/slip_taper_geo_(\d+)_hysea"
    # print(pattern)
    # ############# get this parameters from pkls     
    # subdirs = sorted(
    # [d for d in os.listdir(gf_path) if os.path.isdir(os.path.join(gf_path, d))],
    # key=lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else float('inf')
    # )
    # print(subdirs)

    big_eta = np.zeros((n_subfaults, n_timesteps+1), dtype=object)

    # for ii, subdir in enumerate(subdirs):
    for i in range(n_subfaults):

        isb = str(i+1).zfill(3)

        full_path = os.path.join(gf_path, gf_subfaults_path, 'slip_taper_geo_' + isb + '_hysea')
        
        if os.path.isdir(full_path):
            # str_tmp = f"{(m := re.search(pattern, full_path)).group(1)}_{m.group(2)}"
            str_tmp = subduction + '_' + isb
            big_eta[i, 0] = str_tmp
            file_path = os.path.join(full_path, 'out_ts_new_th.pkl')

            if os.path.isfile(file_path):
                 
                with open(file_path, 'rb') as handle:
                    sample_pickle = pickle.load(handle)

                valid_id = sample_pickle['valid_ids']
                mask = np.isin(poi_idx, valid_id)

                if mask == True:
                    crop_time_name = 'crop_time_' + str(poi_idx)
                    crop_time = sample_pickle[crop_time_name]
                    eta_crop = sample_pickle[poi_idx]
                    lec = len(eta_crop)
          
                    ##### offset 
                    if lec ==17281:
                      hoff=eta_crop[0]
                      eta_crop=eta_crop-hoff

                    big_eta[i,-lec:] = (eta_crop)

                  
    row_sums = np.sum(big_eta[:, 1:], axis=1)  # shape (17000,)
    mask = row_sums != 0
    big_eta_redux = big_eta[mask, :]
    n_redux, _ = big_eta_redux.shape
    n_removed = n_subfaults - n_redux
    # name_matricione='GF_db/gf_matrix_'+str(my_name_poi)+'.npy'
    # np.save(name_matricione,big_eta3)
    # print(f'N. of removed subfaults with no contribution: {n_removed}/{n_subfaults}')

    return big_eta_redux


def plot_hazard_curve(mih_thresholds, all_hc, poi, output_path):
    '''
    Plots hazard curves and uncertainty bounds directly from the all_hc dictionary.
    '''
    plt.figure(figsize=(10, 7))

    # Track labels so the legend only shows main groups once
    added_labels = set()

    for rate_name, hc in all_hc.items():
        if rate_name == "mih_thresholds":
            continue
        # 1. Determine group color and legend label
        if rate_name.endswith("_gl"):
            color = "b"  # Blue for GTM-wGL
            group_label = "GTM-wGL"
            is_main_line = (rate_name == "rate_annual_gl")
        else:
            color = "m"  # Magenta for GTM-woGL
            group_label = "GTM-woGL"
            is_main_line = (rate_name == "rate_annual")

        # 2. Assign styles (Solid for main, dashed/translucent for bounds)
        if is_main_line:
            linestyle = "-"
            linewidth = 5
            alpha = 1.0
            label = group_label if group_label not in added_labels else None
            if label:
                added_labels.add(group_label)
        else:
            linestyle = "--"
            linewidth = 3
            alpha = 0.5
            label = None  # Suppress individual bounds from legend

        # 3. Plot curve
        plt.plot(
            mih_thresholds,
            hc,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            alpha=alpha,
            label=label
        )

    # Axis and aesthetic formatting
    plt.xlabel('Tsunami Maximum Amplitude (m)', fontsize=14, labelpad=10)
    plt.ylabel('Mean Annual Rate of Exceedance', fontsize=14, labelpad=10)
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim(10e-6, 1)
    plt.grid(True, which="both", ls=":", alpha=0.6)
    plt.title(f'Hazard Curves - POI: {poi}', fontsize=14)
    plt.legend(fontsize=12, loc="best", frameon=True)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()

    # Save output plot
    filename = os.path.join(output_path, f'HC_{poi}.png')
    plt.savefig(filename, format='png', dpi=300, bbox_inches='tight')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.close()
    return

def plot_dist_max(mag_rup, h_max_list, poi, output_path):
    
    plt.figure(figsize=(10, 7))

    plt.scatter(mag_rup, h_max_list, c=h_max_list, cmap='viridis', edgecolors='k', s=50)
    plt.colorbar(label='Max Wave Height $H_{max}$ (m)')
    plt.xlabel('Magnitude ($M_w$)')
    plt.ylabel('Max Wave Height $H_{max}$ (m)')
    plt.title(f'Magnitude vs ($H_{max}$)')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()

    filename = os.path.join(output_path, f'Mw_vs_MaxHeight_{poi}.png')
    plt.savefig(filename, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    return 



def calculate_hazard_curve(mih_thresholds, h_max_lst, exc_rates, poi):
    '''
    '''
    
    h_max_lst = np.transpose(h_max_lst)
    mask_exceed = h_max_lst[:, None] >= mih_thresholds[None, :]
    weighted_exceed = mask_exceed * exc_rates[:, None]
    hc = np.sum(weighted_exceed, axis=0)
    

    return hc 



