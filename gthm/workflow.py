import numpy as np
import netCDF4
import os
import sys
import re
import matplotlib.pyplot as plt
import pandas as pd
# import argparse
# import pygthm modules
#sys.path.append('src')
from . import utilities 

def greens_law(poi_dep, nomin_d=None):
    '''
    '''
    if nomin_d is None:
        nomin_d = 1.

    return (abs(poi_dep)/nomin_d)**0.25


def remap(elem, nrows, ncols):
    '''
    id in Davies to Remap
    '''
    elem0 = elem - 1    # elem starts from 0
    r = elem0 % nrows   # posizione lungo dip
    c = elem0 // nrows  # posizione lungo strike
    return r * ncols + c


def safe_float(val):
    '''
    conversion to float or None when not possible
    '''
    try:
        return float(val)
    except ValueError:
        return None


def read_mesh_file_dims(subduction):
    '''
    '''
    mesh_file = os.path.join('data', 'meshes/' + subduction + '_Ntot_Nstrike_Ndip.dat')
    n_tot, n_strike, n_dip = np.loadtxt(mesh_file, unpack=True, dtype='int')
    return n_strike, n_dip


def read_rupture_file(rupt_file):
    '''
    '''
    with netCDF4.Dataset(rupt_file) as rupture_file:
       
        # mag_rup = np.array(rupture_file.variables['Mw'])
        mag_rup = rupture_file.variables['Mw'][:]
        # mesh_index = np.array(rupture_file.variables['event_index_string'])
        mesh_index = rupture_file.variables['event_index_string'][:].filled(b'')
        # slip_string = np.array(rupture_file.variables['event_slip_string'])
        slip_string = rupture_file.variables['event_slip_string'][:].filled(b'')
        all_rate = {
            var: rupture_file.variables[var][:]
            for var in rupture_file.variables
            #if "rate_annual" in var and "variable" not in var
            if "rate_annual" in var
        }
    return mag_rup, mesh_index, slip_string, all_rate


def mesh_index_conversion(mesh_index, n_strike, n_dip):
    '''
    '''
    mesh_index_converted = []
    for row in mesh_index:
        # line = ''.join(c.decode('utf-8') for c in row).strip()
        line = b''.join(row).decode('utf-8').strip()

        remapped_row = [
            remap(int(x), n_dip, n_strike) + 1 
            for x in line.split('-') if x.strip()
        ]

        mesh_index_converted.append(remapped_row)

    return mesh_index_converted


def slip_conversion(slip_string):
    '''
    '''
    split_pattern = re.compile(r'[\s,_]+')
    slip = []
    for row in slip_string:
        # line = ''.join(c.decode('utf-8') for c in row).strip()
        line = b''.join(row).decode('utf-8').strip()
        # split using re by spaces, commas, or underscores without destroying '-'
        # raw_tokens = re.split(r'[\s,_]+', line)
        raw_tokens = split_pattern.split(line)
        slip_in_line = [
            # safe_float(x) for x in raw_tokens if safe_float(x) is not None
            val for x in raw_tokens if (val := safe_float(x)) is not None
        ]
        slip.append(slip_in_line)

    return slip    


def define_mih_thresholds(mih_thresholds_file=None):
    '''
    Threshold for hc
    '''
    if mih_thresholds_file is not None:
        mih_thresholds = np.load('mih_thresholds.npy', allow_pickle=True)
    else:
        ##TODO improve this part of mih definition
        #mih_thresholds = np.linspace(0.001,50)
        start = 0.02
        end = 100.
        step = (np.log(20.0) - np.log(start)) / 99
        # Generate log-spaced values from start up to log(100)
        log_array = np.arange(np.log(start), np.log(end) + step, step)
        mih_thresholds = np.exp(log_array)
        # print(mih_thresholds, len(mih_thresholds))
        # mih_thresholds = np.geomspace(0.02, 100.0, num=100)
    return mih_thresholds


def get_scenarios_mih(subduction, scenario_id, mesh_index_converted, slip, big_eta_redux, gl):
    '''
    '''
    str_eta = big_eta_redux[:,0]

    h_max_list = np.zeros(len(scenario_id))
    h_max_list_gl = np.zeros(len(scenario_id))

    for i in range(len(scenario_id)):
        ts_combined = None

        for ii, idx2 in enumerate(mesh_index_converted[i]):
            compl_name = f"{subduction}_{idx2:03d}"
            exst_condition = np.where(str_eta == compl_name)[0]
            if exst_condition.size == 0:
                continue

            # Correct element-wise slip retrieval and scaling
            #print(slip[i], i, ii, len(scenario_id), len(mesh_index_converted[i]))
            scaled_slip = slip[i][ii] / 10.0

            my_ts_ok = big_eta_redux[exst_condition[0], 1:]
            scaled_ts = my_ts_ok * scaled_slip

            if ts_combined is None:
                ts_combined = scaled_ts.copy()
            else:
                ts_combined += scaled_ts

        if ts_combined is not None:
            #h_max_list[i] = np.max(ts_combined * gl)
             h_max_list[i] = np.max(ts_combined)
             h_max_list_gl[i] = np.max(ts_combined * gl)

    return h_max_list, h_max_list_gl



def run(args):

    # variables from command line parser
    subduction = args.subd
    poi = args.poi
    subduction_id = args.subd_id
    rupt_file = args.rupt_file
    output_path = args.output_path
    gf_path = args.gf_path
    plot = args.plot

    # create output folder
    os.makedirs(output_path, exist_ok=True)

    # ==============================================================================
    # 1. Read input files (mesh, ruptures, pois)
    # ==============================================================================
    print('Reading data....')
    n_strike, n_dip = read_mesh_file_dims(subduction)
    pois_file = os.path.join('data', 'POI_lon_lat_ID_depth.txt')
    poi_idx, poi_lon, poi_lat, poi_dep = utilities.get_idx_poi(poi, pois_file)
    print(f'Selected POI: {poi}, {poi_idx}, {poi_lon}, {poi_lat}, {poi_dep}')
    mag_rup, mesh_index, slip_string, all_rate = read_rupture_file(rupt_file)

    # creation of scenarios id
    scenario_id = [f"{subduction_id}{str(i+1).zfill(7)}" for i in range(len(mesh_index))]

    # ==============================================================================
    # 2. Davies data conversion
    # ==============================================================================
    print('Conversion of strings to values....')
    mesh_index_converted = mesh_index_conversion(mesh_index, n_strike, n_dip)
    slip = slip_conversion(slip_string)

    # ==============================================================================
    # 3. GFs
    # ==============================================================================
    print('Generation of GF Matrix....')

    # calculate GL
    gl = greens_law(poi_dep)

    # mih thresholds
    mih_thresholds = define_mih_thresholds()

    n_timesteps = 17281

    big_eta_redux = utilities.greens_functions_db(poi_idx, gf_path, subduction, n_timesteps)

    h_max_list, h_max_list_gl = get_scenarios_mih(subduction, scenario_id, mesh_index_converted, slip, big_eta_redux, gl)

    # ==============================================================================
    # 4. Compute hazard curves using pre-loaded arrays
    # ==============================================================================
    print('Hazard curve computation....')

    all_hc = {'mih_thresholds': mih_thresholds}

    for rate_name, rates in all_rate.items(): 
        hc = utilities.calculate_hazard_curve(mih_thresholds, h_max_list, rates, poi)
        hc_gl = utilities.calculate_hazard_curve(mih_thresholds, h_max_list_gl, rates, poi)

        all_hc[rate_name] = hc
        all_hc[rate_name+ '_gl'] = hc_gl
        
       # print(f"Processed: {rate_name}")

    # ==============================================================================
    # 5. PLOTTING
    # ==============================================================================

    if plot:
        utilities.plot_hazard_curve(mih_thresholds, all_hc, poi, output_path)
        utilities.plot_dist_max(mag_rup, h_max_list, poi, output_path)


    # ==============================================================================
    # 6. Save hc to CSV
    # ==============================================================================
    print('Saving hazard curve....')

    df = pd.DataFrame(all_hc)
    csv_filename = f"all_hc_{subduction}_{poi}.csv"
    df.to_csv(os.path.join(output_path, csv_filename), index=False)

    # ==============================================================================
    # 7. Save Scenarios, H_max, and Rates to CSV
    # ==============================================================================
    print('Saving rates and hmax....')

    # Construct the output dictionary
    output_dict = {
        'scenario_id': scenario_id,
        'h_max_list': h_max_list,
        'h_max_list_gl': h_max_list_gl
    }

    # Add each annual rate array as a separate column
    for rate_name, rates in all_rate.items():
        # np.ravel ensures multi-dimensional rate arrays match length
        output_dict[rate_name] = np.ravel(rates)

    # Convert dictionary to DataFrame
    df_scenarios = pd.DataFrame(output_dict)

    # Define output file path & save
    csv_filename = f"scenarios_hmax_{subduction}_{poi}.csv"
    full_csv_path = os.path.join(output_path, csv_filename)

    df_scenarios.to_csv(full_csv_path, index=False)
    print(f"Saved scenario data successfully to {full_csv_path}")

