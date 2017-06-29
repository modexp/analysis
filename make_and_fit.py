#!/project/datagrid/anaconda/envs/joran/bin/python
# -*- coding: utf-8 -*-

############################################
# written by
# Joran Angevaare 4-5-2017
# jorang@xs4all.nl
# This program opens real data and makes MC datasets based on the real data. The fit results are saved. 
# The likelihoods of a modulation hypothesis and the no modulation hypothesis are compared. All this is
# interpreted by the logL.py
# See http://www.physics.purdue.edu/darkmatters/doku.php?id=modulation:an:binnedlikelihood for an explanation
############################################

############################################
# Definitions used in the program
############################################

data_dir = "/data/modulation/Raw_Data/combined/"
# anaf_dir = "/dcache/xenon/jorana/Modulation/processed_use_thesis/analysis_used/"
# cali_dir = "/dcache/xenon/jorana/Modulation/processed_use_thesis/calibration_used/"
anaf_dir = "/dcache/xenon/jorana/Modulation/processed/analysis/"
cali_dir = "/dcache/xenon/jorana/Modulation/processed/calibration/"
stbc_dir = "/user/jorana/Modulation/stoomboot/"
out_dir  = "/data/xenon/mod_admin/data/logL/"  # logL_noTcor_10jun logL29May

PERIOD = 1 # year. This gives which modulation periodicity is fitted

# For the temperature correction we use we need to get rid of wrong readings of the temperature 
# sensor, for our setup it is impossible to have temperatures above 40 cel (safety switch) or below
# 20 cel (due to the AC). If you want to run without the temperature correction use TCOR = False.
THIGH = 40 
TLOW  = 20
TCOR  = True
############################################

############################################
# Very important parameter that takes care of what is being opened:
fromdate = '201511' 
############################################


import ROOT                         # To open root files
import numpy as np                  # Fast
import matplotlib.cm as cm          # colorbars 
import matplotlib                   # Plotting
import datetime                     # For dates on axis 
from os import listdir              # To see what is in a folder
import os
import emcee                        # fitting
from iminuit import Minuit, describe, Struct
from matplotlib.ticker import MaxNLocator
matplotlib.use('Agg')
import matplotlib.pyplot as plt     # Plotting
matplotlib.rc('font', size = 24)                   # Use big fonts...
import argparse


############################################
# Parse the arguments
############################################

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('-ch','--chan', help='channel', required=True)
parser.add_argument('-pk','--peak', help='peak', required=True)
parser.add_argument('-tau','--tau', help='fake halflife', required=True)
parser.add_argument('-amp','--amp', help='fake modulation amplitude', required=True)
parser.add_argument('-phi','--phi', help='fake phi', required=True)
parser.add_argument('-num','--num_fake', help='number of fake datastets', required=True)
parser.add_argument('-fixa','--fixa', help='Fit for this fixed modulation amplitude', required=False)

args = vars(parser.parse_args())

CH  = int(args['chan'])
PK  = int(args['peak'])
TAU = float(args['tau'])
AMP = float(args['amp'])
PHI = float(args['phi'])
NUM = int(args['num_fake'])


if args['fixa'] != -1:
    FIX_A = float(args['fixa'] )
    if FIX_A != -1. : FIT_NDIM = 3
    else: FIT_NDIM = 4
else: FIT_NDIM = 4

print("make_and_fit::\tRunning for NDIM = %i"%FIT_NDIM)

############################################


max_fakes = 500              # There are only so many jobs a CPU can handle since the fake_data_lists get very large.
y2s     = 365.25 * 24 * 3600 # years to seconds
s2y     = 1. / y2s           # seconds to years
toffset = 70 * y2s           # the time stamp on the time from labview is 70 years off

# Energies of the sources
sourceE = [[1460], [1460],
           [511,   1157.020, 511 + 1157.020],
           [511,   1157.020, 511 + 1157.020],
           [1173.2,1332.5,   1173.2 + 1332.5],           
           [1173.2,1332.5,   1173.2 + 1332.5],
           [661.7], [661.7]]

# A list of the half-life of the source format half-life[channel] = tau1/2
halflife =  [1e9,   1e9, 59.1,   59.1, 5.2711, 5.2711, 30.08,  30.08]
dhalflife = [1e9,   1e9, 0.3,    0.3,  0.0004, 0.0004,  0.09,   0.09]

# Names of the sources
sourceName = ['Background', 'Background', 'Ti-44', 'Ti-44', 'Co-60', 'Co-60', 'Cs-137', 'Cs-137']
sourceLatex =['Background', 'Background', '$^{44}$Ti', '$^{44}$Ti', '$^{60}$Co', '$^{60}$Co', '$^{137}$Cs', '$^{137}$Cs']

col_lis = ['green', 'red', 'blue', 'cyan', 'brown', 'pink', 'magenta','black']         

# Some datafiles have "problems" these should not be opened
physics_exclude = [
     # Rates in all channel too small. Might be a problem in delta_t or some downtime
    'ANA_mx_n_20160817_1855.root', 
    # Problem with channel 7. Calibration constants changing too rapidly and we get strange rates henceforth
    'ANA_mx_n_20160530_0926.root', 'ANA_mx_n_20160601_0811.root',
    # Problem with Ti-channels, the high end of the energy spectrum is off
    'ANA_mx_n_20151130_0922.root', 'ANA_mx_n_20151202_0906.root', 'ANA_mx_n_20151204_1911.root',
    'ANA_mx_n_20151207_0819.root',
    # Double rates???
    'ANA_mx_n_20161028_0726.root',
    # The HV measurements
    'ANA_mx_n_20170220_0945.root', 'ANA_mx_n_20170220_1130.root', 'ANA_mx_n_20170220_1246.root',
    'ANA_mx_n_20170220_1444.root', 'ANA_mx_n_20170220_1900.root', 'ANA_mx_n_20170220_2027.root',
    'ANA_mx_n_20170221_0714.root', 'ANA_mx_n_20170221_0848.root', 'ANA_mx_n_20170221_1101.root',
    'ANA_mx_n_20170221_1242.root', 'ANA_mx_n_20170221_1426.root', 'ANA_mx_n_20170221_1607.root',
    'CAL_mx_n_20170220_0945.root', 'CAL_mx_n_20170220_1130.root', 'CAL_mx_n_20170220_1246.root',
    'CAL_mx_n_20170220_1444.root', 'CAL_mx_n_20170220_1900.root', 'CAL_mx_n_20170220_2027.root',
    'CAL_mx_n_20170221_0714.root', 'CAL_mx_n_20170221_0848.root', 'CAL_mx_n_20170221_1101.root',
    'CAL_mx_n_20170221_1242.root', 'CAL_mx_n_20170221_1426.root', 'CAL_mx_n_20170221_1607.root',
    'ANA_mx_n_20170307_1346.root', 'ANA_mx_n_20170307_1550.root', 'ANA_mx_n_20170307_1551.root', 
    'ANA_mx_n_20170307_1753.root', 'ANA_mx_n_20170307_2037.root', 'ANA_mx_n_20170308_0709.root',
    'ANA_mx_n_20170308_0930.root', 'ANA_mx_n_20170308_1152.root', 'ANA_mx_n_20170308_1401.root',
    'ANA_mx_n_20170308_1600.root', 'ANA_mx_n_20170308_1759.root',
    'CAL_mx_n_20170307_1346.root', 'CAL_mx_n_20170307_1550.root', 'CAL_mx_n_20170307_1551.root',
    'CAL_mx_n_20170307_1753.root', 'CAL_mx_n_20170307_2037.root', 'CAL_mx_n_20170308_0709.root', 
    'CAL_mx_n_20170308_0930.root', 'CAL_mx_n_20170308_1152.root', 'CAL_mx_n_20170308_1401.root', 
    'CAL_mx_n_20170308_1600.root', 'CAL_mx_n_20170308_1759.root',
    # disk space full, problems with writing files to disk
    'ANA_mx_n_20170413_0720.root', 'CAL_mx_n_20170413_0720.root'
    'ANA_mx_n_20170421_0747.root']

# The files that did not work, trying to open these files will result in a problem later on. There is 
badlist = [#'ANA_mx_n_20151016_1419.root', 
           'ANA_mx_n_20160412_0739.root', 'ANA_mx_n_20160418_0950.root', 'ANA_mx_n_20160527_1421.root', 
           'CAL_mx_n_20151016_1419.root', 'ANA_mx_n_20151016_1419.root', 'CAL_mx_n_20151022_1134.root', 
           'CAL_mx_n_20151022_1143.root', 'CAL_mx_n_20151022_1200.root', 'CAL_mx_n_20151022_1207.root', 
           'CAL_mx_n_20160412_0739.root', 'CAL_mx_n_20160418_0950.root', 'CAL_mx_n_20160527_1421.root',
           'CAL_mx_n_20161215_1047.root', 'ANA_mx_n_20161215_1047.root', 'ANA_mx_n_20170105_0836.root', 
           'CAL_mx_n_20170105_0836.root',
          # We had a power failure therefore exlude (there is also no data so it is not interesting:
          'CAL_mx_n_20170117_1607.root', 'CAL_mx_n_20170117_1620.root', 'CAL_mx_n_20170117_1648.root', 
          'CAL_mx_n_20170117_1653.root', 'CAL_mx_n_20170117_1942.root', 'CAL_mx_n_20170118_1634.root',
          # An I/O error caused LabView to crash
          'CAL_mx_n_20170127_1007.root', 'ANA_mx_n_20170127_1007.root'
          # A startbit was not in the correct place
          'ANA_mx_n_20170508_0709.root', 'CAL_mx_n_20170508_0709.root', 
          # 'ANA_mx_n_20170508_0709.root', 
          'CAL_mx_n_20170509_0743.root',
          'ANA_mx_n_20170510_1503.root',
          'ANA_mx_n_20150701_1303.root',
          'ANA_mx_n_20150702_1150.root',
          'ANA_mx_n_20150702_1349.root',
          'ANA_mx_n_20150707_1248.root',
          'ANA_mx_n_20150708_0755.root',
          'ANA_mx_n_20150708_0947.root',
          'ANA_mx_n_20150709_0727.root',
          'ANA_mx_n_20160617_0958.root', # All rates get to small in this file
          'ANA_mx_n_20160616_0951.root'  # All rates get to small in this file
          'ANA_mx_n_20161122_0758.root',
          'ANA_mx_n_20161121_1433.root',
          'ANA_mx_n_20161123_1218.root',]

badlist += physics_exclude 

############################################


###################### THE INITIALIZER ##########################
# This parts opens the analyzed files of the 
#################################################################
def files_to_open(set_date, i_want='ANA'):
    '''Opens files in this directory. I have CAL files and ANA files. I start from the first date as given
    by the first argument'''
    fnames,calnames, start, to_open = [], [], False, i_want + '_mx_n'
    
    # Loop through all the files
    for file in all_files:
        if i_want in file and set_date in file:
            start = True # Start when we want to start
        if any(file in b for b in badlist):
            continue
        if to_open in file and start == True:
            if i_want == 'ANA' : fnames.append(anaf_dir + file)                
            if i_want == 'CAL' : fnames.append(cali_dir + file)                
    return fnames

def file_opener(fnames, calnames):
    badfiles = []
    for k in in_ana: # For all properties from in_ana make a list of this property
        globals()[k] = []
    # For each file in fnames open it and append the property to the right list
    for i in range(len(fnames)):
        f= ROOT.TFile.Open(fnames[i]) # Open the ANA file
        try: # Sometimes the anafiles are corrupted, there will be no ana;1. 
            tree = f.Get("ana;1")
        except ReferenceError: # If this happens, a reference error occurs, add this name to a list of bad files
            badfiles.append(fnames[i]) 
            continue
        # A file can also be corrupted such that it is a TObject instead of a tree with leaves
        if 'TObject' in str(type(tree)): 
            badfiles.append(fnames[i])
            continue
        # For each event in the tree, get its properties and append to the right list
        for j , event in enumerate( tree) :
            for k in in_ana:
                if k == 'time': # Time is saved wrt t0, therefore need to be added to get the absolute time
                    if j == 0: dt = (getattr(event, k))
                    timestamp = getattr(event, k) + getattr(event, 't0')
                    globals()[k].append(timestamp)
                    time_start.append(timestamp - dt)
                    time_end.append(  timestamp + dt)
                elif k == 'frac': # for frac some events do not have this property, add a value of None here.
                    try: globals()[k].append( getattr(event, k) )
                    except AttributeError: globals()[k].append( None )
                elif k == 'bgrate' or k == 'bgdrate': # for chi2ndf some events do not have this property, add a value of None here.
                    try: globals()[k].append( getattr(event, k) )
                    except AttributeError: globals()[k].append( -1 )
                elif k == 'chi2ndf': # for chi2ndf some events do not have this property, add a value of None here.
                    try: globals()[k].append( getattr(event, k) )
                    except AttributeError: globals()[k].append( -1 )
                elif k == 'tot_error': # for tot_error some events do not have this property, add a value of None here.
                    try: globals()[k].append( getattr(event, k) )
                    except AttributeError: globals()[k].append( -1 )
                elif 'hv' in k: # for hv0-hv7 some events do not have this property, add a value of None here.
                        try: globals()[k].append( getattr(event, k) )
                        except AttributeError: globals()[k].append( None )    
                elif 'drate' in k:
                    if getattr(event, k) == 'nan' or getattr(event, k) == 'NaN':
                        print('Initializer::file_opener::A "None" occured')
                        globals()[k].append( -10 ) 
                    else: globals()[k].append( getattr(event, k) )       
                else: # Most poperties get added here to the right list
                    globals()[k].append( getattr(event, k) )
                
                    
    for k in in_cal: # For all properties from in_cal make a list of this property
        globals()[k] = []
    
    # For each file in calnames open it and append the property to the right list
    for i in range(len(calnames)):
        f = ROOT.TFile.Open(calnames[i]) # open the CAL file
        try: # Sometimes the calibration files are corrupted, there will be no cal;1. 
            tree = f.Get("cal;1")
        except ReferenceError: #If this happens, a reference error occurs, add this name to a list of bad files
                badfiles.append(calnames[i])
                continue
        # For each event in the tree, get its properties and append to the right list
        for j , event in enumerate(tree) :
            for k in in_cal:
                # Some of the properties are not a simple leaf in the (cal) root file but are root vectors
                # these are especially tricky to open. These can be opened by treating the rootvector as a 
                # variable and than converting it to a list before adding it to the list of this property
                if k == 'c0' or k == 'c1' or k == 'c2' or k =='chi2':
                    vec1 = getattr(event, k)
                    list1 = [value for value in vec1] # This 
                    globals()[k].append(list1)
                else: 
                    # Sometimes strange things happened with the calibration time therefore I check if the 
                    # timestamp does not have strange values.
                    if k == 'cal_tmin':
                        if getattr(event, k) < 1:
                            print('cal time is negative')
                            badfiles.append(calnames[i])
                        elif getattr(event, k) > 4e9: 
                            print('cal time to big in file %s' %calnames[i])
                            if calnames[i] not in badfiles:
                                badfiles.append(calnames[i])
                        else: globals()[k].append( getattr(event, k) )
    # As written below, we now have some properties from the calibration that are in a vector. To convert these
    # to a list we first create some lists:
    for string in ['c0','c1', 'c2', 'chi2']:
        for j in range(8):
            name = string + '_chan_' + str(j)
            globals()[name] = []
            for i in range(len(c0)):
                ccc = globals()[string]
                globals()[name].append(ccc[i][j])

    # And then convert these to np arrays.   
    for cx in ['c0','c1', 'c2', 'chi2']:
        for ch in range(8):
            name = cx + '_chan_' + str(ch)
            in_cal.append(name)

    for k in in_cal:
        if k == 'id':
            id_array = np.array(globals()[k])
        else:
            globals()[k] = np.array(globals()[k])
    in_ana.append('time_start')
    in_ana.append('time_end')
    for k in in_ana:
        globals()[k] = np.array(globals()[k])

    if not badfiles == []: print ('Initializer::The files that have a problem are:\t', badfiles)
    return badfiles

print("make_and_fit::\tLoad data")
all_files = sorted(listdir(anaf_dir)) # All files in this directory
all_files+= sorted(listdir(cali_dir)) # All files in this directory

fnames   = files_to_open(fromdate)
calnames = files_to_open(fromdate,'CAL')

print "make_and_fit::\tInitializer::I've opened the files:",fnames[0].split('/')[-1],', ... ,',fnames[-1].split('/')[-1], '(total:',    len(fnames),') and', calnames[0].split('/')[-1],', ... ,',calnames[-1].split('/')[-1], '(total:',len(calnames),')' 

in_ana = ['t0', 'time', 'channel', 'peak', 'rate', 'drate', 'e', 'res', 'temp', 'pres', 'bx', 'by',
          'bz', 'btot', 'humid','frac','hv0', 'hv1','hv2','hv3','hv4', 'hv5','hv6','hv7','chi2ndf', 
          'bgrate', 'bgdrate', 'tot_error']
in_cal = ['id', 'cal_tmin', 'cal_tmax', 'c0', 'c1', 'c2', 'chi2' ]         

time_start = []
time_end   = []
print("make_and_fit::\topen files")
badfiles = file_opener(fnames, calnames[-2:-1])

############################################
# Let's make sure we can actually store in the folders that we want to save figures in.
if not os.path.exists(out_dir+"/errors/"): 
    cmd = "mkdir " + out_dir+"/errors/"
    print(cmd)
    os.system(cmd)

figfolder = out_dir+"/fits/"
if not os.path.exists(figfolder): 
    cmd = "mkdir " + figfolder
    print(cmd)
    os.system(cmd)
for string in ['2','4','6']:
    if not os.path.exists(figfolder+string+"/"): 
        cmd = "mkdir " + figfolder +string+"/"
        print(cmd)
        os.system(cmd)

############################################


############################################ Fitter ############################################
# The fitter works as follows:
# 0. load the data (above)
# 1. make a fake data set based on the real data
# 2. fit that dataset
# 3. write the fit results to a txt file that is separately processed by logL.py
# 4. make a fit of the real data and save that in the same file

def cutter(t_start = -1, t_end = -1, ch_0 = -1, ch_1 = -1, pk_0 = -1, 
    pk_1 = -1, pk_2 = -1, T_low = -1, T_high = -1, P_min = -1, P_max = -1,
    chi_max = -1):
    '''Feed cuts to cut a ANA file and will return the right list. Note that we don't want -1's'''
    cuts = [t0>0]
    if t_start != -1: cuts = np.all([cuts, time > t_start], axis = 0)
    if t_end   != -1: cuts = np.all([cuts, time < t_end  ], axis = 0)
    if ch_1 < ch_0: ch_1, ch_0 = ch_0, ch_1
    if ch_1 != -1:
        if ch_0 != -1: cuts = np.all([cuts, np.any([channel == ch_0, channel == ch_1],axis = 0)], axis = 0)
        else: cuts = np.all([cuts, channel == ch_1], axis = 0)
    if pk_0   != -1: cuts = np.all([cuts, peak == pk_0], axis = 0)
    if pk_1   != -1: cuts = np.all([cuts, peak == pk_1], axis = 0)
    if pk_2   != -1: cuts = np.all([cuts, peak == pk_2], axis = 0)
    if T_low  != -1: cuts = np.all([cuts, T_low < temp], axis = 0)
    if T_high != -1: cuts = np.all([cuts, T_high > temp],axis = 0)
    if P_min  != -1: cuts = np.all([cuts, P_min < pres], axis = 0)
    if P_max  != -1: cuts = np.all([cuts, P_max > pres], axis = 0)
    if chi_max!= -1: cuts = np.all([cuts, chi2ndf < chi_max], axis = 0)
    return cuts

def binned_rate_xy(chanlist, 
        pk,    # Specify source (or channel) and pk
        plot_x,
        plot_y,
        plot_dy,
        fit_y,
        fit_res,               # format: [[A,da], [tau, dtau], [phi, dphi], [a, da]]
        fitted_pars = 2,       # 2-> A & t12 3 +> phi 4+>a
        sub = False,           # Plot the residuals of the rate - exponential fit
        save = False,          # Save the figure
        savename = False,      # Specify the name of the figure
        doplot = True,         # Show the plot
        incr = 1,              # The number of fits per bin
        jenkins = False,       # Only works if sub = True, plots the jenkins hypothesis in the residual
        bgrateplot = False,    # Very specific option, used to see the fitted BG rate in this source
        exp_text = False,      # Show the results of the exponential fit
        lit_val = False,       # Show the literature value of the halflife
        legend = True,         # If true, plot the legend
        coloritem = False,      # An ana_file propperty can be used as a colorbar
        colorbar_name = 'Temperature $[^\circ C]$',# The label that goes with the colorbar
        extra_text = False     # If some extra note needs to be placed in the figure, put it here.
        ):
    '''Plots the rate of a given source and peak number using the provided rate over time. It will return a figure of the 
    rate of this source over time. '''
    plt.rcParams['figure.figsize'] = (12.0, 12.0)    # Use big plots
    # Check the source name, if a channel is given, interpreted the corresponding source
    chan = chanlist[0]
    print('Rateplotter::\tbinned_rate_xy::plot')
    if type(coloritem) == np.ndarray: 
        std_col_bar  = np.std(coloritem)
        mean_col_bar = np.mean(coloritem)        
        col_bar      = np.clip(coloritem, mean_col_bar - 2 * std_col_bar, mean_col_bar + 2 * std_col_bar)
        numcut       = sum(np.any([coloritem < mean_col_bar - 2 * std_col_bar, coloritem > mean_col_bar + 2 * std_col_bar], axis = 0))
        print('Rateplotter::\tbinned_rate::\tClipped %i items for the colorbar = %.0e of total'%(numcut, numcut/len(col_bar)))          
    
    
    xs, ys, dys, xlist, cdata, plot_fit = [], [], [], plot_x, [], []
    
    # Effectively re-bin with a number of 'incr' measured rates per bin
    for i in range(0, len(xlist), incr):
        xmax = min(len(xlist), i + incr) 
        xs.append(np.average(plot_x[i:xmax]))          
        ys.append(np.average(plot_y[i:xmax]))
        plot_fit.append(np.average(fit_y[i:xmax]))
        if type(coloritem) == np.ndarray: cdata.append(np.average(col_bar[i:xmax]))
        # now the errors
        plot_dysum = 0
        for sig in plot_dy[i:xmax]: plot_dysum += (sig ** 2) 
        plot_dysum  = np.sqrt(plot_dysum) / len(plot_dy[i:xmax])
        dys.append(plot_dysum)
        
    plot_x, plot_y, plot_dy, plot_fit = np.array(xs), np.array(ys), np.array(dys), np.array(plot_fit)
    xsyears = (plot_x - toffset) * s2y

    # Open the fit results
    fit_res = np.transpose(fit_res)
    amp_fit, tau_fit, phi_fit, a_fit = fit_res[0]
    damp_fit, dtau_fit, dphi_fit, da_fit = (fit_res[1] + fit_res[2]) / 2
    fit_res = np.transpose(fit_res)
                
    # The plotting part
    plot_x_date = [datetime.datetime.fromtimestamp(np.int(x - toffset)) for x in plot_x]
     
    fig, ax1 = plt.subplots() 
    
    # Plot the rate subtracted from the fit and express in %'s
    if sub == True:    
        ax1.set_ylabel('$\Delta$ Rate [%]')
        ax1.errorbar(plot_x_date, plot_y, yerr = plot_dy, linestyle = 'None',  label = 'Data - exponential decay', color = 'black', capsize = 5, elinewidth = 2, zorder = 0 )
        
    # Plot the rates and fit it with
    elif sub == False:          
        ax1.errorbar(plot_x_date, plot_y, yerr = plot_dy, linestyle = 'None', marker = "o", markersize = 5, label = 'Data', color = 'b',  capsize = 5, elinewidth = 2, zorder = 0)
        ax1.set_ylabel('Rate [Hz]')
    if type(coloritem) == np.ndarray: plt.scatter(plot_x_date, plot_y, marker = "o", s = 75, c = cdata, cmap = plt.cm.rainbow, zorder = 3)    
        # The fit
    if not bgrateplot and not sub: ax1.plot(plot_x_date, plot_fit, color = 'red', linewidth = 2.0, linestyle = 'dashed', label = 'Fitted exponential decay')

    props = dict(boxstyle = 'round', facecolor = 'white', alpha = 1)
    if jenkins:
        if fitted_pars > 2 and sub: ax1.plot(plot_x_date, 100 * a_fit * np.cos(2 * np.pi * (xsyears-xsyears[0])/1 + phi_fit), label = 'MCMC sinus for T = 1 y', color='red', linewidth = 2, linestyle = '--')
        # Now compare this to jenkins findings (an amplitude of 0.034-0.015 % and a maximum at feb 1
        # +/- 14 day <- plot_dy or not taken into account here)
        jen_xs = np.linspace(int(xsyears[0]) - 1, xsyears[-1], num = 100)
        jen_xs -= int(xsyears[0]) - 1 # Shift so xs start at 0 i.e. at 1st of january
        jen_plot_y_min = sin_function(jen_xs, 0.034, (- 0.25 + 29 / 365.25) * 2 * np.pi ) # Small sin (note last arg. is phi)
        jen_plot_y_max = sin_function(jen_xs, 0.15 , (- 0.25 + 29 / 365.25) * 2 * np.pi ) # Large sin
        jen_xs += int(xsyears[0]) - 1 # Shift back
        # Cut anything that is calculated for a date smaller that the first (real) date
        jen_plot_y_min = jen_plot_y_min[jen_xs > xsyears[0]] 
        jen_plot_y_max = jen_plot_y_max[jen_xs > xsyears[0]]
        jen_xs     = jen_xs[jen_xs > xsyears[0]]    
        jen_exp_fit = np.power(2, -((jen_xs - jen_xs[0]) / tau_fit)) * amp_fit
        jen_xs_yr = jen_xs
        
        # Convert to nice dates
        jen_xs = [datetime.datetime.fromtimestamp(np.int(x * y2s)) for x in jen_xs]        
        if not sub: jen_plot_y_min, jen_plot_y_max = jen_exp_fit * (1 + 0.01 * jen_plot_y_min), jen_exp_fit * (1 + 0.01 * jen_plot_y_max )
        if fitted_pars > 2 and not sub: 
            ax1.plot(plot_x_date, 
                     np.power(2, -((xsyears - xsyears[0]) / tau_fit)) * amp_fit * 
                     (1 + a_fit * np.cos(2 * np.pi * xsyears/1 + phi_fit)), label = 'MCMC sinus for T = 1 y', color='red', linewidth = 2, linestyle = '--')
        # And plot
        ax1.plot(jen_xs, jen_plot_y_min, label = "O'Keefe et al. 2013, amplitude \n$0.034 - 0.15\%$ max. at 29 jan", color = 'g', linestyle = '-', linewidth = 2)
        ax1.plot(jen_xs, jen_plot_y_max, color = 'g', linestyle = '-', linewidth = 2)
        ax1.fill_between(jen_xs, jen_plot_y_min, jen_plot_y_max, color= 'g', alpha = 0.2) # nice region
        if sub: ax1.plot(jen_xs, jen_plot_y_min * 0, c = 'r', label = 'No modulation', linewidth = 2)
        
        # Add a little textbox about the fit results
        if fitted_pars > 2 :
            textstr = "fitted model: $\propto 1 + a \cdot\cos(2\pi t/T-\phi)$\n$a = %.2G \pm %.2G$ percent\n$\phi= %.2G \pm %.2G$\n$T=1$ yr" %(100 * a_fit, 100 * da_fit, phi_fit, dphi_fit)
            ax1.text(1.05, 0.1, textstr, color = 'r', transform = ax1.transAxes, fontsize = 24, verticalalignment='top', bbox=props)
            
        if legend: lgnd1 = ax1.legend( bbox_to_anchor = (1.04, 1), loc = 2, borderaxespad = 0., markerscale = 1, numpoints = 1)
    elif legend:   lgnd1 = ax1.legend(bbox_to_anchor = (1.04, 1), loc = 2, borderaxespad = 0., markerscale = 1, numpoints = 1)
    
    if type(coloritem) == np.ndarray:
        cb1 = plt.colorbar(ax=[ax1] , orientation = 'horizontal')
        cb1.set_label(colorbar_name)
    
    plt.xticks(rotation=45)
    if np.max(abs(plot_y) + abs(plot_dy)) < 0.15: plt.ylim(-0.15, 0.15)
    
    fig.autofmt_xdate()
          
    textstr = 'Source: %s at %.1f keV \n$t_{1/2}=%.5g\pm %.2g^{stat} $ year' % (sourceLatex[chan], sourceE[chan][pk], tau_fit, dtau_fit)
    if lit_val: textstr += '\nLiterature: $%.5g\pm%.2g$ year'%(halflife[chan], dhalflife[chan])
    if exp_text and chan >= 2: ax1.text(1.05, 0.4, textstr,  transform = ax1.transAxes, fontsize = 24, verticalalignment='top', bbox=props)
    if extra_text: ax1.text(0.05, 0.95, extra_text, transform = ax1.transAxes, fontsize = 24, verticalalignment='top', bbox=props)
    if save: fig.savefig(out_dir + savename +".png", dpi = 50, bbox_extra_artists=(lgnd1,), bbox_inches='tight')
    if doplot: plt.show()
    plt.clf()
    fig.clf()


def sin_function(x, A, phi):    
    '''returns A sin(2pi (x-phi)) with a period of 1 in the same units as x'''
    f = A * np.sin(2 * np.pi * (x) - phi / 1.)
    return f


def gaussian(x, mu, sigma):
    return np.exp((-(x-mu)**2)/(2 * sigma ** 2)) / (np.sqrt(2*np.pi) * sigma)


def lnlike(theta, x0, x1, y, yerr, tau12, dtau12, T, fixed_a):
    '''The likelihood function of the data. returns the total log-likelihood. This is calculated by
    evaluating the total probability. This can work for 1,2,3 or 4 parameters given by theta. I 
    mostly use 2 (overall rate and half-life) and 4 (overall rate, half-life, modulation phase and
    modulation amplitude.'''
    if len(theta) == 1:
        A = theta[0]
        tau, a, phi = tau12, 0, 0
    if len(theta) == 2: 
        [A, tau] = theta
        a, phi = 0, 0
    if len(theta) == 3: 
        [A, tau, phi] = theta
        if fixed_a: a = fixed_a
        else: a = 0
    if len(theta) == 4: [A, tau, phi, a] = theta

    # Terms without modulation
    model = A * tau * (np.power(2, - x0 / tau) - np.power(2, - x1 / tau)) / np.log(2)

    # Terms with modulation integrated over time between x0 and x1
    if len(theta) > 2:
        argx0 = 2 * np.pi * x0 / T + phi
        argx1 = 2 * np.pi * x1 / T + phi
        model+= (np.power(2, - (x0-x1) / tau) * a * A * T * tau * (
            np.power(2, x1 / tau) * ( T * np.cos(argx0) * np.log(2) -
            2 * np.pi * tau * np.sin(argx0)) + 
            np.power(2, x0 / tau) * (- T * np.cos(argx1) * np.log(2) +
            2 * np.pi * tau * np.sin(argx1)))/((2 * np.pi * tau) ** 2 + (T * np.log(2)) ** 2))

    inv_sigma2 = 1.0 / ((yerr * (x1-x0)) ** 2)


    return np.sum(-0.5 * (( y * (x1 - x0) - model) ** 2 * inv_sigma2) - np.log(np.sqrt(2 * np.pi * inv_sigma2))) 



def use_minuit(x0, x1, y, yerr, tau12, dtau12, chan, pk, ndim = 4, T = PERIOD, fixed_a= False,
        nwalkers = 200, nsteps = 300, burnin = 50, plot_level = [],
        save = False, extra_text = False,print_level = 0, ncall = 10000):
    '''Given a data with x_start, x_end (units years) and a rate (with error), an expected half-life tau 
    (and error) optimizes Rate0, and possible modulations (with amplitude a, phase phi and period T, 
    with a default period of 1 year. Specify the number of parameters with ndim (2 = exp, 3 = exp &
    phi, 4 = exp, phi & a))'''

    # Estimate the parameters
    A_est   = np.mean(y[0:10])
    A_max   = 10 * A_est
    tau_est = halflife[chan]
    phi_est = 0.0
    
    if fixed_a: a_est = fixed_a
    if ndim <= 2: a_est = 0
    if ndim == 4: a_est = 0.001
        
    est_list= [A_est, tau_est, phi_est, a_est]
    label_list = ['A', 'tau', 'phi', 'a']
    
    
    if ndim == 1: 
        def f(amp): return -lnlike([amp], x0, x1, y, yerr, tau12, dtau12, T, fixed_a)
        m = Minuit(f, 
                amp       = A_est,      limit_amp = (0, A_max),          error_amp = 1,  
                errordef  = 1.0, print_level = 0)      
    if ndim == 2: 
        def f(amp, t12): return -lnlike([amp,t12], x0, x1, y, yerr, tau12, dtau12, T, fixed_a)
        m = Minuit(f, 
            amp       = A_est,      limit_amp = (0, A_max),          error_amp = 1,  
            t12       = tau_est,    limit_t12 = (0, 1000),          error_t12 = dtau12, 
            errordef  = 1.0, print_level = 0)
    if ndim == 3: 
        def f(amp, t12, phi): return -lnlike([amp,t12, phi], x0, x1, y, yerr, tau12, dtau12, T, fixed_a)
        m = Minuit(f, 
                amp       = A_est,      limit_amp = (0, A_max),          error_amp = 5,  
                t12       = tau_est,    limit_t12 = (5, 80),            error_t12 = dtau12, 
                phi       = phi_est,    limit_phi = (-2 * np.pi, 2 * np.pi),    error_phi = 0.1,  
                errordef  = 1.0, print_level = 0)
    if ndim == 4: 
        def f(amp, t12, phi, a): return -lnlike([amp,t12, phi, a], x0, x1, y, yerr, tau12, dtau12, T, fixed_a)
        m = Minuit(f, 
                amp       = A_est,      limit_amp = (0, A_max),                 error_amp = 5,  
                t12       = tau_est,    limit_t12 = (5, 80),                    error_t12 = dtau12, 
                phi       = phi_est,    limit_phi = (-2 * np.pi, 2 * np.pi),    error_phi = 0.1,  
                a         = a_est,      limit_a   = (-0.1, 0.1),                error_a   = 0.0001,              
                errordef  = 1.0, print_level = print_level)

    m.migrad(ncall = ncall, resume = False)
    m.hesse()
    try: m.minos()
    except RuntimeError: 
        # If this happens, the fit had a hard time to find the optimum parameters in phi and a presumably because if a is 0 there is no change in 
        # likelihood for phi so every phi yields the same likelihood and makes it very hard to obtain accurate uncertainties
        print("\t\tDo miniuit.migrad to obtain the errors on the best fit")
    m.migrad()

    res_list = [[m.values[label], m.errors[label], m.errors[label]] for label in ['amp','t12', 'phi', 'a'][0:ndim]]
    for i in range(ndim + 1, 4 + 1):
        res_list.append([-1, 0, 0])    

    # Save some of the fits to be able to check the results
    if ndim >= 3 and save == True:
        # Open the fit results
        if ndim == 4: A, tau, phi, amp = [m.values[label] for label in ['amp','t12', 'phi', 'a'][0:ndim]]
        if ndim == 3: 
            amp = fixed_a 
            A, tau, phi = [m.values[label] for label in ['amp','t12', 'phi', 'a'][0:ndim]]
            res_list = [[m.values[label], m.errors[label], m.errors[label]] for label in ['amp','t12', 'phi', 'a'][0:ndim]]
            res_list.append([amp, 0 , 0])
        
        # For plotting the results
        xlist = (x1+x0)/2
        T = 1
        rate_exp = (A * np.power(2, - (xlist) / tau)) * (1 + amp * np.cos(2 * np.pi * xlist/T + phi) )
        
        # Plot the rate over time
        binned_rate_xy([chan, chan + 1], pk, xlist*y2s + time[0], y, yerr, rate_exp, res_list, fitted_pars = 4, save = True,
            savename = "/fits/%i/pk%i_t12_%i_a%.4f_phi%.2f_fix_%.4f_n_%i"%(CH, PK, TAU, AMP, PHI, FIX_A, FIT_NDIM) + datetime.datetime.today().isoformat(),
            doplot = False, incr = 100, exp_text = True, lit_val = False, legend = True, coloritem = False, extra_text = extra_text)
        
        # Plot the residuals over time
        rate_exp = (A * np.power(2, - (xlist) / tau)) * (1 +0 * np.cos(2 * np.pi * xlist/T + phi) )
        if type(extra_text) == bool: extra_text =""
        binned_rate_xy([chan, chan + 1], pk, xlist*y2s + time[0], 100 * (y-rate_exp)/rate_exp, 100 * yerr/rate_exp, rate_exp, res_list, fitted_pars = 4, sub = True, save = True,
            savename = "/fits/%i/sub_pk%i_t12_%i_a%.4f_phi%.2f_fix_%.4f_n_%i"%(CH, PK, TAU, AMP, PHI, FIX_A, FIT_NDIM) + datetime.datetime.today().isoformat(), 
            doplot = False, incr = 100, jenkins = True, bgrateplot = False, exp_text = True, lit_val = False, legend = True, coloritem = False, extra_text = extra_text + " Errors too small")
    
    # Return the fit results and the likelihood
    theta = [m.values[label] for label in ['amp','t12', 'phi', 'a'][0:ndim]]
    lnL   = lnlike(theta, x0, x1, y, yerr, tau12, dtau12, T, fixed_a)
    plt.close('all')
    return res_list, lnL





print("make_and_fit::\tLoaded data")

def make_fake_data(channel, pk, tau, amp, phi, dates = [time[0], time[-1]], T = 1, num_fake = 1):
    '''Generate a fake (random) dataset with a specified half-life, modulation- amplitude, phase and period
    reruns num_fake datasets with fake data'''
    print("\tmake_fake_data::\tdataset for source %s"%(sourceName[channel]))

    # Add the data from the channels
    chanlist = [channel, channel + 1]
    cuts = [cutter(ch_0 = ch, pk_0 = pk, t_start = dates[0], t_end = dates[1], T_high = THIGH, T_low = TLOW) for ch in chanlist]
    
    # The time-bin width
    delta_t = np.mean((time_end[cuts[0]] - time_start[cuts[0]]))
    # assert 1800 < delta_t < 3650, "The analyzer should produce a rate every hour"

    # The time-bin edges (in years)
    x0list, x1list = s2y * time_start[cuts[0]], s2y * time_end[cuts[0]]
    # And centers
    xlist = ((x0list + x1list) / 2 - x0list[0])
    # Remove t_initial from the data.
    x0list, x1list = x0list - x0list[0], x1list - x0list[0]
    # The rates and uncertainties
    ratelist = np.sum(rate[cuts[i]] for i in range(len(chanlist)))
    dratelist= np.sqrt(np.sum(drate[cuts[i]] ** 2 for i in range(len(chanlist))))

    # Temperature correction
    if TCOR: ratelist, dratelist = T_correct(ratelist, dratelist, temp[cuts[0]], 30, chanlist, pk)

    # Fit the overall amplitude of the data
    fit = use_minuit(x0list, x1list, ratelist, dratelist, tau, dhalflife[channel], channel, pk, ndim = 1, save = False)
    fit_A = fit[0][0][0]
    
    # Header of the MC/fake data
    fake_data_list = [['''this is a fake data list for channel %i, pk %i, tau = %g, mod_amp = %g, mod_phi = %g, mod_T = %g,
        between times'''%(channel, pk, tau, amp, phi, T), datetime.datetime.fromtimestamp(time[0])], ["t_start", "t_end", "rate", "drate"]]
    print("\t\tfitted A is %g Hz, delta T is %g" %(fit_A, delta_t))

    # The model
    rate_exp = (fit_A * np.power(2, - (xlist) / tau)) * (1 + amp * np.cos(2 * np.pi * xlist/T + phi) )
    
    # Make fake datasets. For every bin make a rate from a Poisson distribution, with an error given by the sqrt(N) statistics
    for i in range(num_fake):
        print("\tmake_fake_data::\tdataset\t%i"%i)
        rate_fake  = np.random.poisson(delta_t * rate_exp)
        drate_stat = np.sqrt(rate_fake)
        rate_fake  = rate_fake / delta_t
        drate_stat = drate_stat / delta_t
        fake_data_list.append([x0list, x1list, rate_fake, drate_stat])
    return fake_data_list

def T_correct(R, dR, T, Tset, ch, pk = -1):
    '''For a specified channel returns the parameters of a linear 
    correction to the data based on earlier experiments.'''
    dT = 0.5
    
    if type(ch) == list and len(ch) == 2:
        if 2 in ch and pk == 0:
            a  = -0.0158593141593
            da = 0.0068639952349
        elif 2 in ch and pk == 1:
            a  = -0.00469920847945
            da = 0.00306539446059
        elif 2 in ch and pk == 2:
            a  = -0.0030329862872
            da = 0.00227559699723
        elif 4 in ch and pk == 0:
            a  = -0.0123813164354
            da = 0.00292304568909
        elif 4 in ch and pk == 1:
            a  = -0.00894765762722
            da = 0.00270212060575
        elif 4 in ch and pk == 2:
            a  = -0.00194859226013
            da = 0.000865550759004
        elif 6 in ch and pk == 0:
            a  = -0.0064197104552
            da = 0.00455859751298
        else: 
            assert False, "channel not correct"

    R_correct  = R + a * (Tset - T)
    dR_correct = np.sqrt(dR ** 2 + (da * (Tset - T)) ** 2 + (dT * a) ** 2)
    
    return R_correct, dR_correct

def fit_real_data(ch, pk, dates = [time[0], time[-1]]):
    '''Takes the data list where each of the entries of that list contains four numpy arrays with
    t_start, t_end, the rate and drate. Returns the fitted value'''
    print("\tfit_real_data::\tStart")

    # Add the data from the channels
    chanlist = [ch, ch + 1]
    cuts = [cutter(ch_0 = chan, pk_0 = pk, t_start = dates[0], t_end = dates[1]) for chan in chanlist]
    

    # The time-bin edges (in years)
    x0list, x1list = s2y * time_start[cuts[0]], s2y * time_end[cuts[0]]
    # Remove t_initial from the data.
    x0list, x1list = x0list - x0list[0], x1list - x0list[0]
    # And centers
    xlist = (x0list + x1list) / 2
    # The rate and uncertainty
    ratelist = np.sum(rate[cuts[i]] for i in range(len(chanlist)))
    dratelist= np.sqrt(np.sum(drate[cuts[i]] ** 2 for i in range(len(chanlist))))
        
    # Temperature correction 
    if TCOR: ratelist, dratelist = T_correct(ratelist, dratelist, temp[cuts[0]], 30, chanlist, pk)
    
    result_list = []
    # The results from H1
    result = use_minuit(x0list, x1list, ratelist, dratelist, halflife[ch], dhalflife[ch], ch, pk, 
        ndim = FIT_NDIM, fixed_a= FIX_A, extra_text ="REAL DATA", save = True)

    LLH1 = result[1] 

    # The results from H0
    LLH0 = use_minuit(x0list, x1list, ratelist, dratelist, halflife[ch], dhalflife[ch], ch, pk, 
        ndim = 2, extra_text ="REAL DATA", save = True)[1]
    
    print("\t\tLLR result REAL DATA LL-H1: %f\tLL-H0: %f\tLLR %f"%(LLH1, LLH0, (LLH1-LLH0)))
    result_list.append(["DATA LLR ", LLH1-LLH0])

    print("\tmake_and_fit::\tfit_real_data\tDONE")
    return result_list

def fit_fake_data(ch, pk, tau, amp, phi, dates = [time[0], time[-1]], T = 1, num_fake = 1):
    '''Takes a fake data list where each of the entries of that list contains four numpy arrays with
    t_start, t_end, the rate and drate. Returns the fitted value'''
    print("\tfit_fake_data::\tStart")

    # Load fake datasets
    fake_data_list = make_fake_data(ch, pk, tau, amp, phi, dates = dates, T = T, num_fake = num_fake)
    

    # Add channels to get the rate
    chanlist = [ch, ch + 1]
    cuts = [cutter(ch_0 = chan, pk_0 = pk, t_start = dates[0], t_end = dates[1]) for chan in chanlist]
    ratelist = np.sum(rate[cuts[i]] for i in range(len(chanlist)))

    # Keep track of the results
    result_list = [[ch, pk, ratelist[0], tau, phi, amp]]
    
    names = {'A':ratelist[0], 'tau':tau, 'a':amp, 'phi':phi}
    for i in range(2, len(fake_data_list)):
        print("\tfit_fake_data::\tfit dataset\t%i"%(i-2))
        # The results from H0
        result = use_minuit(fake_data_list[i][0], fake_data_list[i][1], fake_data_list[i][2], 
            fake_data_list[i][3], halflife[ch], dhalflife[ch], ch, pk, ndim = 2, save = False)
        result, LLH0 = result

        # Print the result
        for j, name in enumerate(['A', 'tau']):
            print("\t\t%s \tlit=\t%g\t(%.3g\t+/-%.3g\t%.3g)"%(name, names[name], result[j][0], -result[j][1], + result[j][2]))
        
        # The results from H1 (and save some the first fit for a quality check)
        result = use_minuit(fake_data_list[i][0], fake_data_list[i][1], fake_data_list[i][2], 
            fake_data_list[i][3], halflife[ch], dhalflife[ch],
            ch, pk, ndim = FIT_NDIM, fixed_a= FIX_A, save = True if i == 2 else False)
        result, LLH1 = result
        
        # Add to what is being saved
        result_list.append(["\tLLR result ", LLH1-LLH0])

        # Print fit result
        for j, name in enumerate(['A', 'tau', 'phi', 'a']):
            print("\t\t%s \tlit=\t%g\t(%.3g\t+/-%.3g\t%.3g)"%(name, names[name], result[j][0], -result[j][1], + result[j][2]))
        print("\t\tLLR result fake LL-H1: %f\tLL-H0: %f\tLLR %f\n"%(LLH1, LLH0, (LLH1-LLH0)))
        
        # Save the fit result
        result_list.append(result)
    print("\tfit_fake_data::\tDone")    
    return result_list

def save_fitted_data():
    print("make_and_fit::\tsave_fitted_data\tStart")  
    fake_data_sets = []

    for number_times in range(0, NUM, max_fakes):
        do_x_times = max_fakes if number_times + max_fakes <= NUM else NUM%max_fakes
        print("make_and_fit::\tProcessing\t%i-%i (of %i total)"%(number_times, number_times + do_x_times, NUM))

        # Make datasets and save the results in fake_data_sets
        fake_data_sets = fit_fake_data(CH, PK, TAU, AMP, PHI, num_fake = do_x_times)
        # Also fit the real data
        fake_data_sets.append(fit_real_data(CH, PK))

        # Save the results to a text file that can be read bt logL.py
        print("make_and_fit::\tsave_fitted_data\tmade fits")
        textname = out_dir + 'fit_res_ch_%i_pk_%i_tau_%.1f_phi_%.1f_a_%.7f_fix_%.5f_ndim_%i.txt'%(CH, PK, TAU, PHI, AMP, float(FIX_A), FIT_NDIM)
        print("make_and_fit::\tsave_fitted_data\twrite data")
        file = open(textname, 'a')
        file.write('fitparameters:,')
        for data in fake_data_sets:
            for item  in data:
                if not type(item) == list: file.write(str(item) + ",")
                else: 
                    for subitem in item:
                        file.write(str(subitem) + ",")
            file.write("\n")
        file.close()
    print("make_and_fit::\tsave_fitted_data\tData written")

print("make_and_fit::\tRun code")
save_fitted_data()
print("make_and_fit::\tDONE, bye bye")
