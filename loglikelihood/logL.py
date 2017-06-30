#!/project/datagrid/anaconda/envs/joran/bin/python
# -*- coding: utf-8 -*-

############################################
# written by
# Joran Angevaare 4-5-2017
# jorang@xs4all.nl
############################################

############################################
# Interprets the results produced by make_and_fit.py for the log-likelihood hypthesis test
# Expects that there is already processed data on the out_dir that has been written there by the 
# make_and_fit.py
############################################


############################################
# Run options
# My username on the Nikhef cluster 
username = 'jorana'
institute='nikhef'

############################################


import time as pythontime         # to wait time
import numpy as np                  # Fast
import ROOT                         # To open root files
import matplotlib.cm as cm          # colorbars 
import matplotlib                   # Plotting

import datetime                     # For dates on axis 
from os import listdir              # To see what is in a folder
import os
from matplotlib.ticker import MaxNLocator
matplotlib.use('Agg')
import matplotlib.pyplot as plt     # Plotting
matplotlib.rc('font', size = 24)                   # Use big fonts...


# Directories
data_dir = "/data/modulation/Raw_Data/combined/"
anaf_dir = "/dcache/xenon/jorana/Modulation/processed_use_thesis/analysis_used/"
cali_dir = "/dcache/xenon/jorana/Modulation/processed_use_thesis/calibration_used/"
stbc_dir = "/data/xenon/mod_admin/Modulation/stoomboot/"
logL_dir = "/data/xenon/mod_admin/Modulation/loglikelihood/"
out_dir  = "/data/xenon/mod_admin/data/logL/" 

# Conversions
y2s     = 365.25 * 24 * 3600 # years to seconds
s2y     = 1. / y2s         # seconds to years
toffset = 70 * y2s       # the timestamp from labview is 70 years off

# Source energies sorted by channel, peak
sourceE = [[1460], # ch0 BG
           [1460], # ch1 BG 
           [511,   1157.020, 511 + 1157.020],  # ch2 Ti
           [511,   1157.020, 511 + 1157.020],  # ch3 Ti
           [1173.2,1332.5,   1173.2 + 1332.5], # ch4 Co          
           [1173.2,1332.5,   1173.2 + 1332.5], # ch5 Co
           [661.7], # ch6 Cs
           [661.7]] # ch7 Cs

# A list of the halflife of the source format halflife[channel] = tau1/2
halflife = [1e9,    1e9,    59.1,   59.1,   5.2711, 5.2711,     30.08,  30.08]
dhalflife= [1e9,    1e9,    0.3,    0.3,    0.0004, 0.0004,     0.09,   0.09]

# Names of the sources, the index of the list corresponds to the channel
sourceName = ['Background', 'Background',   'Ti-44',    'Ti-44',    'Co-60', 'Co-60', 'Cs-137', 'Cs-137']
sourceLatex =['Background', 'Background', '$^{44}$Ti', '$^{44}$Ti', '$^{60}$Co', '$^{60}$Co', '$^{137}$Cs', '$^{137}$Cs']


col_lis = ['green', 'red', 'blue', 'cyan', 'brown', 'pink', 'magenta','black']         


############################################
def gaussian(x, mu, sigma):
    return np.exp((-(x-mu)**2)/(2 * sigma ** 2)) / (np.sqrt(2*np.pi) * sigma)

def num(s):
    '''Converts string to int/float'''
    if s == '\n': return None
    try: return int(s)
    except ValueError:
        try: return float(s)
        except ValueError:
            print("logL.py::\tnum::\tcould not convert '%s' to float"%s)
            return None

def fit_hist(entry, nbins, xlabel):
    ''''Makes a very simple histogram of the entries with nbins and a label'''
    bins = np.linspace(entry.min(), entry.max(), nbins)
    hist, bin_edges = np.histogram(entry, bins=nbins)
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:]) 
    # Plot the histogram as points with an errorbar 
    dx = bin_edges[1]-bin_edges[0]
    plt.errorbar(bin_centers, hist, yerr = np.sqrt(hist), xerr = [dx for every in bin_centers],
        linestyle = 'None', marker = 'o', color = 'r', 
        label = 'histogram $\mu\pm\sigma =%.2g\pm%.2g$'%(np.mean(entry), np.std(entry)))
    plt.ylabel('Counts')
    plt.xlabel(xlabel)
    lgnd2 = plt.legend(bbox_to_anchor=(1,1), loc=2, borderaxespad=0., markerscale = 1, numpoints = 1)


def check_qstat():
    '''Checks how many files are running on the cluster for user: username'''
    # make a temporary file where the output of the linux command is written.
    TMPFILE = out_dir + '/qlist.tmp' 
    cmd = 'qstat -u '+ username +' | grep mc_ | grep -c R > ' + TMPFILE
    os.system(cmd)
    f = open(TMPFILE, 'r')
    n = 0
    for line in f:
        n = int(line)
    cmd = 'rm ' + TMPFILE
    os.system(cmd)
    return n


def cor_phi(phi,sign):
    '''The fit results from make_and_fit.py are degenerate in phi (the modulation phase). This corrects for that.'''
    returnlis = []
    for i in range(len(phi)):
        retphi = phi[i]
        if sign[i] == -1: retphi += + np.pi
        while retphi> 2*np.pi: retphi -= 2 * np.pi
        while retphi< 0: retphi += 2 * np.pi
        returnlis.append(retphi)
    return np.array(returnlis)

def show_fitted_data(data):
    '''This opens the data that contains the headers and fit results and makes a histogram of the 
    values that are fitted an compares that to the value that was added to the Monte Carlo'''

    # fake_data_sets = data
    # the first line in data contains the values that went into the MC
    ch, pk, A_real, tau_real, phi_real, a_real = data[0]       
    # everything else are the fit results
    toplot = data[1:len(data)]

    ch, pk = int(ch), int(pk)
    
    # Sometimes the datafile is corrupt and the structure is incorrect. If that is the case we do 
    # not proceed if it works it gives us a list of the sign of the fitted modulation amplitude.
    try: signs = np.array([toplot[i][3][0]/abs(toplot[i][3][0]) for i in range(len(toplot))])
    except IndexError:
        print("logL.py::\tshow_fitted_data::\tIndexError for, num %i. The line reads:"%(i))
        print(toplot[i])
        return None

    # In the data there are 
        # A:   overall amplitude
        # tau: the half-life
        # phi: the modulation phase
        # a:   the modulation amplitude
    # we loop over these to make a plot of the fit results
    for j, name in enumerate(['A', 'tau', 'phi', 'a']):

        ys = np.linspace(0,len(toplot), len(toplot))   
        try: xs = np.array([toplot[i][j][0] for i in range(len(toplot))])
        except IndexError:
            print("logL.py::\tshow_fitted_data::\tIndexError for, num %i. The line reads:"%(i))
            print(toplot[i])
            return None
        # The fitted values of phi are degenerate, correct for that degeneracy
        if name == 'phi': xs = cor_phi(xs, signs)            
        # The negative values of the modulation amplitude are meaningless, plot the absolute value    
        if name == 'a': xs = abs(xs)
        
        plt.xlim(min(xs), max(xs))
        xerr = np.transpose([toplot[i][j][1:3] for i in range(len(toplot))])

        plt.title("%s, peak %f"%(sourceLatex[ch], sourceE[ch][pk]))
        
        # make an arbitrary amount of bins
        nbins = 5 * int(np.sqrt(len(xs)))
        while nbins > len(xs): nbins = int(np.sqrt(nbins))
        
        # the plotting happens in the fit_hist
        fit_hist(xs, nbins / 2, name)            

        # Add a mean and the assumed value of the MC
        plt.axvline(data[0][j + 2], label = "real value of %s"%name, c = 'black', linewidth = 4)
        plt.axvline(np.average(xs), label = "average of %s"%name,    c = 'blue',  linewidth = 4)

        lgnd1 = plt.legend(bbox_to_anchor=(1,1), loc=2, borderaxespad=0., markerscale = 1, numpoints = 1)
        plt.savefig(out_dir + '/figures/%s/fit_res_ch_%i_pk_%i_tau_%.3f_phi_%.1f_a_%.7f_%s.png'%(name, ch, pk, tau_real, phi_real, a_real, file_name), dpi = 300, bbox_extra_artists=(lgnd1,), bbox_inches='tight')
        plt.clf()

def interpretLLR(LLR_results, LLR_data):
    '''Determines the p-value of the H0 hypothesis and plots the LLR (log-likelihood-ratio) of a given simulation'''
    LLR_results = np.array([num(LLR) for LLR in LLR_results])
    LLR_data    = num(LLR_data)
    
    nbins = 250
    values, bin_edges, some = plt.hist(LLR_results, bins = nbins)
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])

    # The p-value of the H0-hypothesis (no modulation) is the probability of the simulated data having a LLR that is
    # greater or equal to the LLR of the real data
    tot_num = float(len(LLR_results))
    big_num = float(len(LLR_results[LLR_results>LLR_data])) # the number of simulations that have a smaller LLR than the data
    pvalue  = big_num / tot_num
    dpvalue = np.sqrt(pvalue/tot_num + pvalue * pvalue / tot_num)

    plt.ylabel("Counts/bin")
    plt.xlabel("Log-likelihood ratio")
    plt.axvline(LLR_data, lw = 5, color = 'black')
    
    plt.title("H_0\nP-value = %g +/- %g"%(pvalue, dpvalue)) 
    plt.yscale("log")
    plt.savefig(out_dir + '/figures/LLR/LLR_res_%s.png'%(file_name), dpi = 300, bbox_inches='tight')
    plt.clf()

    print("logL.py::\tinterpretLLR::\tmean LLR \t%.2f" %np.mean(LLR_results) + "\tLLR data \t%.2f\tp-value\t%.2g"%(LLR_data, pvalue))
    return pvalue, dpvalue

def plot_pval(a, p, dp, name):
    '''Save the found p-values in a small text file and make a plot of the p-value as function of the modulation amplitude'''
    textname = out_dir + name +".txt"
    file = open(textname, 'w')
    file.write(name + " a =\n")
    # Write to this file the the modulation amplitude (a), the p-value (p) and the uncertainty on p (dp)
    for i in a: file.write("%.8f"%i + ", ")
    file.write("\np =\n")
    for i in p: file.write("%.8f"%i + ", ")
    file.write("\ndp =\n")
    for i in dp: file.write("%.8f"%i + ", ")
    file.close()

    # Plotting details
    plt.title("p-value of H_0")
    plt.errorbar(a,p, yerr= dp, linestyle = "None", marker = "o")
    plt.ylabel("p-value")
    plt.xlabel("Simulated value of a")
    plt.savefig(out_dir + '/figures/%s.png'%(name), dpi = 300, bbox_inches='tight')
    plt.clf()



############################################
# Start of the main program
############################################

# Wait until all jobs are finished running 
if institute == 'nikhef':
    for i in range(10):
        if check_qstat() == 0: continue
        else: print("logL.py::\twaited \t%i hours"%i); pythontime.sleep(60*60)

print('logL.py::\tstart')

############################################
# Check that all folders exist
figfolder = out_dir+"/figures/"
if not os.path.exists(figfolder): 
    print('logL.py::\tmake folder')    
    cmd = "mkdir " + figfolder
    print(cmd)
    os.system(cmd)

for string in ['A', 'tau', 'phi', 'a', 'LLR']:
    if not os.path.exists(figfolder+string+"/"): 
        print('logL.py::\tmake folder')    
        cmd = "mkdir " + figfolder +string+"/"
        print(cmd)
        os.system(cmd)
############################################


############################################
# Open the files in out_dir that have been produced by make_and_fit.py
fit_files = sorted(listdir(out_dir))
# Keep a dictionary of which p-values belong to which value of the simulated modulation amplitude
Pdict = {}
for file_name in fit_files:
    
    LLRdirct = {}
    if 'fit_res_ch_' in file_name and '.txt' in file_name:
        print("\nlogL.py::\tstart with file " + file_name)
        file = open(out_dir + file_name, 'r')
        bool_header = False
        data        = []
        LLR_results = []

        # Loop over the lines in the file, these van contain either the input parameters of the MC
        # the result of the LLR for the simulated data or the LLR for the real data or it are the 
        # fit results of the simulated data
        for line in file:
            line_list = line.split(',')
            if 'fitparameters:' in line: 
                if not bool_header: data.append([num(line_list[i]) for i in range(1, len(line_list)-1)])
                bool_header = True
                plist = [num(line_list[1]),num(line_list[2]), num(line_list[4])]

                if not None in plist: 
                    thisP = "Ps_ch_%i_%i_t12_%.f"%(plist[0], plist[1], plist[2])
                    this_a = num(line_list[6])
                else: print("logL.py::\tSkipped this line. There is a 'None' value in the LLR file")
                
            elif 'LLR result' in line: LLR_results.append(line_list[1])
            elif "DATA LLR"   in line: LLR_data = ["DATA LLR", line_list[1]]
            else:  
                # it are the fit results of the simulated data
                this_list = []
                for i in range(0, len(line_list), 3):
                    if len(line_list[i:-1]) < 3: continue
                    toappend = [num(line_list[j]) for j in range(i, i + 3)]
                    if not None in toappend: this_list.append(toappend)
                data.append(this_list)
                    
        file.close()

        # Make a plot of the fit results
        show_fitted_data(data)

        # Get the p-value of the data
        pval, dpval = interpretLLR(LLR_results, LLR_data[1])
        try: 
            Pdict[thisP][0].append(this_a)
            Pdict[thisP][1].append(pval)
            Pdict[thisP][2].append(dpval)
        except KeyError:
            Pdict[thisP] = [[this_a], [pval], [dpval], thisP]

for pitem in Pdict:
    plot_pval(*Pdict[pitem])

print('logL.py::\tDone! Bye bye')






















# print("LogL::\t run code")
# show_fitted_data()

# print("LogL::\tDONE, bye bye")
