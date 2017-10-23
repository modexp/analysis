#!/usr/bin/env python
#
# Adam Brown, December 2016
# abrown@physik.uzh.ch

print("Welcome to the slow control plot updater")

import os
import sys
from datetime import datetime, timedelta
import time
import json
import configparser

import matplotlib
matplotlib.use('Agg') # Don't try to use X forwarding for plots
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import root_pandas

# How far back to plot in hours
config = configparser.ConfigParser()
config.read('config.ini')
plotrange1 = config.getint('SlowControl', 'plotrange1')
plotrange2 = config.getint('SlowControl', 'plotrange2')

# These to show whether data is good or bad
def getconfigrange(section, value):
    return eval(config[section][value])

pressrange = getconfigrange('AlarmRanges', 'pressure')
temprange = getconfigrange('AlarmRanges', 'temperature')
humidrange = getconfigrange('AlarmRanges', 'humidity')
magrange = getconfigrange('AlarmRanges', 'magfield')
radonrange = getconfigrange('AlarmRanges', 'radon')
datadelayrange = getconfigrange('AlarmRanges', 'datadelay')
maxUpdateDelay = getconfigrange('AlarmRanges', 'maxupdatedelay')

xFormatter  = matplotlib.dates.DateFormatter('%d-%m\n%H:%M')
yFormatter = matplotlib.ticker.ScalarFormatter(useOffset=False)

plotDirectory = config['SlowControl']['plotoutdir']
warningFilename = plotDirectory + '/warning.json'
processedDataDir = os.environ['MODEXP_PROCESSED_DATA_DIR']
rawDataDir = os.environ['MODEXP_RAW_DATA_DIR']
anaDataDir = os.environ['MODEXP_ANALYSIS_DATA_DIR']
run_dir = os.environ['MODEXP_TEMP_SCRIPTS_DIR']

timeRange1 = (datetime.utcnow() - timedelta(hours = plotrange1), datetime.utcnow())
timeRange2 = (datetime.utcnow() - timedelta(hours = plotrange2), datetime.utcnow())

# Tests if a timestamp is with the useful range
def isRecent(timestamp, timeLimit):
    diff = datetime.now() - datetime.fromtimestamp(timestamp)
    return diff < timedelta(hours = timeLimit)

# See if directory exists and make it if not
def ensureDir(d):
    if not os.path.exists(d):
        os.makedirs(d)

# Find recent files with particular extension
def getRecentFiles(parentDir, fileExtension, timeLimit):
    retFileList = []
    for (dir, _, fileList) in os.walk(parentDir):
        for file in fileList:
            _, extension = os.path.splitext(file)
            fullPath = dir + '/' + file
            if extension == '.' + fileExtension and isRecent(os.path.getmtime(fullPath), timeLimit):
                retFileList.append(fullPath)
    return retFileList

# Get most recent file in a directory tree
def getMostRecentFile(parentDir, fileExtension):
    retFile = None
    retFileTime = None
    for (dir, _, fileList) in os.walk(parentDir):
        for file in fileList:
            _, extension = os.path.splitext(file)
            fullPath = dir + '/' + file
            if extension == '.' + fileExtension and (retFile == None or os.path.getmtime(fullPath) > retFileTime):
                retFile = fullPath
                retFileTime = os.path.getmtime(fullPath)
    return retFile

# Find all recent slow data (within longest plot range)
slowDataFiles = getRecentFiles(processedDataDir, 'sroot', plotrange2)

# If no slow data files found exit
if len(slowDataFiles) == 0:
    print("No slow data found, nothing to do")
    sys.exit(0)

# Load slow data into pandas dataframe
slowdata = pd.concat([root_pandas.read_root(file) for file in slowDataFiles])
slowdata['unix_time'] = slowdata['stime'] - 2208988800 # Convert from Mac to UNIX time
slowdata['pandas_time'] = pd.to_datetime(slowdata['unix_time'], unit = 's')

# Trim slow data to only 'interesting' times (within one plot range)
sd_filter = (slowdata['pandas_time'] > timeRange1[0]) | (slowdata['pandas_time'] > timeRange2[0])
slowdata = slowdata[sd_filter]

# What to plot and what to call it
dataIndices = ['temp', 'pres', 'humid', 'btot', 'radon']
axisLabels = ["Temperature [Deg]", "Pressure [Pa]", "Humidity [%]", "Magnetic field [gauss?]", "Radon"]
plotTitles = ["Temperature", "Pressure", "Humidity", "Magnetic field", "Radon"]
fileNameStems = ['temperature', 'pressure', 'humidity', 'magfield', 'radon']

# Loop through variables and plot
for (index, yLabel, title, filename) in zip(dataIndices, axisLabels, plotTitles, fileNameStems):
    print("Plotting %s..." % title)

    fig = plt.figure(figsize=(10,3))
    plt.plot_date(slowdata['pandas_time'], slowdata[index], 'b.')
    
    # Format axes
    ax = plt.gca()
    ax.xaxis.set_major_formatter(xFormatter)
    ax.yaxis.set_major_formatter(yFormatter)

    # Labels
    plt.ylabel(yLabel)
    plt.xlabel("Time")
    
    # For each time range save the plot with correct x-axis range
    plt.title(title + " " + str(plotrange1) + " hours")
    plt.xlim(timeRange1)
    plt.tight_layout()
    fig.savefig(plotDirectory + '/' + filename + '1.png')
    
    plt.title(title + " " + str(plotrange2) + " hours")
    plt.xlim(timeRange2)
    fig.savefig(plotDirectory + '/' + filename + '2.png')


# Plot HV seperately so they can all be together on one graph
print("Plotting high voltages")
fig = plt.figure(figsize=(10, 4))
for i in range(8):
    dataIndex = 'hv%d' % i
    plt.plot_date(slowdata['pandas_time'], slowdata[dataIndex], '.', label=("HV %d" % i))

# Format axes
ax = plt.gca()
ax.xaxis.set_major_formatter(xFormatter)
ax.yaxis.set_major_formatter(yFormatter)

# Labels
plt.ylabel("Voltage [V]")
plt.xlabel("Time")

# For each time range save the plot with correct x-axis range
filename = 'highvoltage'

plt.xlim(timeRange1)

# Resize axes and put legend in the space created
plt.tight_layout()
box = ax.get_position()
ax.set_position([box.x0, box.y0,
                 box.width, box.height * 0.75])
plt.legend(loc='upper center', ncol=4, borderaxespad=0.0,
           bbox_to_anchor=(0.0, 1.0, 1.0, 0.333), mode='expand')

plt.xlim(timeRange1)
fig.savefig(plotDirectory + '/' + filename + '1.png')

plt.xlim(timeRange2)
fig.savefig(plotDirectory + '/' + filename + '2.png')


### Make status table in html
print("Making status table")

# Helper functions
def goodMsg(status):
  if status:
    return "Good"
  else:
    return "Bad"

def goodColor(status):
  if status:
    return "#c1ffd5"
  else:
    return "#ffc1c8"

# Store warnings while we create the table to save later to json
warnlist = []

def tableRow(name, value, allowedrange, precision):
  status = (value > allowedrange[0]) and (value < allowedrange[1])
  if not status:
      warnlist.append({'category': 'warning',
                       'type': 'range',
                       'variable': name,
                       'value': value,
                       'allowed_range': allowedrange})
  return "<tr bgcolor=%s><td>%s</td><td>%.*f</td><td>%.*f to %.*f</td><td>%s</td></tr>" \
                 % (goodColor(status), name, precision, value, precision, allowedrange[0], precision, allowedrange[1], goodMsg(status))

# Get the latest row
numRows = slowdata.shape[0]
last_row = slowdata.iloc[0]
last_time = slowdata.max()['stime']
last_row = slowdata.loc[slowdata['stime'] == last_time].iloc[0]
lastupdatedelay = (datetime.utcnow() - last_row['pandas_time']) / np.timedelta64(1, 's')

# Save the latest values to a file
latestfile = open(plotDirectory+'/datatable.php', 'w')
latestfile.write(tableRow("Time since last slow data [s]", lastupdatedelay, datadelayrange, 0))
latestfile.write(tableRow("Pressure [Pa]", last_row['pres'], pressrange, 0))
latestfile.write(tableRow("Temperature [deg C]", last_row['temp'], temprange, 1))
latestfile.write(tableRow("Magnetic field [???]", last_row['btot'], magrange, 1))
latestfile.write(tableRow("Humidity [%]", last_row['humid'], humidrange, 1))

# Use PHP to see whether this script itself ran recently enough
latestfile.write("<?php if(time() < %f): ?>" % (time.time() + maxUpdateDelay * 60))
latestfile.write("<tr bgcolor=%s><td>Last slow control update</td><td>%s</td><td>Within %.0f minutes</td><td>%s</td></tr>"
                 % (goodColor(True), datetime.utcnow().strftime("%d.%m.%y %H:%M"), maxUpdateDelay, goodMsg(True)))
latestfile.write("<?php else: ?>")
latestfile.write("<tr bgcolor=%s><td>Last slow control update</td><td>%s</td><td>Within %.0f minutes</td><td>%s</td></tr>"
                 % (goodColor(False), datetime.utcnow().strftime("%d.%m.%y %H:%M"), maxUpdateDelay, goodMsg(False)))
latestfile.write("<?php endif ?>")
latestfile.close()

# Make new table with most recently modified data files
datastatfile = open(plotDirectory + '/datastat.php', 'w')
recent_bin = getMostRecentFile(rawDataDir, 'bin')
recent_slo = getMostRecentFile(rawDataDir, 'slo')
recent_ana = getMostRecentFile(anaDataDir, 'root')
datastatfile.write("<p>Most recent raw data file: %s<br />" % recent_bin)
datastatfile.write("Most recent slow data file: %s<br />" % recent_slo)
datastatfile.write("Most recent analysed data: %s</p>" % recent_ana) 
datastatfile.close()

### Save any warnings to the json
with open(warningFilename, 'w') as warnfile:
    warnfile.write(json.dumps(warnlist))

### Update spectra for most recent analysed data
print("Updating spectra")
for channel in range(0, 8):
    rootcommand = 'monitor.C("%s", "spectrum", %d, true, true)' % (recent_ana, channel)
    commandstring = 'cd $MODEXP_ANALYSIS_DIR/monitor; '
    commandstring += 'root -q -b \'%s\'' % rootcommand
    os.system(commandstring)

cmd = 'cp $MODEXP_ANALYSIS_DIR/monitor/plots/spectrum_channel*_log.png ' + plotDirectory
os.system(cmd)

print("Plots updated. Done.")
