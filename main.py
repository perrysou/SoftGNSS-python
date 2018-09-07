import initialize

# ./init.m

# --------------------------------------------------------------------------
#                           SoftGNSS v3.0
# 
# Copyright (C) Darius Plausinaitis and Dennis M. Akos
# Written by Darius Plausinaitis and Dennis M. Akos
# --------------------------------------------------------------------------
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
# USA.
# --------------------------------------------------------------------------

# Script initializes settings and environment of the software receiver.
# Then the processing is started.

# --------------------------------------------------------------------------


# Clean up the environment first =========================================
# clear
# close_('all')
# clc
# format('compact')
# format('long','g')
# --- Include folders with functions ---------------------------------------
# addpath('include')
# addpath('geoFunctions')
# Print startup ==========================================================
print '\n', \
    'Welcome to:  softGNSS\n\n', \
    'An open source GNSS SDR software project initiated by:\n\n', \
    '              Danish GPS Center/Aalborg University\n\n', \
    'The code was improved by GNSS Laboratory/University of Colorado.\n\n', \
    'The software receiver softGNSS comes with ABSOLUTELY NO WARRANTY;\n', \
    'for details please read license details in the file license.txt. This\n', \
    'is free software, and  you  are  welcome  to  redistribute  it under\n', \
    'the terms described in the license.\n\n', \
    '                   -------------------------------\n\n'
# Initialize settings class=========================================
settings = initialize.Settings()

# Generate plot of raw data and ask if ready to start processing =========
try:
    print 'Probing data "%s"...' % settings.fileName
    settings.probeData()
    settings.probeData('/Users/yangsu/Downloads/GNSS_signal_records/GPS_and_GIOVE_A-NN-fs16_3676-if4_1304.bin')
finally:
    pass

print '  Raw IF data plotted '
print '  (run setSettings or change settings in "initialize.py" to reconfigure)'
print ' '
gnssStart = True
# gnssStart = int(raw_input('Enter "1" to initiate GNSS processing or "0" to exit : ').strip())

if gnssStart:
    print ' '
    settings.postProcessing()
