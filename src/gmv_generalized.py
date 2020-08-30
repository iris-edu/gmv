#!/usr/bin/env python

import sys
import os
import time
import collections
import numpy as np
from collections import OrderedDict
import datetime
import operator

import getopt

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import animation
from scipy.interpolate import griddata

from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics.base import kilometers2degrees
from obspy.clients.fdsn import Client
import json

import gmv_utils as utils

from PIL import Image
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

from mpl_toolkits.basemap import Basemap, maskoceans

import gmv_param as param

"""
    Name:
    gmv_generalized.py
    
    Description:
    
    This is a Python version of the  Ground Motion Visualization (http://ds.iris.edu/ds/products/gmv/) code in MATLAB:
    Credits
          - Chuck Ammon, Professor of Geosciences at Penn State’s original concept and visualizations.
          - Bob Woodward at IRIS – adapted the visualization code to MATLAB.
          - IRIS DMC Data products expanded and enhanced the MATLAB code.

    This script is capable of producing GMVs for different geographical regions, different sensor technologies and 
    combined or super GMVs (http://ds.iris.edu/ds/products/usarraygmv-super/).

    Script can be configured via its parameter file gmv_param.par or via the command line arguments. Currently
    parameters are optimized for use with  Lambert conformal map projection and seismic channels. Change in projection
    and/or technology may require parameter tuning.

    Copyright:
    Copyright (C) 2020  Product Team, IRIS Data Management Center

    This is a free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 3 of the
    License, or (at your option) any later version.

    This script is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License (GNU-LGPL) for more details.  The
    GNU-LGPL and further information can be found here:
    http://www.gnu.org/

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    History:
        2020-09-03 Manoch: V.2020.247 R1.1 Public release.
        2020-07-27 Manoch: V.2020.209 added chunking option for data request
        2020-07-08 Manoch: V.2020.191 added script version to production label
        2020-07-06 Manoch: V.2020.189 introduced North Polar region and STD water level 
        2020-05-30 Manoch: V.2020.152 R1 release
        2020-05-25 Manoch: V.2020.146, Added channel label instead of Primary, ...
        2020-05-11 Manoch: V.2020.132, Initial release.
"""
script_version = 'V.2020.247'
script = sys.argv[0]
script = os.path.basename(script)


def usage():
    """The usage message.
    """
    new_line = '\n'
    print(f'{new_line}{new_line}{script} ({script_version}):')
    print(f'{new_line}This is a Python version of the Ground Motion Visualization (GMV) code{new_line}'
          f'that creates visualizations of real data showing how seismic waves from earthquakes sweep across the '
          f'network of seismic stations for which data are openly available and are collected using the federated '
          f'data access.http://ds.iris.edu/ds/products/gmv/'
          f'Credits:{new_line}'
          f'\t- Chuck Ammon, Professor of Geosciences at Penn State’s original concept and visualizations{new_line}'
          f'\t- Bob Woodward at IRIS – adapted the visualization code to MATLAB{new_line}'
          f'\t- IRIS DMC Data products expanded and enhanced the MATLAB code.{new_line}{new_line}'
          f'This script is capable of producing GMVs for different geographical regions, sensors, technologies'
          f'{new_line}'
          f'and combined or super GMVs (http://ds.iris.edu/ds/products/usarraygmv-super/).{new_line}'
          f'Script can be configured via its parameter file (gmv_param.par) or command line arguments. '
          f'{new_line}'
          f'Currently parameters are optimized for use with the Lambert conformal map projection and seismic '
          f'channels.{new_line}'
          f'Change in projection and/or technology may require additional parameter tuning.{new_line}'
          f'{new_line}')
    print(f' -h  --help\t\tOutput this message.{new_line}',
          f'-v, --verbose\t\t[default: {param.verbose}] Turn on the verbose mode.{new_line}',
          f'-b, --band\t\t[default: {param.channel_band}] Channel bands to use for GMV production (separate by comma, '
          f'for example: '
          f'LH OR LH,BH OR  LH,BH,HH).{new_line}',
          f'-c, --comp\t\t[default: {param.request_comp}] Number of channel components (1 or 3).{new_line}',
          f'-d, --dur\t\t[default: {param.animation_duration}] Duration of the GMV animation in seconds.{new_line}',
          f'-D, --std\t\t[default: {param.std_max}] Maximum acceptable standard '
          f'deviation for trace vetting.{new_line}',
          f'-e, --eloc\t\t[*required] Single quoted event location as \'lat,lon\' (example:  \'10.779,-62.907\'). For '
          f'multiple events., '
          f'{new_line}\t\t\tseparate each set by a space (example: \'24.760,-109.890 25.200,-109.870\'. '
          f'Use of parentheses is optional: \'(24.760,-109.890) (25.200,-109.870)\'.{new_line}',
          f'-g, --gain\t\t[default: {param.gain}]Trace amplification to generate GMVs. For 3C GMVs, gain is only '
          f'applied to the Z-component.{new_line}',
          f'-G, --gc\t\t[default: {param.draw_great_circle}] raw the great circle path between the event location and '
          f'the reference station.{new_line}',
          f'-l, --tstep\t\t[default: {param.time_step}] Time step in seconds to use for sampling '
          f'traces and create the video. '
          f'Example: -l 2 samples traces every 2 seconds.{new_line}',
          f'-m, --emag\t\t[*required] Event magnitude. For multiple events, single quoted magnitude list and separate '
          f'them with  a space (example: -m \'7.3 6.2\'.{new_line}',
          f'-n, --net\t\t[default: {param.network.replace("*", "all")}] Network(s) to request data from. For '
          f'multiple networks, separate them with comma. '
          f'Use all to request from all networks (example: -n TA OR -n US,TA -n all).{new_line}',
          f'-N, --rnet\t\t[default: {param.proj_regions[param.region]["ref_sta"][0]}] Network of the reference station.'
          f'Example -N TA.{new_line}',
          f'-o, --output\t\tThe output video file name (by default the video will be MP4 format. '
          f'Example: -o Test_video will output Test_video.mp4.{new_line}',
          f'-p, --delay\t\t[default: {param.animation_delay}] Delay in seconds from the event origin time '
          f'to start the video. Example: -p 120.{new_line}',
          f'-P, --phase\t\t[default: {param.phase_spacing}] Seconds between Phases marked on the trace. '
          f'This is used to avoid overprinting phase labels. '
          f'Example -P 30.{new_line}',
          f'-q, --qscale\t\t[default: {param.proj_regions[param.region]["quiver_scale"]}] Quiver scale for 3C plots. '
          f'Number of data units per arrow length unit, e.g., m/s per '
          f'plot width; a smaller scale parameter makes the arrow longer.{new_line}',
          f'-r, --region\t\t[default: {param.region}] The region to create the GMV for. The selected region '
          f'must be a key of the proj_regions. The acceptable regions are: {param.proj_regions.keys()}.'
          f'{new_line}',
          f'-s, --sizem\t\t[default: {param.proj_regions[param.region]["marker_size"]}] Marker '
          f'(https://matplotlib.org/api/markers_api.html) size in points.{new_line}',
          f'-S, --rsta\t\t[default: {param.proj_regions[param.region]["ref_sta"][1]}] Reference station code '
          f'to plot.{new_line}',
          f'-t, --etime\t\t[*required] The event time as YYYY-MM-DDTHH:MM:SS.{new_line}',
          f'-T, --title\t\t[default: based on the event magnitude (-m) and location (-e)] '
          f'Single quoted GMV title.{new_line}',
          f'-z, --depth\t\t[*required] Event depth in km.{new_line}',
          f'{new_line}For additional configuration see the parameter file (gmv_param.py)')
    print("\n\nExample:\n"
          "\n\t1. Sample request with the least number of arguments:\n\t\tgmv_generalized.py -e 55.1046,-158.4725 "
          "-z 10.0 -m 7.8 -t 2020-07-22T06:12:42 -o GMV_Example_Default"
          "\n\t2. Sample complete request:\n\t\tgmv_generalized.py --band=LH,BH --comp=1 -n all -t "
          "2020-07-22T06:12:42 -T "
          "'July 22, 2020, Alaska Peninsula, M 7.8' -m 7.8 -z 10.0 -e 55.1046,-158.4725 -r ak -d 1200 "
          "-s 6.0 -p -180 -q 3.5 -g 3 -D 0.05 -G -o GMV_Example_Custom"
          )
    print('\n\n\n')


def get_map_params(llc_point, urc_point):
    """Compute center point, width and height for the base map based on the LL and UR coordinates
    """

    _mid_latitude = (llc_point[0] + urc_point[0]) / 2.0
    _mid_longitude = (llc_point[1] + urc_point[1]) / 2.0
    _width, _azim, _back_azim = gps2dist_azimuth(float(llc_point[0]), float(llc_point[1]),
                                                 float(llc_point[0]), float(urc_point[1]))
    _height, _azim, _back_azim = gps2dist_azimuth(float(llc_point[0]), float(llc_point[1]),
                                                  float(urc_point[0]), float(llc_point[1]))

    return (_mid_latitude, _mid_longitude), _width, _height


# Set up formatting for the movie files.
t_run = time.time()

# ----------------- PARAMETERS START -----------------

url_data_centers = param.url_data_centers

# These parameters are used for subplots layout. The user is expected to provide map width and hight. Using the below
# values, we can calculate the image size based on the map size.

metadata_dir = param.metadata_dir
metadata = list()
metadata_filter = False

video_dir = utils.mkdir(param.video_dir)
log_dir = utils.mkdir(param.log_dir)
metadata_dir = utils.mkdir(metadata_dir)

figure_width = param.figure_width
subplot_columns = param.subplot_columns
trace_height = param.trace_height
subplot_rows = param.subplot_rows

logo_dir = utils.mkdir(param.image_dir)
logo_file = os.path.join(logo_dir, param.logo_file)
add_log = os.path.isfile(logo_file)
logo_zoom = param.logo_zoom
logo_location = param.logo_location
frame_seconds = param.frame_seconds

proj_regions = param.proj_regions

output_type = param.output_type
output_factor_symbol = param.output_factor_symbol
output_units = param.output_units
channel_comp_index = param.channel_comp_index
channel_band_list = param.channel_band_list
channel_comp_list = param.channel_comp_list

production_label = param.production_label
production_label = f'{production_label} {script_version}'

no_filter_color = param.no_filter_color
no_filter_text = param.no_filter_text
horizontal_motion_text = param.horizontal_motion_text
horizontal_motion_color = param.horizontal_motion_color

request_bands = channel_band_list
request_comp = param.request_comp

chan_list = list()
for band in request_bands:
    for _comp in channel_comp_list[request_comp]:
        chan_list.append(f'{band}{_comp}')
request_chans = ','.join(chan_list)

request_net = param.request_net

gain = param.gain

subplot_map_rows = param.subplot_map_rows
subplot_trace_rows = param.subplot_trace_rows
subplot_rows = param.subplot_rows

technology = param.technology
output_type = param.output_type
draw_great_circle = param.draw_great_circle
great_circle_line_color = param.great_circle_line_color

zero_base = param.zero_base
scale_plot = param.scale_plot

ref_sta_marker_color = param.ref_sta_marker_color
ref_sta_marker_size = param.ref_sta_marker_size
ref_sta_marker_edge_color = param.ref_sta_marker_edge_color
ref_sta_marker_edge_width = param.ref_sta_marker_edge_width

region = param.region
marker_size = proj_regions[region]['marker_size']
quiver_scale = proj_regions[region]['quiver_scale']
quiver_scale_units = proj_regions[region]['quiver_scale_units']

figure_aspect_ratio = float(proj_regions[region]['size'][0] / proj_regions[region]['size'][1])

marker_event = param.marker_event
marker_track = param.marker_track
marker_size_event = param.marker_size_event
marker_size_track = param.marker_size_track

# Sign of the delay: -X start X seconds sooner; +x start X seconds later
animation_delay = param.animation_delay
animation_duration = param.animation_duration
phase_spacing = param.phase_spacing

chunk_count = param.chunk_count


contour_handle = list()
author = param.author
bit_rate = param.bit_rate
time_step = param.time_step

# Must have a function for this
p_delay_time = [0, -4, 0, -2,
                -4, -15, -2, -11,
                -22, -7]

top_color_label = param.top_color_label
bottom_color_label = param.bottom_color_label

# Extra data request in seconds to make sure we have enough data for all event(s).
request_padding = param.request_padding

frames_per_second = param.frames_per_second
insert_aux = False
selected_station = {'ref': {'set': False, 'net': str(), 'sta': str(), 'max': float, 'tr': None,
                            'times': list(), 'output_factor': float},
                    'aux': {'set': False, 'net': str(), 'sta': str(), 'max': float, 'tr': None,
                            'times': list(), 'output_factor': float}
                    }

codec = param.codec

contour_alpha = param.contour_alpha

travel_time_model = param.travel_time_model

# If not None, trace is normalized by dividing by specified value norm instead of dividing by its absolute maximum.
# If a negative value is specified then its absolute value is used. If it is zero (either through a zero array or
# by being passed), nothing will happen and the original array will not change.
norm = None

std_check = param.std_check
std_water_level = param.std_water_level
std_max = param.std_max
std_window = param.std_window
if not std_check:
    std_window = 0

# Width of a cosine taper to apply to the waveform data in time domain prior to deconvolution.
# Specified as a fraction of the trace length from 0 to 0.5.
do_taper = param.do_taper
taper_fraction = param.taper_fraction

# dc_to_exclude = ['NCEDC', 'SCEDC']
dc_to_exclude = param.dc_to_exclude

# do_event = sys.argv[1]
do_event = ''
# plot_type = sys.argv[2]
plot_type = param.plot_type
net_to_exclude = param.net_to_exclude

if not scale_plot:
    gain = 1.0
if plot_type == 'contour':
    plot_contours = True
    plot_markers = False
else:
    plot_contours = False
    plot_markers = True
if plot_contours:
    utils.print_message('ERR', 'Plot type contour is not supported under current version!', None)
    sys.exit(2)

show_station_locations = param.show_station_locations

if plot_markers:
    show_station_locations = False
if not plot_markers and not plot_contours:
    utils.print_message('ERR', 'Both markers and contours are turned OFF!', None)
    sys.exit(2)

filter_freqmin = param.filter_freqmin
filter_freqmax = param.filter_freqmax

output_factor = None

output_file = None

animation_start_time = None

event_marker = None
marker_lines = None
time_guide = None
time_guide_label = None
sta_type = None

low_water_level = None
high_water_level = None

request_max_lat = None
request_min_lon = None
request_max_lon = None
request_min_lat = None

title_text = None

reference_net = None
reference_sta = None

event_time = [None]
event_lat_lon = [None]
event_mag = [None]
event_depth = [None]
event_id = [None]

tt_url = None
event_lat = None
event_lon = None

contour_lat = [(25, 50), (52, 72)]
contour_lon = [(-125, -65), (-170, -130)]
# contour_lat = [(30, 50)]
# contour_lon = [(-129, -107)]
meridian_inc = 10
parallel_inc = 10

timer_is_on = param.timer_is_on

min_trace_length = param.min_trace_length

spatial_resolution = param.spatial_resolution
grid_method = param.grid_method

# List of seismic phases to calculate travel times for
phase_list = param.phase_list

# Number of seconds between two consecutive phases to avoid over plots
# Rayleigh slowness (s/degree); originally 29.835 but using slowness of 28.5 to provides a better phase display
r_slowness = param.r_slowness

# Iris Travel Time service URL.
iris_traveltime_url = param.iris_traveltime_url

location_order = param.location_order

exclude_temporary_networks = True

# The number of rgb quantization levels
color_levels = param.color_levels
sampling_rate = param.sampling_rate

# Define a filter band to prevent amplifying noise during the deconvolution.
pre_filter = param.pre_filter

verbose = param.verbose

source_track_file = param.source_track_file
track_records = list()

if source_track_file is not None:
    track_file = open(source_track_file, 'r')
    data = track_file.read()
    track_file.close()
    track_records = data.split('\n')


# ----------------- PARAMETERS END -----------------

def divide_to_chunks(long_list, chunks):
    """Yield successive chunk_count-sized chunks from long_list."""
    # For item i in a range that is a length of l,
    for i in range(0, len(long_list), chunks):
        # Create an index range for l of n items:
        yield long_list[i:i + chunks]


def chan_label(band_list, comp):
    """Create a label for band designation"""
    band_count = len(band_list)
    if band_count <= 0:
        utils.print_message('ERR', f'No bands given {band_list}!\n', log_file)
        return None
    elif band_count == 1:
        chan_text = f'band: {band_list[0]}'
    elif band_count == 2:
        chan_text = f'band: {band_list[0]} or {band_list[1]}'
    else:
        chan_text = ', '.join(band_list[:len(band_list) - 1])
        chan_text = f'band: {chan_text}, or {band_list[-1]}'

    if comp == '1':
        chan_text = f'{chan_text}; orientation: Z'
    else:
        chan_text = f'{chan_text}; orientation: 3-C'
    return chan_text


def select_unit_factor(data_list, unit_symb_dict):
    """Check a dictionary of unit symbols to find the best unit factor to use
       such that the maximum output value is < 1000
    """
    _max = max(abs(data_list))
    use_factor = None
    for _key in sorted(unit_symb_dict.keys()):
        if use_factor is None:
            use_factor = _key
        else:
            if _max * _key < 1000.0:
                use_factor = _key
    return use_factor


def sign(input_number, plus=True):
    """Returns sign of a number as a string/"""
    if float(input_number) < 0.0:
        return '-'
    else:
        if plus:
            return '+'
        else:
            return ''


def is_number(n):
    """Check if the input string input is a number.
    """
    try:
        float(n)
    except ValueError:
        return False
    return True


def zero_to_one(in_value):
    """Convert a given value between -1 and +1 to zero to +1 range.
    """
    if in_value is None:
        return None
    elif 1.0 < in_value < -1.0:
        utils.print_message('ERR', f'Bad value to index ({in_value}). Must be between -1 and +1!\n', log_file)
        return None

    # Normalize from -1/+1 to 0/+1
    norm_value = 0.5 * in_value + 0.5
    return norm_value


def time_it(in_t0, forced=False):
    """Compute the elapsed time since  last call (t0)"""
    if not forced:
        if not timer_is_on:
            return in_t0
    t1 = time.time()
    dt = t1 - in_t0
    if dt >= 60.0:
        utils.print_message('TIME', 'Elapsed time: {:0.1f}m'.format(dt / 60.0), log_file)
    else:
        utils.print_message('TIME', 'Elapsed time: {:0.2f}s'.format(dt), log_file)

    return time.time()


def set_marker(map_handle, use_color):
    """Create a map marker to plot on the base map"""
    this_marker = map_handle.plot([], [], marker=marker[technology], color=use_color,
                                  markersize=marker_size, markeredgewidth=0.05, markeredgecolor=(0, 0, 0, 1),
                                  linewidth=0, zorder=5)
    return this_marker


def value_to_color(use_cmap, this_value, scale, low_level=None, high_level=None):
    """Convert a given value  to a color value from selected color map.
    """
    if abs(this_value) <= zero_base:
        this_value = 0.0
    elif low_level is not None and this_value <= low_level:
        this_value *= scale

    elif high_level is not None and this_value >= high_level:
        this_value *= scale

    elif low_level is None and high_level is None:
        this_value *= scale

    if abs(this_value) <= zero_base:
        this_value = 0.0
    elif this_value < -1.0:
        this_value = -1.0
    elif this_value > 1.0:
        this_value = 1.0

    rgba_color = use_cmap(zero_to_one(this_value))
    return rgba_color


def get_request_items(req_line):
    """Split a request line to its components."""
    (req_net, req_sta, req_loc, req_chan, req_start, req_end) = req_line.strip().split()
    return req_net, req_sta, req_loc, req_chan, req_start, req_end


def get_service_url(ws_catalog, ws_dc):
    """Extract the service URL from dataselect service URL."""
    ws_service_url = ws_catalog[ws_dc]['dataselect_service'].split('/fdsn')[0]
    return ws_service_url


def get_dc(dc_url):
    """Get list of the FDSN data centers and the associate information."""

    if verbose:
        utils.print_message('INFO', f'requesting data center information from: {dc_url}', log_file)

    _dc_info = {}
    _dc_data = json.loads(utils.read_url(dc_url, log_file))

    for _dc in _dc_data:
        _dc_info[_dc['name']] = utils.ObjDict(_dc)

    return utils.ObjDict(_dc_info)


def is_net_temporary(net):
    """Exclude temporary networks."""
    if len(net) <= 2:
        if net[0].isdigit():
            return True
        if net[0].lower() in ['x', 'y', 'z']:
            return True
    return False


def get_chan_band(chan):
    """Extract band of a given channel code."""
    if len(chan) != 3:
        utils.print_message('ERR', f'Invalid channel code {chan}!\n', log_file)
        return None
    else:
        return chan[0:2]


def get_chan_comp(chan):
    """Extract component of a given channel code."""
    if len(chan) != 3:
        utils.print_message('ERR', f'Invalid channel code {chan}!\n', log_file)
        return None
    else:
        return chan[2]


def make_cmap():
    """Create and register a custom colormap.
    """

    cmap_name = 'bwr'  # 'bwr_r'#''coolwarm'
    cmap = None
    try:
        cmap = cm.get_cmap(cmap_name, lut=color_levels)
    except Exception as _er:
        utils.print_message('ERR', f'problem with registering the {cmap_name} colormap\n{_er}', log_file)

        sys, exit(1)
    cmap.set_bad(alpha=0.0)

    return cmap


def get_travel_times(req_url):
    """Get travel time for the selected phases or return HTML code if code=True"""
    req = Request(req_url)
    travel_time_records = dict()
    sorted_travel_time_records = OrderedDict()

    utils.print_message(f'INFO', f'requesting\n{req_url}', log_file)

    try:
        _response = urlopen(req)
    except HTTPError as _er:
        utils.print_message(f'ERR', f'response:\n{_er}', log_file)
        return travel_time_records
    except URLError as _er:
        utils.print_message(f'ERR', f'response:\n{_er}', log_file)
        return travel_time_records
    except Exception as _er:
        utils.print_message(f'ERR', f'response:\n{_er}', log_file)
        return travel_time_records

    mybytes = _response.read()
    _response.close()
    _lines = mybytes.decode('utf-8').split('\n')
    dist_degrees = None
    for _line in _lines:
        if not _line.strip():
            continue

        # Find the phases.
        _values = _line.strip().split()

        if _values:
            _phase = _values[2]
            dist_degrees = float(_values[0])
            _time = float(_values[3])

            if _phase not in travel_time_records.keys():
                travel_time_records[_phase] = _time
    travel_time_records['R1'] = dist_degrees * r_slowness
    travel_time_records['R2'] = (360.0 - dist_degrees) * r_slowness
    sorted_tt = sorted(travel_time_records.items(), key=operator.itemgetter(1))
    for _item in sorted_tt:
        sorted_travel_time_records[_item[0]] = _item[1]
    print(f'[INFO] TT: {sorted_travel_time_records}, delta:{dist_degrees}', file=log_file, flush=True)
    return sorted_travel_time_records, dist_degrees


def get_fedcatalog_stations(req_url, req_start, req_end, req_band, req_comp,
                            do_exclude_temporary_networks, run_verbose=False):
    """Get station list from fedcatalog service."""

    # This dictionary stores all the fedcatalog information.
    fedcatalog_info = dict()

    # This dictionary provides a template for fetdatalog creation.
    catalog_info = dict()

    bulk_list = collections.OrderedDict()
    dc_chunk_list = dict()

    # Since requests are made in the order of band priority, we create a done_list and
    # once a station is inserted in the done _list, it is complete.
    done_list = list()
    for band_index, _band in enumerate(req_band):

        # Build the channel list based on the band an the component list. Then send the request out.
        _chan = list()
        for this_comp in channel_comp_list[req_comp]:
            _chan.append(f'{_band}{this_comp}')
        _url = req_url.replace('START', req_start).replace('END', req_end).replace('CHAN', ','.join(_chan))
        utils.print_message('INFO', f'sending request to fedcatalog for {_band} band: {_url}', log_file)

        try:
            content = utils.read_url(_url, log_file)
        except Exception as _er:
            utils.print_message('ERR', 'Request  {}: {}'.format(_url, _er), log_file)
            continue

        # Go through the station list and see if they qualify.
        _lines = content.split('\n')

        _line_index = -1
        previous_dc = None
        dc_name = None
        for _line in _lines:
            _line_index += 1

            # Skip the blank and the comment lines.
            if not _line.strip() or _line.startswith('#'):
                continue

            # From the parameter=value lines, we are interested in the DATACENTER and DATASELECTSERVICE lines.
            elif '=' in _line:
                _par, _value = _line.split('=')

                # Found the data center name.
                if _par == 'DATACENTER':
                    if dc_name is not None:
                        previous_dc = dc_name
                    utils.print_message('INFO', f'from the {_value} data center', log_file)
                    dc_name, dc_url = _value.strip().split(',')

                    # Initialize the data center information, create chunk_count containers for chunked requets.
                    if dc_name not in catalog_info.keys():
                        utils.print_message('INFO', f'Initiating fedcatalog request for {dc_name}', log_file)
                        catalog_info[dc_name] = utils.ObjDict(
                            {'url': dc_url, 'dataselect_service': '', 'bulk': []})

                    # if this is not the first data center, save the previous data center's bulk list (only the
                    # complete list for the requested components (1C or 3C) are saved.
                    if len(bulk_list.keys()) > 0:
                        this_dc_list = list()
                        for _key in bulk_list.keys():
                            # 1-C list.
                            if req_comp == '1':
                                if bulk_list[_key][0]:
                                    done_list.append(_key)
                                    if run_verbose:
                                        utils.print_message('INFO', f'{_key} in {req_comp}C request list', log_file)
                                    this_dc_list.append(bulk_list[_key][0])
                            elif bulk_list[_key][0] and bulk_list[_key][1] and bulk_list[_key][2]:
                                # 3-C list
                                for i in range(3):
                                    if run_verbose:
                                        utils.print_message('INFO', f'{_key} in {req_comp}C request list', log_file)
                                    this_dc_list.append(bulk_list[_key][i])
                                done_list.append(_key)

                        # Break the  list into chunks and add it to fedcatalog_info. We incorporate band_index,
                        # in case multiple bands are requested. Otherwise, cunk_index of the next band will overwrite
                        # chunk_index of this band.
                        for chunk_index, chunk in enumerate(divide_to_chunks(this_dc_list, chunk_count)):
                            chunk_dc = f'{previous_dc}_{band_index}_{chunk_index}'

                            # Keep track of chunks for each DC for later use.
                            if previous_dc not in dc_chunk_list.keys():
                                dc_chunk_list[previous_dc] = list()
                            dc_chunk_list[previous_dc].append(chunk_dc)

                            fedcatalog_info[chunk_dc] = catalog_info[previous_dc].copy()
                            fedcatalog_info[chunk_dc]['bulk'] = chunk

                        # The list is saved. Now, reset the bulk_list.
                        bulk_list = collections.OrderedDict()

                    continue
                # Found the dataselect service address.
                elif _par == 'DATASELECTSERVICE':
                    # Save the dataselect service address for all chunks.
                    if dc_name in dc_chunk_list.keys():
                        for chunk_dc in dc_chunk_list[dc_name]:
                            fedcatalog_info[chunk_dc]['dataselect_service'] = _value.strip()

                    # Save the dataselect service address in the catalog for this DC,
                    catalog_info[dc_name]['dataselect_service'] = _value.strip()
                    utils.print_message('INFO', f'dataselect service is {_value.strip()}', log_file)
                    continue
                else:
                    # Ignore the other definitions.
                    continue

            # The rest are the station lines.
            if True:
                # Skip the blank lines.
                if not (_line.strip()):
                    continue

                # Get the station information.
                net, sta, loc, chan, sta_start, sta_end = get_request_items(_line)

                # See if network is acceptable.
                if net in net_to_exclude:
                    utils.print_message('WARN', f'Skipped {net}.{sta}.{chan}, {net} is in the network exclude list',
                                        log_file)
                    continue
                # Exclude the temporary networks.
                if do_exclude_temporary_networks and is_net_temporary(net):
                    utils.print_message('WARN', f'{net}.{sta}.{chan} rejected, the {net} temporary networks not valid ',
                                        log_file)
                    continue

                # This should not happen but still chack and only accept channels that are in the request_bands.
                if get_chan_band(chan) not in _band:
                    utils.print_message('WARN', f'{net}.{sta}.{chan} rejected, {chan} band not in {_band} ', log_file)
                    continue

                # If this net.sta combination has not been selected before, accept it.
                _net_sta_key = f'{net}_{sta}'
                if _net_sta_key in done_list:
                    utils.print_message('WARN', f'{net}.{sta}.{chan} rejected, station already in the request list ',
                                        log_file)
                    continue

                # Accept the first channel of a new station.
                if _net_sta_key not in bulk_list.keys():
                    bulk_list[_net_sta_key] = [[], [], []]
                    chan_index = channel_comp_index[get_chan_comp(chan)]
                    if run_verbose:
                        utils.print_message('INFO', f'channel {chan} from {_net_sta_key} in {req_comp}C request list',
                                            log_file)
                    bulk_list[_net_sta_key][chan_index] = (net, sta, loc, chan, req_start, req_end)
                # We already have a record for this sta.chan, see if the new location or channel has priority.
                else:
                    # First check on the band, if the same band, it must be a different component, insert it.
                    existing_loc = None
                    existing_band = None
                    for _i in range(3):
                        if bulk_list[_net_sta_key][_i]:
                            existing_band = get_chan_band(bulk_list[_net_sta_key][_i][3])
                            existing_loc = bulk_list[_net_sta_key][_i][2]
                            break
                    chan_index = channel_comp_index[get_chan_comp(chan)]
                    # If already got all the components, no need to consider the new channel/location.
                    if (len(bulk_list[_net_sta_key][0]) and req_comp == '1') or \
                            (len(bulk_list[_net_sta_key][0]) and len(bulk_list[_net_sta_key][1]) and
                             len(bulk_list[_net_sta_key][2]) and req_comp == '3'):
                        utils.print_message('INF', f'skipped {net}.{sta}.{loc}.{chan} because already got all '
                                                   f'the channels for {net}.{sta}.{existing_loc}', log_file)
                    # If it is a different channel from the same location, save it.
                    elif loc == existing_loc and not bulk_list[_net_sta_key][chan_index]:
                        if run_verbose:
                            utils.print_message('INFO', f'channel {chan} from '
                                                        f'{_net_sta_key} in {req_comp}C request list', log_file)
                        bulk_list[_net_sta_key][chan_index] = (net, sta, loc, chan, req_start, req_end)
                    elif loc != existing_loc:
                        utils.print_message('INFO', f'Not enough components, replacing the existing band '
                                                    f'{existing_band} with {get_chan_band(chan)} for {net}.{sta}',
                                            log_file)
                        bulk_list[_net_sta_key] = [[], [], []]
                        if run_verbose:
                            utils.print_message('INFO', f'channel {chan} from '
                                                        f'{_net_sta_key} in {req_comp}C request list', log_file)
                        bulk_list[_net_sta_key][chan_index] = (net, sta, loc, chan, req_start, req_end)
                        continue

        # Save the last data center's bulk list.
        if len(bulk_list.keys()) > 0:
            this_dc_list = list()
            for _key in bulk_list.keys():
                # 1-C list.
                if req_comp == '1':
                    if bulk_list[_key][0]:
                        done_list.append(_key)
                        if run_verbose:
                            utils.print_message('INFO', f'{_key} in {req_comp}C request list', log_file)
                        this_dc_list.append(bulk_list[_key][0])
                # 3-C list.
                elif bulk_list[_key][0] and bulk_list[_key][1] and bulk_list[_key][2]:
                    done_list.append(_key)
                    if run_verbose:
                        utils.print_message('INFO', f'{_key} in {req_comp}C request list', log_file)
                    for i in range(3):
                        this_dc_list.append(bulk_list[_key][i])

            # Break the  list into chunks and add it to fedcatalog_info.
            for chunk_index, chunk in enumerate(divide_to_chunks(this_dc_list, chunk_count)):
                chunk_dc = f'{dc_name}_{chunk_index}'

                # Keep track of chunks for each DC for later use.
                if dc_name not in dc_chunk_list.keys():
                    dc_chunk_list[dc_name] = list()
                dc_chunk_list[dc_name].append(chunk_dc)

                fedcatalog_info[chunk_dc] = catalog_info[dc_name].copy()
                fedcatalog_info[chunk_dc]['bulk'] = chunk

            # Reset the bulk_list.
            bulk_list = collections.OrderedDict()

    return utils.ObjDict(fedcatalog_info)


# ----- MAIN -----
q_set = False

try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hHvb:c:d:D:e:g:Gi:l:m:n:N:o:p:P:q:r:s:S:t:T:z:',
                                       ['help', 'header', 'verbose', 'band=', 'comp=', 'dur=', 'std=', 'eloc=', 'gain=',
                                        'gc', 'eid=', 'emag=', 'output=', 'phase=', 'qscale=', 'tstep=',
                                        'net=', 'rnet=', 'delay=', 'region=', 'sizem=', 'rsta=', 'etime=', 'title=',
                                        'depth='])
    for opt, arg in options:
        if opt in('-h', '--help'):
            usage()
            sys.exit(3)
        elif opt in ('-v', '--verbose'):
            verbose = True
        elif opt in ('-b', '--band'):
            request_bands = arg.replace(' ', '').strip().split(',')
        elif opt in ('-c', '--comp'):
            request_comp = arg
            if request_comp not in channel_comp_list.keys():
                print(f'\n[ERR] Invalid component, "{request_comp}"!')
                print(f'Valid components are: {channel_comp_list.keys()}')
                for key, value in proj_regions.items():
                    print(f'{key}\u2014{value["name"]}')
                sys.exit(3)
        elif opt in ('-d', '--dur'):
            animation_duration = float(arg)
        elif opt in ('-D', '--std'):
            std_value = float(arg)
            if std_value < 0:
                std_check = False
                std_window = 0
            else:
                std_check = True
                std_max = std_value
        elif opt in ('-e', '--eloc'):
            event_lat_lon = arg.replace('[', '').replace(']', '').replace(', ', ',').split()
            if len(event_lat_lon) <= 0 or ',' not in event_lat_lon[0]:
                usage()
                utils.print_message('ERR', f'Bad event location, -e or --eloc {arg}', None)
                exit(3)

        elif opt in ('-i', '--eid'):
            event_id = arg.replace('[', '').replace(']', '').replace(', ', ',').split()
        elif opt in ('-m', '--emag'):
            event_mag = arg.replace('[', '').replace(']', '').replace(', ', ',').split()
        elif opt in ('-n', '--net'):
            request_net = arg
            if request_net.lower() == 'all':
                request_net = '*'
        elif opt in ('-l', '--tstep'):
            time_step = float(arg.strip())
        elif opt in ('-N', '--rnet'):
            reference_net = arg
        elif opt in ('-o', '--output'):
            output_file = arg.strip()
        elif opt in ('-p', '--delay'):
            animation_delay = float(arg)
        elif opt in ('-P', '--phase'):
            phase_spacing = float(arg)
        elif opt in ('-r', '--region'):
            region = arg.strip().lower()
            if region not in proj_regions.keys():
                print(f'\n[ERR] Invalid region name, "{region}"!')
                print('Valid regions are:')
                for key, value in proj_regions.items():
                    print(f'{key}\u2014{value["name"]}')
                sys.exit(4)
            figure_aspect_ratio = float(proj_regions[region]['size'][0] / proj_regions[region]['size'][1])
            if region == 'gl':
                meridian_inc = 60
                parallel_inc = 20
                logo_zoom *= 0.7
            elif region == 'np':
                meridian_inc = 20
                parallel_inc = 20
                logo_zoom *= 0.7
            marker_size = proj_regions[region]['marker_size']
        elif opt in ('-s', '--sizem'):
            marker_size = float(arg)
        elif opt in ('-q', '--qscale'):
            quiver_scale = float(arg)
            q_set = True
        elif opt in ('-g', '--gain'):
            gain = float(arg)
            scale_plot = True
        elif opt in ('-G', '--gc'):
            draw_great_circle = True
        elif opt in ('-S', '--rsta'):
            reference_sta = arg
        elif opt in ('-t', '--etime'):
            event_time = arg.replace('[', '').replace(']', '').replace(', ', ',').split()
        elif opt in ('-T', '--title'):
            title_text = arg
        elif opt in ('-z', '--depth'):
            event_depth = arg.replace('[', '').replace(']', '').replace(', ', ',').split()
except getopt.GetoptError as er:
    usage()
    utils.print_message('[ERR]', er, None)
    sys.exit(3)

log_file = sys.stdout
log_to_screen = param.log_to_screen
if not log_to_screen:
    if log_dir is None:
        log_file = open(sys.stdout, 'a')
        utils.print_message('WARN', f'Output to screen', log_file)
    else:
        log_file_name = os.path.join(param.log_dir, script.replace('.py', ''))
        log_file_name = f"{log_file_name}_{datetime.datetime.now().strftime('%Y-%m-%d')}"
        log_file_name = '.'.join([log_file_name, region, request_comp, 'log'])
        log_file = open(log_file_name, 'a')
        sys.stdout = log_file
print(f"\n\n............{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} START............", file=log_file,
      flush=True)

if event_lat_lon[0] is None:
    usage()
    utils.print_message('ERR', f'Event location (-e, --eloc) not defined', log_file)
    exit(3)

if event_depth[0] is None:
    usage()
    utils.print_message('ERR', f'Event Depth (-z, --depth) not defined', log_file)
    exit(3)

default_title = False
if title_text is None:
    default_title = True
    title_text = f'M {event_mag[0]} Event located at {event_lat_lon[0]}'
if reference_net is None or reference_sta is None:
    reference_net, reference_sta = proj_regions[region]['ref_sta']
if not q_set:
    quiver_scale = proj_regions[region]['quiver_scale']
    quiver_scale_units = proj_regions[region]['quiver_scale_units']

utils.print_message('INFO', f'Reference station {reference_sta} from network {reference_net}', log_file)

chan_list = list()
for band in request_bands:
    for _comp in channel_comp_list[request_comp]:
        chan_list.append(f'{band}{_comp}')
request_chans = ','.join(chan_list)

utils.print_message('INFO', f'Running for region {region}.\n', log_file)

(request_min_lat, request_max_lat), (request_min_lon, request_max_lon) = proj_regions[region]['limits']
region_name = proj_regions[region]['name']
(llcrnrlat, llcrnrlon), (urcrnrlat, urcrnrlon) = proj_regions[region]['corners']
(lat_1, lat_2), (lon_1, lon_2) = proj_regions[region]['limits']
(lat_0, lon_0) = proj_regions[region]['center']

metadata.append(f"movie_lat_min={llcrnrlat}")
metadata.append(f"movie_lat_max={urcrnrlat}")
metadata.append(f"movie_lon_min={llcrnrlon}")
metadata.append(f"movie_lon_max={urcrnrlon}")

figure_height = float(figure_width / figure_aspect_ratio)
# Increase the figure height to allow room for the trace based on the total/map rows ratio.
if trace_height is None:
    map_rows_ratio = float(subplot_rows / subplot_map_rows)
    figure_height *= map_rows_ratio
# Increase the figure height to allow room for the trace based on the given trace height.
else:
    figure_height += trace_height
    subplot_trace_rows = round(subplot_rows * float(trace_height / figure_height))
    subplot_map_rows = subplot_rows - subplot_trace_rows

# This is a test query! In the actual query we must be able to ask for multiple channels and internally use a
# channel priority system to select channels and discard the rest so we do not spend time to do pre-processing.
url_fedcatalog_request = f'http://service.iris.edu/irisws/fedcatalog/1/query?net={request_net}&sta=*&loc=*&' \
    f'cha=CHAN&targetservice=dataselect&level=channel&format=request&' \
    f'startbefore=START&endafter=END&' \
    f'maxlat={request_max_lat}&minlon={request_min_lon}&maxlon={request_max_lon}&' \
    f'minlat={request_min_lat}&' \
    f'includeoverlaps=false&nodata=404'

dc_info = get_dc(url_data_centers)

eq_size = list()

utils.print_message('INFO', f'Magnitude list {event_mag}.\n', log_file)
try:
    for mag in event_mag:
        if not is_number(mag):
            eq_size.append(mag)
        elif float(mag) >= 6.5:
            eq_size.append('large')
        elif 6.5 > float(mag) >= 6.0:
            eq_size.append('medium')
        else:
            eq_size.append('small')
except Exception as ex:
    usage()
    utils.print_message('ERR', f'A valid magnitude (-m or --emag) is required\n{ex}', None)
    exit(3)


event_timestamp = list()
animation_timestamp = list()
reftime_delays = list()

for e_index, e_time in enumerate(event_time):
    try:
        try:
            e_time = e_time.replace(' ', 'T')
            t0 = UTCDateTime(event_time[0].replace(' ', 'T'))
            event_timestamp.append(UTCDateTime(e_time))
        except Exception as ex:
            usage()
            utils.print_message('ERR',
                                f'Invalid event_time (-t or -etime options) {e_time}\n'
                                f'{ex}', log_file)
            sys.exit(3)

        reftime_delays.append(event_timestamp[-1] - t0 + p_delay_time[e_index])
        animation_timestamp.append(event_timestamp[-1] + animation_delay)
    except Exception as ex:
        utils.print_message('ERR', f'Could not convert event_time {e_time} with animation delay of {animation_delay}\n'
                                   f'{ex}', log_file)
        sys.exit(3)

request_start_time = list()
request_end_time = list()
animation_start_time = list()
animation_end_time = list()
for a_index, a_time in enumerate(animation_timestamp):
    animation_start_time.append(a_time.strftime('%Y-%m-%dT%H:%M:%S'))
    animation_end_time.append(a_time + animation_duration - p_delay_time[a_index])

    # Always request std_window seconds before the animation time for STD check.
    request_start_time.append((a_time - std_window).strftime('%Y-%m-%dT%H:%M:%S'))
    utils.print_message('INFO', f'Adjusting data request start time from {animation_start_time[-1]} to '
                                f'{request_start_time[-1]} for STD check\n', log_file)

    """
    if abs(animation_delay) < std_window and std_check:
        request_start_time.append((a_time - std_window).strftime('%Y-%m-%dT%H:%M:%S'))
        utils.print_message('INFO', f'Adjusting data request start time from {animation_start_time[-1]} to '
                                    f'{request_start_time[-1]} for STD check\n', log_file)
    else:
        request_start_time.append(animation_start_time[-1])
    """
    request_end_time.append((animation_end_time[-1] + request_padding).strftime('%Y-%m-%dT%H:%M:%S'))

frame_count = int(animation_duration / time_step) + 1

# Manoch 2019-09-16 This assumes that events are sorted, make it more general
catalog_list = list()
catalog_times = list()
for index, t_value in enumerate(request_start_time):
    cat = get_fedcatalog_stations(url_fedcatalog_request, request_start_time[index], request_end_time[index],
                                  request_bands, request_comp, exclude_temporary_networks, run_verbose=verbose)
    if cat is not None:
        catalog_list.append(cat)
        catalog_times.append((request_start_time[index], request_end_time[index]))

if not catalog_list:
    sys.exit(2)

# Set the default map parameters.
url_data_centers = param.url_data_centers
url_fedcatalog = param.url_fedcatalog

colors = param.colors

network = param.network
station = param.station

basemap_type = param.basemap_type

line_width = param.line_width
color = param.color

dpi = param.dpi
frames_per_second = param.frames_per_second

# technology = param.technology
channels = param.channels
marker = param.marker
frame_interval = param.frame_interval
author = param.author
bit_rate = param.bit_rate

current_cmap = make_cmap()

start = param.start
end = param.end
include_restricted = param.include_restricted

title_style = param.title_style
title_font_size = param.title_font_size
title_color = param.title_color
bbox_alpha = param.bbox_alpha
bbox_padding = param.bbox_padding

fig = plt.figure(figsize=(figure_width, figure_height), dpi=dpi)

utils.print_message('INFO', f'Figure Width: {figure_width}, Height: {figure_height}, DPI: {dpi}', log_file)

# Trace subplot with shape(r,c), location(r,c), rowspan, colspan. Here we reduce the trace column length by 1 to
# leave space for the legend.
ax2 = plt.subplot2grid((subplot_map_rows + subplot_trace_rows, subplot_columns), (subplot_map_rows, 0),
                       rowspan=subplot_trace_rows, colspan=subplot_columns - 1)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_position(('data', 0))
ax2.spines['bottom'].set_linewidth(0.3)
ax2.spines['left'].set_visible(False)
ax2.set_xticklabels([])
ax2.ticklabel_format(axis='y', style='sci')
ax2.tick_params(right=False, top=False, left=True, bottom=False)
ax2.yaxis.set_tick_params(labelsize=4)

# GMV subplot with shape(r,c), location(r,c), rowspan, colspan.
ax1 = plt.subplot2grid((subplot_map_rows + subplot_trace_rows, subplot_columns), (0, 0),
                       rowspan=subplot_map_rows, colspan=6)

# To make sure the images are tight and space for y-labels are not included.
y_label = ax1.set_ylabel('')
y_label.set_in_layout(False)

if default_title:
    title_string = title_text
else:
    title_string = title_text.title()

if len(event_time) != 1:
    plt.title(title_string)
else:
    plt.title(f'{title_string}\nOrigin Time (OT) = {event_time[0].replace(" ", "T").split("T")[1]} UTC')

# Place the tight_layout after most of the elements are in, so the layout can be configured properly.
# Info - pad: Padding between the figure edge and the edges of subplots, as a fraction of the font size.
# Info - h_pad, w_pad : Padding (height/width) between edges of adjacent subplots, as a fraction of the font size.
#        Defaults to pad.
plt.tight_layout(pad=2, h_pad=0, w_pad=0, rect=None)

# Must specify lat/lon values of corners (llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat) in degrees.
# rsphere=(6378137.00,6356752.3142) specifies WGS84 ellipsoid
# area_thresh=1000 means don't plot coastline features less
# than 1000 km^2 in area.
# lon_0,lat_0 is central point.
proj_type = proj_regions[region]['projection']
if proj_type == 'lcc':
    bg_map = Basemap(projection=proj_type, resolution='l',
                     rsphere=(6378137.00, 6356752.3142), area_thresh=1000.0,
                     llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                     lat_1=lat_1, lat_2=lat_2, lon_1=lon_1, lon_2=lon_2, lon_0=lon_0, lat_0=lat_0)
elif proj_type == 'robin':
    bg_map = Basemap(projection=proj_type, resolution='l',
                     area_thresh=1000.0, lat_0=proj_regions[region]['center'][0],
                     lon_0=proj_regions[region]['center'][1])
elif proj_type == 'npaeqd':
    bg_map = Basemap(projection=proj_type, resolution='l',
                     boundinglat=10, lon_0=270)
    """
    if region == 'am':
        bg_map = Basemap(llcrnrlon=-170., llcrnrlat=-60., urcrnrlon=-30., urcrnrlat=80.,
                         projection='lcc', lat_1=-40., lat_2=45., lon_0=-60., width=11798593.783438614,
                         height=16777242.309268944, resolution='l', area_thresh=1000.)
    else:
        bg_map = Basemap(projection=projection, resolution='l',
                         rsphere=(6378137.00, 6356752.3142), area_thresh=1000.0,
                         width=width, height=height,
                         lat_0=lat_0, lon_0=lon_0)  # lat_1=lat_1, lat_2=lat_2,
    # llcrnrlon=-174.0, llcrnrlat=25.0, urcrnrlon=-65.0, urcrnrlat=73.5, lon_0=-113.0, lat_0=50.0)
    """
else:
    bg_map = Basemap(projection=proj_type, resolution='l', area_thresh=1000.0, lat_0=0, lon_0=0)

if basemap_type == 'etopo':
    bg_map.etopo()
elif basemap_type == 'releaf':
    bg_map.shadedrelief()
else:
    bg_map.fillcontinents(color=color['land'], lake_color=color['lake'])
    # draw a boundary around the map, fill the background.
    # this background will end up being the ocean color, since
    # the continents will be drawn on top.
    bg_map.drawmapboundary(fill_color=color['lake'])

# Draw meridians and parallels (labels = [left,right,top,bottom]).
bg_map.drawmeridians(np.arange(0, 360, meridian_inc), linewidth=0.1, labels=[True, False, False, True], fontsize=8,
                     zorder=1000)
bg_map.drawparallels(np.arange(-90, 90, parallel_inc), linewidth=0.1, labels=[False, True, False, False], fontsize=8,
                     zorder=1000)

bg_map.drawcoastlines(linewidth=line_width['coastlines'])
bg_map.drawcountries(linewidth=line_width['countries'])
bg_map.drawstates(linewidth=line_width['states'])

# Mark the event locations.
for loc_index, loc_value in enumerate(event_lat_lon):
    event_loc = loc_value.strip().replace('(', '').replace(')', '').replace('', '').split(',')
    print(f'[INFO] Event location: {event_loc}', file=log_file, flush=True)
    event_text = 'event location'
    if len(event_lat_lon) > 1:
        event_text = f'{len(event_lat_lon)} event locationss'
    event_lat, event_lon = event_loc
    event_x, event_y = bg_map(event_lon, event_lat)
    if event_x > bg_map.xmax or event_x < bg_map.xmin or event_y < bg_map.ymin or event_y > bg_map.ymax:
        utils.print_message('INFO', 'Event location outside the view port\n', log_file)
        if loc_index == 0:
            bg_map.plot([event_x], [event_y], marker=marker_event, color='white', linestyle='None', markeredgewidth=0.0,
                        markersize=10, zorder=500, label=f'{event_text}\nnot in the viewport')
        else:
            bg_map.plot([event_x], [event_y], marker=marker_event, color='white', linestyle='None', markeredgewidth=0.0,
                        markersize=10, zorder=500)
    else:
        utils.print_message('INFO', '\nMarking the event location\n', log_file)
        if loc_index == 0:
            bg_map.plot([event_x], [event_y], marker=marker_event, color='y', linestyle='None', markeredgecolor='k',
                        markeredgewidth=0.5,
                        markersize=10, zorder=500, label=f'{event_text}')
        else:
            bg_map.plot([event_x], [event_y], marker=marker_event, color='y', linestyle='None', markeredgecolor='k',
                        markeredgewidth=0.5,
                        markersize=10, zorder=500)

# Add logo to the base map
if add_log:
    logo_image = np.array(Image.open(logo_file))
    im = OffsetImage(logo_image, zoom=logo_zoom)
    b_box = AnnotationBbox(im, logo_location, box_alignment=(-0.2, -0.2), frameon=False)

    # Get the axes object from the base map and add the AnnotationBbox artist
    bg_map._check_ax().add_artist(b_box)

# The map limits (corner points of the map).
x1, x2, y1, y2 = bg_map.xmin, bg_map.xmax, bg_map.ymin, bg_map.ymax

# Set the position for the title that we use to display time slider's information.
if region not in ['gl', 'np']:
    title_x, title_y = x2 * 0.99, y1 * 1.01
    title = ax1.text(title_x, title_y, ' ', verticalalignment='bottom', horizontalalignment='right', style=title_style,
                     color=title_color, size=title_font_size, zorder=100000,
                     bbox={'alpha': bbox_alpha, 'pad': bbox_padding})
else:
    title_x, title_y = (x1 + x2) / 2.0, y1 * 1.01
    title = ax1.text(title_x, title_y, ' ', verticalalignment='bottom', horizontalalignment='center', style=title_style,
                     color=title_color, size=title_font_size, zorder=100000,
                     bbox={'alpha': bbox_alpha, 'pad': bbox_padding})

first_time = UTCDateTime(animation_start_time[0]) - time_step
last_time = UTCDateTime(animation_end_time[0]) - time_step

frames = list()
frame_lines = list()

station_coordinates = collections.OrderedDict()
full_color_list = list()

t0 = time.time()
contour_data = list()

first = False
# Get data from each Data Center.
for cat_index, catalog in enumerate(catalog_list):
    start, end = catalog_times[cat_index]
    for dc in catalog.keys():
        st = None
        if dc in dc_to_exclude:
            utils.print_message('WARN', 'skipped data from {} because it is in the exclude list'.format(dc), log_file)
            continue
        else:
            if verbose:
                print('\n\n[INFO] Sending requests for:', file=log_file, flush=True)
                for line in catalog[dc]['bulk']:
                    print(line, file=log_file, flush=True)
                print('\n', file=log_file, flush=True)
            print('\n', file=log_file, flush=True)
            if len(catalog[dc]['bulk']) <= 0:
                utils.print_message('WARN', f'Skipping DC {dc}, no stations to request!\n', log_file)
                continue

            utils.print_message('INFO', 'Requesting data from {} via {}\n'.format(dc, get_service_url(catalog, dc)),
                                log_file)

        # Set the client up for this data center.
        try:
            client = Client(get_service_url(catalog, dc))
            st = client.get_waveforms_bulk(catalog[dc]['bulk'], attach_response=True)

            # For technology = 'seismic', we want to process Z channels first since that is the one we use to reject
            # or accept channel.
            if technology == 'seismic':
                st.sort(keys=['channel'], reverse=True)
        except Exception as er:
            utils.print_message('ERR', f'Request:\n\n{catalog[dc]["bulk"]}\n\n from {dc} '
                                       f'({get_service_url(catalog, dc)}) failed: {er}', log_file)
            t0 = time_it(t0)
            continue
        t0 = time_it(t0)

        rejection_list = list()
        # Work on individual traces in the stream.
        for index, tr in enumerate(st):

            # The trace time axis is based on the animation time. So, we get the trace times accordingly.
            trace_times = list(tr.times(reftime=animation_timestamp[cat_index]) - reftime_delays[cat_index])

            net_sta_key = '{}.{}'.format(tr.stats.network, tr.stats.station)

            # NOTE: Revise later for a better check.
            if trace_times[-1] - trace_times[0] < min_trace_length:
                utils.print_message('INFO', 'Skipped, channel {} of Net.Sta {} from {} because it is shorter than {} s'
                                    .format(tr.stats.channel, net_sta_key, dc, min_trace_length), log_file)
                continue

            # Get the station coordinates. Here we assume there is a possibility that station contains gaps, so
            # we may get more than one trace. So, we will get the station information from the first segment.
            # NOTE: Revisit this since we are not merging!
            if net_sta_key in station_coordinates.keys():
                utils.print_message('WARN', 'Multiple trace, already have data from {} for channel {}.'.
                                    format(net_sta_key, station_coordinates[net_sta_key][2]), log_file)
            else:
                if verbose:
                    utils.print_message('INFO', 'Getting information for station {}.{}.{}.{} from {} to {}.'.
                                        format(tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel,
                                               start, end), log_file)

            try:
                inventory = client.get_stations(network=tr.stats.network, station=tr.stats.station,
                                                location=tr.stats.location, channel=tr.stats.channel, starttime=start,
                                                endtime=end, level="station")
                station_coordinates[net_sta_key] = (inventory.networks[0].stations[0].longitude,
                                                    inventory.networks[0].stations[0].latitude,
                                                    tr.stats.channel)
            except Exception as er:
                utils.print_message('ERR', 'Request error {}.{} failed: {}'.format(dc, net_sta_key, er), log_file)
                continue

            t0 = time_it(t0)
            # Remove the response.
            if output_type in (None, 'PA', 'C'):
                utils.print_message('WARN', 'Response removal option is off for technology {}, removing trend and '
                                            'correcting for sensitivity ony!'.
                                    format(technology), log_file)
                tr.detrend()
                tr.data = tr.data / float(tr.stats.response.instrument_sensitivity.value)
                output_factor = 1.0
            else:
                try:
                    if verbose:
                        utils.print_message('INFO', 'Removing response from {}.{}'.format(dc, net_sta_key), log_file)
                    tr.remove_response(output=output_type, zero_mean=True, taper=do_taper,
                                       taper_fraction=taper_fraction, pre_filt=pre_filter)
                except ValueError as er:
                    utils.print_message('ERR', 'Removing response from {}.{} failed: {}'.format(dc, net_sta_key, er),
                                        log_file)
                    continue
                except Exception as er:
                    utils.print_message('ERR', 'Removing response from {}.{} failed: {}'.format(dc, net_sta_key, er),
                                        log_file)
                    continue
                t0 = time_it(t0)

            # Detrending.
            # try:
            #    utils.print_message('INFO', 'Detrending {}.{}'.format(dc, net_sta_key))
            #    tr.detrend()
            # except Exception as er:
            #    utils.print_message('ERR', 'Detrending {}.{} failed: {}'.format(dc, net_sta_key, er))
            #    continue
            # t0 = time_it(t0)

            # Filtering
            tr_no_filter = tr.copy()
            if eq_size[cat_index] in filter_freqmin.keys():
                try:
                    if verbose:
                        utils.print_message('INFO', 'Filtering {}.{} between {} and {}'.format(dc, net_sta_key,
                                                                                               filter_freqmin[
                                                                                                   eq_size[cat_index]],
                                                                                               filter_freqmax[
                                                                                                   eq_size[cat_index]]),
                                            log_file)
                    if not metadata_filter:
                        metadata.append(f"movie_filter_tlo={1.0 / filter_freqmin[eq_size[cat_index]]}")
                        metadata.append(f"movie_filter_thi={1.0 / filter_freqmax[eq_size[cat_index]]}")
                        metadata_filter = True

                    tr.filter('bandpass', freqmin=filter_freqmin[eq_size[cat_index]],
                              freqmax=filter_freqmax[eq_size[cat_index]], zerophase=True)
                except Exception as er:
                    utils.print_message('ERR', 'Filtering {}.{} failed: {}'.format(dc, net_sta_key, er), log_file)
                    continue
                t0 = time_it(t0)

            # We keep the reference trace in its true amplitude for plotting. For multiple event case, do it only if
            # ref_max is not defined.
            if tr.stats.station == reference_sta and tr.stats.network == reference_net and \
                    not selected_station['ref']['set']:
                sta_type = 'ref'
                selected_station[sta_type]['set'] = True
                selected_station[sta_type]['net'] = tr.stats.network
                selected_station[sta_type]['sta'] = tr.stats.station
                selected_station[sta_type]['tr'] = tr.copy()
                selected_station[sta_type]['tr'].data = selected_station[sta_type]['tr'].data
                selected_station[sta_type]['tr_no_filter'] = tr_no_filter.copy()
                selected_station[sta_type]['times'] = trace_times.copy()
                selected_station[sta_type]['output_factor'] = select_unit_factor(
                    selected_station[sta_type]['tr'].data, output_factor_symbol)
                selected_station[sta_type]['max'] = max(max(
                    selected_station[sta_type]['tr'].data),
                    max(selected_station[sta_type]['tr_no_filter'].data)) * selected_station[sta_type]['output_factor']
                utils.print_message('INFO', f'Normalizing {dc}.{net_sta_key}  as the reference station', log_file)
            # Keep the first trace as auxiliary trace in case we do not find the requested reference station.
            elif not selected_station['aux']['set']:
                sta_type = 'aux'
                selected_station[sta_type]['set'] = True
                selected_station[sta_type]['net'] = tr.stats.network
                selected_station[sta_type]['sta'] = tr.stats.station
                selected_station[sta_type]['tr'] = tr.copy()
                selected_station[sta_type]['tr'].data = selected_station[sta_type]['tr'].data
                selected_station[sta_type]['tr_no_filter'] = tr_no_filter.copy()
                selected_station[sta_type]['times'] = trace_times.copy()
                selected_station[sta_type]['output_factor'] = select_unit_factor(
                    selected_station[sta_type]['tr'].data, output_factor_symbol)
                selected_station[sta_type]['max'] = max(max(selected_station[sta_type]['tr_no_filter'].data), max(
                    selected_station[sta_type]['tr'].data)) * selected_station[sta_type]['output_factor']
                utils.print_message('INFO', f'Normalizing {dc}.{net_sta_key}  as the auxiliary station', log_file)

            # Traces are normalized to their absolute maximum not to the global maximum.
            try:
                if verbose:
                    utils.print_message('INFO', f'Normalizing {dc}.{net_sta_key} to its absolute maximum ', log_file)
                tr.normalize(norm)
            except Exception as er:
                utils.print_message('ERR', f'Normalization for {dc}.{net_sta_key} failed: {er}', log_file)
                continue
            t0 = time_it(t0)

            # If sta was previously rejected based on its STD, do not work on the other channels of this station.
            if net_sta_key in rejection_list:
                continue
            # Use the standard deviation of the first std_window samples to decide if trace is acceptable.
            # The reference station is exempt from this check.
            if std_check:
                if tr.stats.station == selected_station['ref']['sta'] and \
                        tr.stats.network == selected_station['ref']['net']:
                    utils.print_message('INFO', f'NOTE: The reference station {dc}.{net_sta_key} '
                                                f'is exempt from STD check! ', log_file)
                elif tr.stats.station == selected_station['aux']['sta'] and \
                        tr.stats.network == selected_station['aux']['net']:
                    utils.print_message('INFO', f'NOTE: The auxiliary station {dc}.{net_sta_key} '
                                                f'is exempt from STD check! ', log_file)
                elif net_sta_key in rejection_list:
                    continue
                elif channel_comp_index[get_chan_comp(tr.stats.channel)] > 0:
                    utils.print_message('INFO', f'{tr.stats.channel} is '
                                                f'exempt from STD check! ', log_file)
                else:
                    # This distance calculation may be useful for ref sta selection
                    # dist, azim, back_azim = gps2dist_azimuth(float(event_lat), float(event_lon),
                    #                                         float(station_coordinates[net_sta_key][1]),
                    #                                         float(station_coordinates[net_sta_key][0]))
                    # dist /= 1000.0
                    # if dist > std_far_distance_km:
                    #    std_max = std_far_max
                    # else:
                    #    std_max = std_near_max
                    # Perform STD on the channel of zero order (the first channel, or Z for seismic)
                    std_data = tr.data[0:int(std_window * tr.stats.sampling_rate)]

                    # For large events, STD may not be a good indicative of blinking. So, we only use STD criteria if
                    # the maximum amplitude is above the std_water_level.
                    abs_max_std = abs(max(std_data))
                    if abs_max_std > std_water_level or std_water_level <= 0:
                        std = np.std(np.array(std_data))
                        if std > std_max:
                            utils.print_message(f'WARN',
                                                f'Rejected {dc}.{net_sta_key} due to high STD of '
                                                f'{std:0.4f} > {std_max} (abs max STD data {abs_max_std} '
                                                f'std_water_level {std_water_level})', log_file)
                            rejection_list.append(net_sta_key)
                            continue

            utils.print_message('INFO', 'Processed channel {} of {} from {}'.format(tr.stats.channel,
                                                                                    net_sta_key, dc), log_file)
            # NOTES: 1. Later we want to use the raw waveform so we can post the trace height. For now working
            # with the normalized trace is better. For now we are only using BMO station!!
            # 2. Grab the first station in case we cannot find the requested station later.
            if selected_station['ref']['set']:
                sta_type = 'ref'
            else:
                sta_type = 'aux'

            if 'net_sta_key' not in selected_station[sta_type].keys():
                selected_station[sta_type]['net_sta_key'] = net_sta_key
                if verbose:
                    utils.print_message(f'INFO',
                                        f'** Selecting the {sta_type} trace {selected_station[sta_type]["net"]}.'
                                        f'{selected_station[sta_type]["sta"]}\n', log_file)
                if technology == 'seismic':
                    tt_url = iris_traveltime_url
                    tt_url = tt_url.replace('DEPTH', str(event_depth[cat_index]))

                # Mark the reference station location
                if verbose:
                    utils.print_message('INFO', f'** Getting the {start}. station location', log_file)
                selected_station[sta_type]['lon'] = station_coordinates[selected_station[sta_type]['net_sta_key']][0]
                selected_station[sta_type]['lat'] = station_coordinates[selected_station[sta_type]['net_sta_key']][1]

            # Slice the trace at individual frame intervals.
            frame_count = -1
            slice_time = animation_timestamp[cat_index] - reftime_delays[cat_index] - time_step
            skipped = False
            while slice_time <= last_time:

                # Keep track of the frame where the sample belongs.
                frame_count += 1

                # Set up  containers for each frame.
                if frame_count >= len(frames):

                    # Container for the station markers.
                    frames.append(collections.OrderedDict())
                    frame_lines.append(collections.OrderedDict())

                    # Contour  container for each contour area.
                    if plot_contours:
                        contour_data.append(list())
                        for no_mask_index in range(len(contour_lon)):
                            contour_data[frame_count].append({'lon': list(), 'lat': list(), 'x': list(), 'y': list(),
                                                              'value': list()})
                            contour_handle.append(None)
                slice_time += time_step

                # Skip if no data.
                if not len(tr.data):
                    utils.print_message('WARN', 'skipped empty trace from {} for {}'.format(dc, net_sta_key), log_file)
                    break

                # The first sample of the slice is what we want. We need for a slice to
                # have a width (2 * sampling rate) so we capture sufficient samples.
                tr2 = tr.slice(starttime=slice_time + reftime_delays[cat_index],
                               endtime=slice_time + reftime_delays[cat_index] + 2 * tr.stats.sampling_rate,
                               nearest_sample=True)

                # Skip if no data.
                if not len(tr2.data):
                    if not skipped:
                        skipped = True
                        utils.print_message('WARN', 'skipped empty slice(s) from {} for {}'.format(dc, net_sta_key),
                                            log_file)
                    continue

                # Convert station lon, lat to x, y.
                x, y = bg_map(station_coordinates[net_sta_key][0], station_coordinates[net_sta_key][1])

                # Keep the reference_sta color
                if tr2.stats.station == reference_sta and tr2.stats.network == reference_net:
                    this_color = value_to_color(current_cmap, tr2.data[0], gain,
                                                low_level=low_water_level, high_level=high_water_level)

                # If the station is within the selected contour area, then add it to the list.
                if plot_contours:
                    for no_mask_index in range(len(contour_lon)):
                        if min(contour_lat[no_mask_index]) <= station_coordinates[net_sta_key][1] <= \
                                max(contour_lat[no_mask_index]) \
                                and min(contour_lon[no_mask_index]) <= station_coordinates[net_sta_key][0] <= \
                                max(contour_lon[no_mask_index]):
                            contour_data[frame_count][no_mask_index]['lon'].append(station_coordinates[net_sta_key][0])
                            contour_data[frame_count][no_mask_index]['lat'].append(station_coordinates[net_sta_key][1])
                            contour_data[frame_count][no_mask_index]['value'].append(tr2.data[0].copy())
                            contour_data[frame_count][no_mask_index]['x'].append(x)
                            contour_data[frame_count][no_mask_index]['y'].append(y)

                if plot_markers:
                    comp_index = int(channel_comp_index[get_chan_comp(tr2.stats.channel)])
                    if comp_index == 0:
                        # Group station markers based on their color.
                        this_color = value_to_color(current_cmap, tr2.data[0], gain,
                                                    low_level=low_water_level, high_level=high_water_level)
                        # We will create marker plot handle for all colors, so create a list.
                        if this_color not in full_color_list:
                            full_color_list.append(this_color)

                        this_color_index = full_color_list.index(this_color)

                        if this_color_index not in frames[frame_count].keys():
                            frames[frame_count][this_color_index] = \
                                {'x': list(), 'y': list(), 'z': list(), 't': list(), 's': list()}

                        frames[frame_count][this_color_index]['x'].append(x)
                        frames[frame_count][this_color_index]['y'].append(y)
                        frames[frame_count][this_color_index]['z'].append(tr2.data[0])
                        frames[frame_count][this_color_index]['t'].append(slice_time.strftime('%Y-%m-%d %H:%M:%S UTC'))
                        frames[frame_count][this_color_index]['s'].append(frame_count * time_step)
                    else:
                        if net_sta_key in rejection_list:
                            continue
                        if net_sta_key not in frame_lines[frame_count].keys():
                            frame_lines[frame_count][net_sta_key] = \
                                {'lat': None, 'lon': None, 'x': None, 'y': None, 'u': None, 'v': None}
                        frame_lines[frame_count][net_sta_key]['lon'] = station_coordinates[net_sta_key][0]
                        frame_lines[frame_count][net_sta_key]['lat'] = station_coordinates[net_sta_key][1]
                        frame_lines[frame_count][net_sta_key]['x'] = x
                        frame_lines[frame_count][net_sta_key]['y'] = y
                        # Gain only applies to the vertical component only since we already have a quiver scale.
                        if comp_index == 1:
                            frame_lines[frame_count][net_sta_key]['v'] = tr2.data[0]  # * gain
                        else:
                            frame_lines[frame_count][net_sta_key]['u'] = tr2.data[0]  # * gain
                else:
                    # Even if we do not want to show markers, still need frame times
                    frames[frame_count][0] = \
                        {'x': list(), 'y': list(), 'z': list(), 't': list(), 's': list()}
                    frames[frame_count][0]['t'].append(slice_time.strftime('%Y-%m-%d %H:%M:%S UTC'))
                    frames[frame_count][0]['s'].append(frame_count * time_step)

            # NOTES: 1. Later we want to use the raw waveform so we can post the trace height. For now working
            # with the normalized trace is better. For now we are only using BMO station!!
            # 2. Grab the first station in case we cannot find the requested station later.
# Based on the new logic, we always keep the first trace. So we always have a reference trace or the first trace.
if True:
    cat_index = 0
    if verbose:
        utils.print_message(f'INFO',
                            f'Plotting the reference trace {reference_net}.{reference_sta}', log_file)

    ref_sta_x, ref_sta_y = bg_map(selected_station[sta_type]['lon'], selected_station[sta_type]['lat'])

    metadata.append(f"sta_lat={selected_station[sta_type]['lat']}")
    metadata.append(f"sta_lon={selected_station[sta_type]['lon']}")

    event_lat, event_lon = event_lat_lon[cat_index].strip().replace('(', ''). \
        replace(')', '').replace('', '').split(',')

    dist_meter, azim, back_azim = gps2dist_azimuth(float(event_lat), float(event_lon),
                                                   float(selected_station[sta_type]['lat']),
                                                   float(selected_station[sta_type]['lon']))
    dist_km = dist_meter / 1000.0
    dist_degree = kilometers2degrees(dist_km)
    metadata.append(f"sta_dist_deg={dist_degree}")
    metadata.append(f"sta_dist_km={dist_km}")
    metadata.append(f"sta_azim={azim}")

    if show_station_locations:
        bg_map.plot([ref_sta_x], [ref_sta_y], marker='o', markerfacecolor='none', linestyle='None',
                    markeredgewidth=0.3, markeredgecolor='black', alpha=0.3,
                    markersize=3, zorder=3, label='station')
    bg_map.plot([ref_sta_x], [ref_sta_y], marker='o',
                markerfacecolor='white',
                linestyle='None', markeredgewidth=0.3, markeredgecolor='None', alpha=1.0,
                markersize=0, zorder=0, label=f'gain {float(gain): 0.1f}')
    bg_map.plot([ref_sta_x], [ref_sta_y], marker='o',
                markerfacecolor=value_to_color(current_cmap, 1.0, 0.8),
                linestyle='None', markeredgewidth=0.3, markeredgecolor='None', alpha=1.0,
                markersize=3, zorder=3, label=top_color_label)

    bg_map.plot([ref_sta_x], [ref_sta_y], marker='o',
                markerfacecolor=value_to_color(current_cmap, 1.0, -0.8), linestyle='None',
                markeredgewidth=0.3, markeredgecolor='None', alpha=1.0, markersize=3, zorder=3,
                label=bottom_color_label)

    # Draw the 3C horizontal line for the legend.
    if request_comp == '3':
        quiver_line, = bg_map.drawgreatcircle(float(event_lon), float(event_lat), float(event_lon) + 1,
                                              float(event_lat) + 1,
                                              linewidth=0.8, color=horizontal_motion_color,
                                              label=horizontal_motion_text,
                                              zorder=1000, alpha=0)
    # Draw the Great Circle.
    if draw_great_circle:
        gc_line, = bg_map.drawgreatcircle(float(event_lon), float(event_lat),
                                          float(selected_station[sta_type]['lon']),
                                          float(selected_station[sta_type]['lat']),
                                          linewidth=1, color=great_circle_line_color,
                                          label='great circle path', zorder=1000)

    # Mark the reference station again just to get the station legend in.
    bg_map.plot([ref_sta_x], [ref_sta_y], marker='v', color=ref_sta_marker_color, linestyle='None',
                markersize=ref_sta_marker_size, markeredgecolor=ref_sta_marker_edge_color,
                markeredgewidth=ref_sta_marker_edge_width, zorder=500,
                label=f"reference station\n{selected_station[sta_type]['net_sta_key']}")

    metadata.append(f"sta_name={selected_station[sta_type]['net_sta_key']}")

    # Place the legend at loc of the bbox at the end of the trace and right justified..
    leg = plt.legend(numpoints=1, loc='upper left', bbox_to_anchor=(0.80, -0.05), fontsize=6, frameon=False,
                     framealpha=0)
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        if text.get_text().endswith(no_filter_text):
            text.set_color(no_filter_color)
        elif text.get_text().endswith(horizontal_motion_text):
            line.set_alpha(1)
            line.set_color(horizontal_motion_color)

    # Background of the elegend text is white to be readable when it is over other elements.
    # plt.setp(leg.get_texts(), backgroundcolor='w')

    # Swap and right align legend.
    vp = leg._legend_box._children[-1]._children[0]
    for child in vp._children:
        child._children.reverse()
    vp.align = 'right'

    # Location values are latitude,longitude.
    ax2.plot(selected_station[sta_type]['times'], selected_station[sta_type]['tr_no_filter'].data *
             selected_station[sta_type]['output_factor'], color=no_filter_color,
             label=f'{selected_station[sta_type]["tr"].stats.channel} {no_filter_text}', lw=0.5)
    ax2.plot(selected_station[sta_type]['times'], selected_station[sta_type]['tr'].data * float(gain) *
             selected_station[sta_type]['output_factor'], color='black',
             label=f'{selected_station[sta_type]["tr"].stats.channel} '
             f'filtered',
             lw=0.75)
    ax2.set_ylabel(f'{output_factor_symbol[int(selected_station[sta_type]["output_factor"])]}'
                   f'{output_units[output_type]}',
                   fontsize=7, labelpad=0)
    production_date = datetime.datetime.utcnow().replace(microsecond=0).isoformat()
    production_date = production_date.replace('T', ' ')
    metadata.append(f'movie_date={production_date} UTC')
    production_date = f'{production_label} {production_date} UTC - ' \
        f'{chan_label(request_bands, request_comp)}'

    ax2.text(0.0, -0.01, production_date, horizontalalignment='left', fontsize=5,
             verticalalignment='top', transform=ax2.transAxes)
    leg2 = ax2.legend(frameon=False, fontsize=6, loc='lower right', bbox_to_anchor=(1.0, 0.0))
    # Swap and right align legend.
    vp2 = leg2._legend_box._children[-1]._children[0]
    for child2 in vp2._children:
        child2._children.reverse()
    vp2.align = 'right'

    travel_times = dict()

    if technology == 'seismic':
        tt_url = tt_url.replace('EVENTLOC', f'{event_lat},{event_lon}')

        tt_url = tt_url.replace('STATIONLOC',
                                '{},{}'.format(selected_station[sta_type]['lat'],
                                               selected_station[sta_type]['lon']).strip().replace(' ', ''))
        tt_url = tt_url.replace('PHASE', phase_list.strip().replace(' ', ''))
        utils.print_message('INFO', f'Requesting travel times from:\n{tt_url}\n', log_file)
        travel_times, distance = get_travel_times(tt_url)
        prev_phase = 0
        # Zero of the trace time axis is in reference to the animation start time. Need to adjust
        # the time for that.
        for phase in list(travel_times.keys()):
            t_x = travel_times[phase] - animation_delay
            if 0 <= t_x <= animation_duration:
                # R1 and R2 are not considered in phase spacing
                if prev_phase == 0 or t_x - prev_phase >= phase_spacing or phase in ['R1', 'R2']:
                    ax2.axvline(x=t_x, color='blue', lw=0.5)
                    ax2.text(t_x, selected_station[sta_type]['max'] * 0.99, phase,
                             horizontalalignment='center',
                             verticalalignment='center',
                             size=6,
                             color='blue',
                             bbox=dict(facecolor='white', alpha=0.99, lw=0, pad=0))
                    if phase not in ['R1', 'R2']:
                        prev_phase = t_x
    ax2.set_xlim(left=0, right=animation_duration)
    print('\n', file=log_file, flush=True)

# Initialize markers for all colors.
station_markers = list()
if plot_markers:
    t0 = time_it(t0)
    if verbose:
        utils.print_message('INFO', '\n\nInitialize markers', log_file)
    for color_index, color_value in enumerate(full_color_list):
        m, = set_marker(bg_map, color_value)
        station_markers.append(m)

    if verbose:
        utils.print_message('INFO', 'Total of {} colors detected.'.format(len(station_markers)), log_file)
    t0 = time_it(t0)


def draw_contours(map_handle, lon_list, lat_list, value_list):
    """Animation function; this is called sequentially."""

    # The parameters that we need to update for global use
    global t0

    t0 = time_it(t0)
    # Must be np arrays for grid
    lon_list = np.array(lon_list)
    lat_list = np.array(lat_list)
    value_list = np.array(value_list)

    # Find the min and max of coordinates.
    lon_min = lon_list.min()
    lat_min = lat_list.min()
    lon_max = lon_list.max()
    lat_max = lat_list.max()

    # Now let's grid your data.
    # First we'll make a regular grid to interpolate onto.
    lon_new_num = int(((lon_max - lon_min) / spatial_resolution) + 1)
    lat_new_num = int(((lat_max - lat_min) / spatial_resolution) + 1)

    # Create a uniform mesh for contouring.
    lon_new = np.linspace(lon_min, lon_max, lon_new_num)
    lat_new = np.linspace(lat_min, lat_max, lat_new_num)

    # Basic mesh in lon, lat (degrees).
    lon_new, lat_new = np.meshgrid(lon_new, lat_new)

    try:
        # Interpolate at the points in lon_new, lat_new.
        # Method : {'linear', 'nearest', 'cubic'}, optional.
        value_new = griddata((lon_list, lat_list), value_list, (lon_new, lat_new), method=grid_method)

        # if nearest is selected for the method, we want to set the points outside of the convex hull of the
        # input points using linear method, since nearest does not do it.
        if grid_method == 'nearest':
            masked_values = griddata((lon_list, lat_list), value_list, (lon_new, lat_new), method='linear')
            masked_points = np.argwhere(np.isnan(masked_values))
            for mask_point in masked_points:
                value_new[mask_point[0]][mask_point[1]] = masked_values[mask_point[0]][mask_point[1]]
    except Exception as _er:
        utils.print_message('ERR', 'Gridding failed: {}'.format(_er), log_file)
        sys.exit(2)

    # Compute native map projection coordinates of lat/lon grid.
    grid_x, grid_y = map_handle(lon_new, lat_new)

    # Mask  ocean values..
    masked_values = maskoceans(lon_new, lat_new, value_new, resolution='h', grid=1.25, inlands=True)
    # Plot the contours, if there are any points.
    if len(masked_values) <= 3:
        return None

    ch = map_handle.contourf(grid_x, grid_y, masked_values,
                             cmap=current_cmap, alpha=contour_alpha,
                             vmin=-1.0, vmax=1.0, zorder=4)
    # ch = map_handle.contourf(grid_x, grid_y, value_new, colors=np.array(full_color_list),
    #                         levels= np.arange(-1,1,2.0/len(full_color_list)), alpha=.5, zorder=4)
    # ch2 = map_handle.drawlsmask(ocean_color='skyblue', land_color=(0, 0, 0, 0), lakes=True, zorder = 5)
    t0 = time_it(t0)

    return ch


def animate_markers(this_frame_index, these_markers):
    """Animation function; this is called sequentially."""
    if param.log_to_screen:
        print('[INFO] Animating frame {}/{}          '.format(this_frame_index, frame_count), end='\r', file=log_file,
              flush=True)
    else:
        log_file.write('{}/{}  '.format(this_frame_index, frame_count))

    # The parameters that we need to update for global use
    global x, y, interval_title_string, time_guide, time_guide_label, contour_handle, t0, show_station_locations, \
        event_marker, marker_lines

    interval_title_string = None

    # Update markers for each color in this frame.
    if plot_markers:
        for do_color_index in frames[this_frame_index].keys():

            # Set the frame title based on time from event.
            if interval_title_string is None:
                these_seconds = frames[this_frame_index][do_color_index]['s'][0] + animation_delay
                if len(event_time) == 1:
                    interval_title_string = '{} (OT {} {})'. \
                        format(frames[this_frame_index][do_color_index]['t'][0], sign(these_seconds, plus=True),
                               str(datetime.timedelta(seconds=int(round(abs(these_seconds))))))
                else:
                    interval_title_string = 'OT {} {}'. \
                        format(sign(these_seconds, plus=True),
                               str(datetime.timedelta(seconds=int(round(abs(these_seconds))))))
            # if event_marker is not None:
            #    event_marker.remove()
            #    event_marker = None
            this_frame_time = UTCDateTime(frames[this_frame_index][do_color_index]['t'][0].strip().replace(' UTC', ''))
            for _i, _values in enumerate(track_records):
                if not _values.strip():
                    continue

                _date, _time, _lat, _lon = _values.strip().replace(',', ' ').replace('/', '-').split()
                _month, _day, _year = _date.strip().split('-')
                if this_frame_time + time_step >= UTCDateTime('{}-{}-{} {}'.format(
                        _year.strip(), _month.strip(), _day.strip(), _time.strip())) >= this_frame_time:
                    _x, _y = bg_map(-1.0 * float(_lon), float(_lat))
                    for _m in marker_track:
                        event_marker = bg_map.plot([_x], [_y], marker=_m, color='y', linestyle='None',
                                                   markersize=8, zorder=500, markeredgecolor='k', markeredgewidth=0.5, )
                    break

            # Gather marker info for this frame.
            x = frames[this_frame_index][do_color_index]['x'].copy()
            y = frames[this_frame_index][do_color_index]['y'].copy()
            if len(x):
                these_markers[do_color_index].set_data(x, y)

    if request_comp == '3':
        q_x = list()
        q_y = list()
        q_u_in = list()
        q_v_in = list()
        q_lat = list()
        q_lon = list()
        for _key in frame_lines[this_frame_index].keys():

            if (frame_lines[this_frame_index][_key]['x'] is None) \
                    or (frame_lines[this_frame_index][_key]['y'] is None) \
                    or (frame_lines[this_frame_index][_key]['u'] is None) \
                    or (frame_lines[this_frame_index][_key]['v'] is None):
                continue
            q_x.append(frame_lines[this_frame_index][_key]['x'])
            q_y.append(frame_lines[this_frame_index][_key]['y'])
            q_u_in.append(frame_lines[this_frame_index][_key]['u'])
            q_v_in.append(frame_lines[this_frame_index][_key]['v'])
            q_lat.append(frame_lines[this_frame_index][_key]['lat'])
            q_lon.append(frame_lines[this_frame_index][_key]['lon'])

        if marker_lines is not None:
            marker_lines.remove()
        u_rotated, v_rotated = bg_map.rotate_vector(np.array(q_u_in), np.array(q_v_in),
                                                    np.array(q_lon), np.array(q_lat))

        marker_lines = bg_map.quiver(q_x, q_y, u_rotated, v_rotated, color=horizontal_motion_color, zorder=499,
                                     headwidth=0.001,
                                     scale_units=quiver_scale_units, width=proj_regions[region]['quiver_width'],
                                     scale=quiver_scale)

    # Plot contours for all contour regions.
    if plot_contours:
        for this_index in range(len(contour_lon)):
            lon_all = contour_data[this_frame_index][this_index]['lon'].copy()
            lat_all = contour_data[this_frame_index][this_index]['lat'].copy()
            val_all = contour_data[this_frame_index][this_index]['value'].copy()
            x_all = contour_data[this_frame_index][this_index]['x'].copy()
            y_all = contour_data[this_frame_index][this_index]['y'].copy()

            # Mask points outside this region.
            for _i, _lon_val in enumerate(lon_all):
                if min(contour_lat[this_index]) > lat_all[_i] or \
                        max(contour_lat[this_index]) < lat_all[_i] or \
                        min(contour_lon[this_index]) > lon_all[_i] or \
                        max(contour_lon[this_index]) < lon_all[_i]:
                    val_all[_i] = np.nan
                    continue

            if not plot_markers and show_station_locations:
                bg_map.plot(x_all, y_all, marker='o', markerfacecolor='none', linestyle='None',
                            markeredgewidth=0.3, markeredgecolor='black', alpha=0.3,
                            markersize=3, zorder=5, label='station')

            # Remove contours from the previous frame, if any.
            if contour_handle[this_index] is not None:
                for contour_item in contour_handle[this_index].collections:
                    contour_item.remove()
            if len(lon_all):
                contour_handle[this_index] = draw_contours(bg_map, lon_all, lat_all, val_all)
            else:
                contour_handle[this_index] = None
        show_station_locations = False
        if interval_title_string is None:
            use_color_index = 0
            these_seconds = frames[this_frame_index][use_color_index]['s'][0] + animation_delay
            interval_title_string = '{} (OT {} {})'. \
                format(frames[this_frame_index][use_color_index]['t'][0], sign(these_seconds, plus=True),
                       str(datetime.timedelta(seconds=int(round(abs(these_seconds))))))

    # Remove the time slider from the previous frame
    if time_guide is not None:
        time_guide.remove()
        time_guide_label.remove()

    # Add the time slider for this frame.
    use_color_index = 0
    if plot_markers:
        use_color_index = list(frames[this_frame_index].keys())[0]

    # If the requested ref station is found, plot the slider.
    time_guide = ax2.axvline(x=frames[this_frame_index][use_color_index]['s'][0],
                             linewidth=1, color="#FF0000")
    these_seconds = frames[this_frame_index][use_color_index]['s'][0] + animation_delay
    time_guide_label = ax2.text(frames[this_frame_index][use_color_index]['s'][0],
                                selected_station[sta_type]['max'] * 0.85,
                                'OT {} {}'.format(sign(these_seconds, plus=True),
                                                  str(datetime.timedelta(seconds=int(round(abs(these_seconds)))))),
                                horizontalalignment='center',
                                verticalalignment='center',
                                size=6,
                                color='#FF0000',
                                bbox=dict(facecolor='white', alpha=0.8, lw=0, pad=0))

    title.set_text(interval_title_string)
    ax2.set_xlim(left=0, right=animation_duration)
    t0 = time_it(t0)
    return these_markers


if verbose:
    utils.print_message('INFO', '\n\nData collection and processing time:', log_file)
t0 = time_it(t_run)

if verbose:
    utils.print_message('INFO', '\nStarting animation', log_file)

# Call the animator;  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate_markers, fargs=([station_markers]),
                               frames=len(frames), interval=time_step, blit=True, repeat=False)

metadata.append(f"movie_start={animation_start_time[0]}")
metadata.append(f"movie_end={animation_end_time[0]}")
metadata.append(f"movie_channel={request_chans}")

if output_file is None:
    animation_file_name = f'GMV_{request_chans.replace(",", "")}_{animation_start_time[0]}_{script_version}' \
        f'_{region}_{output_type}{do_event}' \
        f'_dt_{time_step}_C_{color_levels}'
    if plot_markers:
        animation_file_name = f'{animation_file_name}_markers'
    if plot_contours:
        animation_file_name = f'{animation_file_name}_contour_{spatial_resolution}_{grid_method}'
    if do_taper:
        animation_file_name = f'{animation_file_name}_T_{taper_fraction}'
    if scale_plot:
        animation_file_name = f'{animation_file_name}_X_{gain}'
    if std_check:
        animation_file_name = f'{animation_file_name}_SN_{std_max}'
    if animation_delay != 0.0:
        animation_file_name = f'{animation_file_name}_delay_{animation_delay}'
    if request_comp == '3':
        animation_file_name = f'{animation_file_name}_Q_{quiver_scale}'
    metadata.append(f"movie_label={animation_file_name}".replace(':', ''))
else:
    animation_file_name = output_file
    metadata.append(f"movie_label={animation_file_name}_{animation_start_time[0]}".replace(':', '').replace('-', '_'))

if metadata_dir is None:
    utils.print_message('WARN', '\nWill not write metadata file', log_file)
else:
    with open(os.path.join(metadata_dir, f'{animation_file_name}.txt').replace(':', ''), 'a') as fp:
        fp.write('\n'.join(metadata))

animation_file_name = os.path.join(video_dir, f'{animation_file_name}.mp4'.replace(':', ''))
utils.print_message('INFO', f'\n\nSaving Animation to {animation_file_name}', log_file)

writer = animation.FFMpegWriter(fps=frames_per_second, codec=codec, metadata=dict(artist=author), bitrate=bit_rate)

anim.save(animation_file_name, writer=writer)
t0 = time_it(t0)
utils.print_message('INFO', 'Total time:', log_file)
t_final = time_it(t_run, forced=True)
print(f"............{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} END............", file=log_file, flush=True)
log_file.close()
sys.exit(0)
