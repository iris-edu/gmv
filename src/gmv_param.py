import os
from datetime import datetime

"""     
Name:
    gmv_param.py
    
Description:
    GMV Parameter file
    
History:
    2020-09-03 Manoch: V.2020.241 R1.1 Public release.
    2020-07-27 Manoch: V.2020.209 introduced chunking for data request to enable smaller request chunks.
    2020-07-06 Manoch: V.2020.189 introduced North Polar region and STD water level
    2020-05-30 Manoch: V.2020.152 R1 release
    2020-05-19 Manoch: changed the ref sta from eur and gl coverages
    2020-05-11 Manoch: Initial release.
"""

version = 'V.2020.241'

# The run settings.
log_to_screen = True
verbose = False


# Set the main GMV script, necessary directories and file.
gmv_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
log_dir = os.path.join(gmv_path, 'log')
metadata_dir = os.path.join(gmv_path, 'metadata')
image_dir = os.path.join(gmv_path, 'image')
video_dir = os.path.join(gmv_path, 'video')

# logo parameters.
logo_file = 'iris_color_screen.png'
logo_zoom = 0.2
logo_location = (-1.0, 2.0)


# Month names.
month_names = ('', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
               'November', 'December')

# For some sources, like hurricanes, it is desired to post the source (storm) location as a function of time.
# Each line of the source_track_file  must contain:
# date, time, latitude, longitude
source_track_file = None

# Get the correct time.
this_date = datetime.now()
time_mark = this_date.strftime('%Y-%m-%dT%H:%M:%S%z')

# Starting magnitude for the large events for which 3C GMVs are created.
large_magnitude = 7.0

# Plot parameters.
# Basemap type to plot on. Options are: [etopo, releaf, solid]
basemap_type = 'solid'

line_width = {'coastlines': 0.01, 'countries': 0.1, 'states': 0.05, 'grids': 0.01}
color = {'lake': '#ADCBE5', 'land': '#FDF8E6', 'seismic': '#7B0007'}

dpi = 250
frames_per_second = 25


# What kind of data we are working with?
echnology = 'seismic'
channels = {'seismic': '?H?,?L?,?N?,?M?,?P?,?DH', 'all': '*', 'infrasound': 'BDF,LDF'}
marker = {'seismic': 'o', 'all': 'o', 'infrasound': 'o'}

# Video settings.
marker_size = 1.25
frame_seconds = 30.0 * 86400.0  # Each frame represents ~ month
frame_interval = 20  # Delay between frames in milliseconds.
author = 'IRIS'

start = '1970-01-02'
end = '2019-05-01'
include_restricted = 'false'
title_text = "Stations  from IRISWS fedcatalog\n"
title_style = 'italic'
title_font_size = 8
title_color = 'black'
title_lat_lon = (7, -139)
bbox_alpha = 0.0
bbox_padding = 2

projection = 'lcc'  # 'robin'

# Some needed URLs.
url_data_centers = 'http://service.iris.edu/irisws/fedcatalog/1/datacenters'

url_fedcatalog = 'http://service.iris.edu/irisws/fedcatalog/1/query?net=NET&sta=STATION&loc=*' \
                 '&cha=CHAN&targetservice=station' \
                 '&level=station&format=text&endafter=END&startbefore=START&includeoverlaps=false&nodata=404'

# The color pallet to use.
colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#c88543', '#563296', '#f032e6', '#bcf60c',
          '#fabebe', '#008080', '#e6beff', '#9a6324', '#46f0f0', '#800000', '#aaffc3', '#808000', '#ffd8b1',
          '#000075', '#000000', '#808080', '#fffac8', '#911eb4']

production_label = 'IRIS DMC GMV'

# Trace plotting parameters.
no_filter_color = '#BEBDB8'
no_filter_text = 'no filter'
horizontal_motion_text = 'horizontal motion'
horizontal_motion_color = 'black'
output_type = 'DISP'
draw_great_circle = False
great_circle_line_color = 'red'
zero_base = 3e-3

# Should we apply gain to the traces?
scale_plot = True
gain = 3.0

# The default network to request data from.
request_net = '_US-ALL'


# Reference station location plot parameters.
ref_sta_marker_color = 'green'
ref_sta_marker_size = 8
ref_sta_marker_edge_color = 'yellow'
ref_sta_marker_edge_width = 1

# Default region for GMV.
region = 'na'

# Individual frame sizes.

figure_width = 6.5
subplot_columns = 6
trace_height = 1.5

subplot_map_rows = 1
subplot_trace_rows = 1
subplot_rows = subplot_map_rows + subplot_trace_rows

if trace_height is None:
    subplot_map_rows = 3
    subplot_trace_rows = 1
    subplot_rows = subplot_map_rows + subplot_trace_rows
else:
    subplot_rows = 20

""" Output units. One of:
    'DISP': displacement, output unit is meters
    'VEL': velocity, output unit is meters/second
    'ACC': acceleration, output unit is meters/second**2
"""
output_type = 'VEL'
output_factor_symbol = {1000: 'm', 1000000: u'\u03bc', 1000000000: u'\u03BC'}
output_units = {'VEL': 'cm/s', 'DISP': 'm', 'ACC': 'm/s^2', 'PA': 'Pascal', 'C': 'Degrees Celsius'}
channel_comp_index = {'Z': 0, '1': 1, 'N': 1, '2': 2, 'E': 2}
channel_band = 'LH,BH'
channel_band_list = channel_band.split(',')
channel_comp_list = {'1': ['Z'], '3': ['1', '2', 'Z', 'N', 'E']}

request_comp = '1'

# To avoid  making one single large request for data to a data center, it is better to make multiple requests.
# The parameter _chunck_count_ in the parameter file determines the number of
# stations per request (chunk) that will be sent to each data center. This number should be adjusted based on the
# number of station-channels involved and the length of each request to avoid memory issues or data center timeouts.
chunk_count = 10

network = '*'
station = '*'

marker_event = '*'
marker_track = ['1', '2', '3', '4']
marker_size_event = 8
marker_size_track = 6

author = 'IRIS'
bit_rate = 1800
time_step = 3.0

# Sign of the delay: -X start X seconds sooner; +x start X seconds later
animation_delay = -20
animation_duration = 2400
phase_spacing = 200

top_color_label = 'motion up'
bottom_color_label = 'motion down'

timer_is_on = False

request_padding = 300

frames_per_second = 25

codec = 'h264'
technology = 'seismic'

contour_alpha = 1.0

travel_time_model = 'iasp91'

# If not None, trace is normalized by dividing by specified value norm instead of dividing by its absolute maximum.
# If a negative value is specified then its absolute value is used. If it is zero (either through a zero array or
# by being passed), nothing will happen and the original array will not change. std_water_level is the normalized
# amplitude which would trigger STD check. If max amplitude within the STD window is below this, then it passes.

std_check = True
std_water_level = 0.1
std_window = 300
std_max = 0.05

# Width of a cosine taper to apply to the waveform data in time domain prior to deconvolution.
# Specified as a fraction of the trace length from 0 to 0.5.
do_taper = True
taper_fraction = 0.5

# dc_to_exclude = ['NCEDC', 'SCEDC']
dc_to_exclude = []

# plot_type = sys.argv[2]
plot_type = 'Marker'
net_to_exclude = ['SY']

show_station_locations = True

filter_freqmin = {'large': 1.0 / 500.0, 'medium': 1.0 / 250.0, 'small': 1.0 / 50.0, 'infrasound': 0.06, 'rocket': 2.0,
                  'explosion': 2.0, 'barom': 1.0 / 6.0 / 3600.0, 'super': 1.0 / 500.0}
filter_freqmax = {'large': 1.0 / 100.0, 'medium': 1.0 / 50.0, 'small': 1.0 / 20., 'infrasound': 0.25, 'rocket': 5.0,
                  'explosion': 5.0, 'barom': 1.0 / 2.0 / 3600.0, 'super': 1.0 / 50.0}

min_trace_length = 300  # 5 minutes

spatial_resolution = 1.0
grid_method = 'nearest'

# List of seismic phases to calculate travel times for
phase_list = 'P,p,Pn,Pg,Pdiff,pPdiff,sPdiff,PcP,PP,PP,PS,PKP,SKS,S,s,Sn,Sg,Sdiff,pSdiff,sSdiff,SP,ScS,SS,SKS'

# Number of seconds between two consecutive phases to avoid over plots
# Rayleigh slowness (s/degree); originally 29.835 but using slowness of 28.5 to provides a better phase display
r_slowness = 28.5

# Iris Travel Time service URL.
iris_traveltime_url = f'http://service.iris.edu/irisws/traveltime/1/query?phases=PHASE&noheader=true&evdepth=DEPTH' \
    f'&evloc=[EVENTLOC]&staloc=[STATIONLOC]&model={travel_time_model.strip()}'

location_order = {'': 1, '  ': 2, '--': 3, '00': 4}
if technology == 'infrasound':
    location_order = {'EP': 1}
elif technology == 'temperature':
    location_order = {'EP': 1}
exclude_temporary_networks = True

# The number of rgb quantization levels
color_levels = 32
sampling_rate = 1.0

# Define a filter band to prevent amplifying noise during the deconvolution.
if technology == 'infrasound':
    pre_filter = (0.0033, 0.004, 0.5, 0.6)
else:
    pre_filter = (0.0008, 0.001, 0.08, 0.1)
verbose = False


"""  Lambert conformal map parameter
     name: name of the region
     corners: region's LL and UR corners as (llcrnrlat, llcrnrlon), (urcrnrlat, urcrnrlon))
     limits: provides limits for the projection box in which area remains the same
     limits sre given as ((lat1, lat2), (lon1, lon2))
     center: map's center point (lat_0, lon0)
     ref_sta: the default reference station for the region
     projection: Lambert conformal (lcc)
     size: width, height of the figure (provides the aspect ration)
     quiver_scale_units: quiver scale units for 3C quiver markers
     quiver_scale: scaling factor for quivers. Larger numbers will make quivers smaller
     quiver_width: width of the quiver to avoid plotting the arrow head
     marker_size: marker size for station markers
"""
proj_regions = dict()

proj_regions['af'] = dict()
proj_regions['af']['name'] = 'Africa'
proj_regions['af']['corners'] = ((-39.0, -20.0), (40.0, 55.0))
proj_regions['af']['limits'] = ((-35, 40), (-20, 55))
proj_regions['af']['center'] = (20.0, 30.0)
proj_regions['af']['ref_sta'] = ('AF', 'POGA')
proj_regions['af']['projection'] = 'lcc'
proj_regions['af']['size'] = (3.4, 3.7)
proj_regions['af']['quiver_scale_units'] = 'inches'
proj_regions['af']['quiver_scale'] = 2.0
proj_regions['af']['quiver_width'] = 0.003
proj_regions['af']['marker_size'] = 6.0

proj_regions['ak'] = dict()
proj_regions['ak']['name'] = 'Alaska'
proj_regions['ak']['corners'] = ((50.0, -170.0), (75.0, -130.0))
proj_regions['ak']['limits'] = ((50.0, 75.0), (-170.0, -130.0))
proj_regions['ak']['center'] = (60.0, -140.0)
proj_regions['ak']['ref_sta'] = ('AK', 'MARJ')
proj_regions['ak']['projection'] = 'lcc'
proj_regions['ak']['size'] = (3.85, 3.7)
proj_regions['ak']['quiver_scale_units'] = 'inches'
proj_regions['ak']['quiver_scale'] = 2.0
proj_regions['ak']['quiver_width'] = 0.003
proj_regions['ak']['marker_size'] = 6.0

proj_regions['am'] = dict()
proj_regions['am']['name'] = 'Americas'
proj_regions['am']['corners'] = ((-60.0, -150.0), (70.0, 20.0))
proj_regions['am']['limits'] = ((-60.0, 80.0), (-170.0, -30.0))
proj_regions['am']['center'] = (0.0, -100.0)
proj_regions['am']['ref_sta'] = ('IW', 'MFID')
proj_regions['am']['projection'] = 'lcc'
proj_regions['am']['size'] = (2.9, 3.7)
proj_regions['am']['quiver_scale_units'] = 'inches'
proj_regions['am']['quiver_scale'] = 2.0
proj_regions['am']['quiver_width'] = 0.003
proj_regions['am']['marker_size'] = 4.0

proj_regions['as'] = dict()
proj_regions['as']['name'] = 'Asia'
proj_regions['as']['corners'] = ((-15.0, 45.0), (42.0, 155.0))
proj_regions['as']['limits'] = ((0.0, 80.0), (25.0, 179.9))
proj_regions['as']['center'] = (40.0, 90.0)
proj_regions['as']['ref_sta'] = ('CB', 'GTA')
proj_regions['as']['projection'] = 'lcc'
proj_regions['as']['size'] = (4.7, 3.7)
proj_regions['as']['quiver_scale_units'] = 'inches'
proj_regions['as']['quiver_scale'] = 2.0
proj_regions['as']['quiver_width'] = 0.003
proj_regions['as']['marker_size'] = 4.0

proj_regions['au'] = dict()
proj_regions['au']['name'] = 'Australia'
proj_regions['au']['corners'] = ((-40.0, 110.0), (-8.0, 152.5))
proj_regions['au']['limits'] = ((-40.0, -10.0), (110.0, 160.0))
proj_regions['au']['center'] = (-25.0, 132.5)
proj_regions['au']['ref_sta'] = ('AU', 'WRKA')
proj_regions['au']['projection'] = 'lcc'
proj_regions['au']['size'] = (4.4, 3.7)
proj_regions['au']['quiver_scale_units'] = 'inches'
proj_regions['au']['quiver_scale'] = 2.0
proj_regions['au']['quiver_width'] = 0.003
proj_regions['au']['marker_size'] = 6.0

proj_regions['aua'] = dict()
proj_regions['aua']['name'] = 'Australasia'
proj_regions['aua']['corners'] = ((-50.0, 90.0), (10.0, 175.0))
proj_regions['aua']['limits'] = ((-60.0, 15.0), (86.0, 179.9))
proj_regions['aua']['center'] = (-25.0, 135.0)
proj_regions['aua']['ref_sta'] = ('AU', 'WRKA')
proj_regions['aua']['projection'] = 'lcc'
proj_regions['aua']['size'] = (4.4, 3.7)
proj_regions['aua']['quiver_scale_units'] = 'inches'
proj_regions['aua']['quiver_scale'] = 2.0
proj_regions['aua']['quiver_width'] = 0.003
proj_regions['aua']['marker_size'] = 6.0

proj_regions['eur'] = dict()
proj_regions['eur']['name'] = 'Europe'
proj_regions['eur']['corners'] = ((32.0, -10.0), (70.0, 60.0))
proj_regions['eur']['limits'] = ((30.0, 85.0), (-25.0, 70.0))
proj_regions['eur']['center'] = (50.0, 10.0)
proj_regions['eur']['ref_sta'] = ('GE', 'STU')
proj_regions['eur']['projection'] = 'lcc'
proj_regions['eur']['size'] = (3.0, 3.7)
proj_regions['eur']['quiver_scale_units'] = 'inches'
proj_regions['eur']['quiver_scale'] = 2.0
proj_regions['eur']['quiver_width'] = 0.003
proj_regions['eur']['marker_size'] = 6.0

proj_regions['gl'] = dict()
proj_regions['gl']['name'] = 'Global'
proj_regions['gl']['corners'] = ((-90, -180.0), (90.0, 180))
proj_regions['gl']['limits'] = ((-90, 90.0), (-180.0, 180))
proj_regions['gl']['center'] = (0, 0)
proj_regions['gl']['ref_sta'] = ('GE', 'STU')
proj_regions['gl']['projection'] = 'robin'
proj_regions['gl']['size'] = (6.709, 3.7)
proj_regions['gl']['quiver_scale_units'] = 'inches'
proj_regions['gl']['quiver_scale'] = 2.0
proj_regions['gl']['quiver_width'] = 0.003
proj_regions['gl']['marker_size'] = 2.0

# setup north polar aimuthal equidistant basemap.
# The longitude lon_0 is at 6-o'clock, and the
# latitude circle boundinglat is tangent to the edge
# of the map at lon_0.
proj_regions['np'] = dict()
proj_regions['np']['name'] = 'North Polar'
proj_regions['np']['corners'] = ((-90, -180.0), (90.0, 180))
proj_regions['np']['limits'] = ((-90, 90.0), (-180.0, 180))
proj_regions['np']['center'] = (0, 0)
proj_regions['np']['ref_sta'] = ('GE', 'STU')
proj_regions['np']['projection'] = 'npaeqd'
proj_regions['np']['size'] = (6.709, 6.709)
proj_regions['np']['quiver_scale_units'] = 'inches'
proj_regions['np']['quiver_scale'] = 2.0
proj_regions['np']['quiver_width'] = 0.003
proj_regions['np']['marker_size'] = 2.0

proj_regions['na'] = dict()
proj_regions['na']['name'] = 'North America'
proj_regions['na']['corners'] = ((0.0, -150.0), (50.0, 0.0))
proj_regions['na']['limits'] = ((10.0, 80.0), (-170.0, -50.0))
proj_regions['na']['center'] = (50.0, -107.0)
proj_regions['na']['ref_sta'] = ('IW', 'MFID')
proj_regions['na']['projection'] = 'lcc'
proj_regions['na']['size'] = (4.8, 3.7)
proj_regions['na']['quiver_scale_units'] = 'inches'
proj_regions['na']['quiver_scale'] = 4.0
proj_regions['na']['quiver_width'] = 0.001
proj_regions['na']['marker_size'] = 4.0

proj_regions['sa'] = dict()
proj_regions['sa']['name'] = 'South America'
proj_regions['sa']['corners'] = ((-55.0, -110.0), (20.0, -35.0))
proj_regions['sa']['limits'] = ((-60, 20), (-85, -30))
proj_regions['sa']['center'] = (-20.0, -60.0)
proj_regions['sa']['ref_sta'] = ('BR', 'PTBL')
proj_regions['sa']['projection'] = 'lcc'
proj_regions['sa']['size'] = (2.9, 3.7)
proj_regions['sa']['quiver_scale_units'] = 'inches'
proj_regions['sa']['quiver_scale'] = 2.0
proj_regions['sa']['quiver_width'] = 0.001
proj_regions['sa']['marker_size'] = 4.0

proj_regions['us'] = dict()
proj_regions['us']['name'] = 'US'
proj_regions['us']['corners'] = ((20.0, -125.0), (50.0, -60.0))
proj_regions['us']['limits'] = ((25.0, 57.0), (-125.0, -66.5))
proj_regions['us']['center'] = (35.0, -107.0)
proj_regions['us']['ref_sta'] = ('IW', 'MFID')
proj_regions['us']['projection'] = 'lcc'
proj_regions['us']['size'] = (4.8, 3.7)
proj_regions['us']['quiver_scale_units'] = 'inches'
proj_regions['us']['quiver_scale'] = 2.0
proj_regions['us']['quiver_width'] = 0.001
proj_regions['us']['marker_size'] = 6.0
