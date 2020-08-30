 Incorporated Research Institutions for Seismology (IRIS)
 Data Management Center (DMC)
 Data Products Team
 Ground Motion Visualization (GMV)

 2020-09-03

------------------------------------------------------------------------------------------------------------------------

 DESCRIPTION:

The Ground Motion Visualization (GMV, http://ds.iris.edu/ds/products/gmv/) is a video-based IRIS DMC data product that 
illustrates how seismic waves travel away from an earthquake location by animating the normalized recorded wave 
amplitudes at each seismometer location using colored markers. Color of each marker depicts amplitude of the vertical 
ground motion, as detected by the station’s seismometer and normalized to its peak amplitude. 

This Python bundle is the main code used to produce such GMVs for different geographical regions, sensors, technologies 
The GMV production script (_gmv\_generalized.py_) can be configured via its parameter file (_gmv\_param.par_) or 
command line arguments. Currently parameters are optimized for use with the Lambert conformal map projection and 
seismic channels but with additional parameter tuning it is possible to change projection and/or technology. The code 
uses the FDSN Web Services (https://www.fdsn.org/webservices/) to retrieve waveform data from different FDSN data 
centers (https://service.iris.edu/irisws/fedcatalog/1/datacenters).


This bundle contains the following files:

     src/
       _gmv\_generalized.py_
           - a Python version of the  original MATLAB GMV code that was originally used in GMV
             production. Calling the code without any parameters displays a list of options available to tune GMV 
             production. It also provides two examples to run.
          
       _gmv\_param.py_
           - a Python file contains all the parameters available to configure GMV production. All 
             parameter definitions must follow Python rules. Each parameter group is commented for clarification.
     
       - _gmv\_utils.py_—
           - a Python utility library used by the main script.


    CHANGES.txt
       - a text file containing the history of changes to this bundle

    INSTALL.txt
       - installation notes

    README.md
       - this file


 INSTALLATION:

    see the INSTALL.txt file


USAGE:
   
    gmv_generalized.py --band=LH,BH --comp=1 -n all -t 2020-07-22T06:12:42 -T 'July 22, 2020, Alaska Peninsula, M 7.8' -m 7.8 -z 10.0 -e 55.1046,-158.4725 -r ak -d 1200 -s 6.0 -p -180 -q 3.5 -g 3 -D 0.05 -G -o GMV_Example_Custom


    -h  --help		Output this message.
    -v, --verbose           [default: False] Turn on the verbose mode.
    -b, --band		[default: LH,BH] Channel bands to use for GMV production (separate by comma, for example: LH OR LH,BH OR  LH,BH,HH).
    -c, --comp		[default: 1] Number of channel components (1 or 3).
    -d, --dur		[default: 2400] Duration of the GMV animation in seconds.
    -D, --std		[default: 0.05] Maximum acceptable standard deviation for trace vetting.
    -e, --eloc		[*required] Single quoted event location as 'lat,lon' (example:  '10.779,-62.907'). For multiple events., 
        			separate each set by a space (example: '24.760,-109.890 25.200,-109.870'. Use of parentheses is optional: '(24.760,-109.890) (25.200,-109.870)'.
    -g, --gain		[default: 3.0]Trace amplification to generate GMVs. For 3C GMVs, gain is only applied to the Z-component.
    -G, --gc		[default: False] raw the great circle path between the event location and the reference station.
    -l, --tstep             [default: 3.0] Time step in seconds to use for sampling traces and create the video. Example: -l 2 samples traces every 2 seconds.
    -m, --emag		[*required] Event magnitude. For multiple events, single quoted magnitude list and separate them with  a space (example: -m '7.3 6.2'.
    -n, --net		[default: all] Network(s) to request data from. For multiple networks, separate them with comma. Use all to request from all networks (example: -n TA OR -n US,TA -n all).
    -N, --rnet		[default: IW] Network of the reference station.Example -N TA.
    -o, --output	        The output video file name (by default the video will be MP4 format. Example: -o Test_video will output Test_video.mp4.
    -p, --delay	        [default: -20] Delay in seconds from the event origin time to start the video. Example: -p 120.
    -P, --phase	        [default: 200] Seconds between Phases marked on the trace. This is used to avoid overprinting phase labels. Example -P 30.
    -q, --qscale	        [default: 4.0] Quiver scale for 3C plots. Number of data units per arrow length unit, e.g., m/s per plot width; a smaller scale parameter makes the arrow longer.
    -r, --region	        [default: na] The region to create the GMV for. The selected region must be a key of the proj_regions. The acceptable regions are: dict_keys(['af', 'ak', 'am', 'as', 'au', 'aua', 'eur', 'gl', 'np', 'na', 'sa', 'us']).
    -s, --sizem	        [default: 4.0] Marker (https://matplotlib.org/api/markers_api.html) size in points.
    -S, --rsta		[default: MFID] Reference station code to plot.
    -t, --etime	        [*required] The event time as YYYY-MM-DDTHH:MM:SS.
    -T, --title	        [default: based on the event magnitude (-m) and location (-e)] Single quoted GMV title.
    -z, --depth	        [*required] Event depth in km.
 
For additional configuration, please see the parameter file (gmv_param.py)


EXAMPLES:

	1. Sample request with the least number of arguments:
		gmv_generalized.py -e 55.1046,-158.4725 -z 10.0 -m 7.8 -t 2020-07-22T06:12:42 -o GMV_Example_Default
	2. Sample complete request:
		gmv_generalized.py --band=LH,BH --comp=1 -n all -t 2020-07-22T06:12:42 -T 'July 22, 2020, Alaska Peninsula, M 7.8' -m 7.8 -z 10.0 -e 55.1046,-158.4725 -r ak -d 1200 -s 6.0 -p -180 -q 3.5 -g 3 -D 0.05 -G -o GMV_Example_Custom

SELECTING PARAMETERS:

 GMVs are event-based animation and therefore you need to start with the event parameter. As a minimum, you need to 
 provide the fallowing event parametrs:
 
 - location in the form of _latitude,longitude_ using the _-e_ option (-e 55.1046,-158.4725)
 - depth in kilometers using the _-z_ option (-z 10.0)
 - magnitude using the _-m_ option (-m 7.8). The magnitude lets the script select an optimum filter band to create GMVs
 - time using the _-t_ option (-t 2020-07-22T06:12:42). Time must follow the  yyyy-mm-ssThh24:mm:ss format and must be in UTC time
 
All other parameters are optional. However, by looking at the default values in the parameter file, you may decide to override 
them. For example change:
 - frequency band using _-b_ option
 - number of channel components (_-c_) used for GMV creations (single, Z or three component, Z and the two horizontals))
 - duration of the video in seconds (_-d_). You do not want to make duration too long, otherwise you will get a large animation 
   file without much ground motion at the later part. This will also affect the processing time for animation creation
   The value will depend on the location of the event with respect to the stations. Always start with a shorter duration and
   increase it if needed
 - trace QC via standard deviation (-D). To enhance animation quality and bring out lower amplitude arrivals, GMV uses a 
   trace normalization  for each station. This allows traces far from event location display arrival time as strong as 
   traces closest to the source. However, for stations with low S/N this will result in noise amplification that results in
   station markers that continuously blink during animation. To avoid this, the code uses a standard deviation algorithm
   that will reject traces that have high standard deviation before event tim. The default value of 0.05 appears to be very
   effective to reduce blinking. However if you notice that you are missing too many stations (increase this value using the _-d_
   option). If you think there are still too many stations that blink all the time, reduce the default value using the _-D_ option.
 - To bring out the weak phases, you can increase the gain applied to traces (gain is only applied to the vertical component) using the _-g_ option.
 - To display the great circle path on the animation panel, use the _-G_ option. No value is needed. Just using the _-G_ option turns the great circle path on
 - The _-l_ option allows you to control the time sampling option. By reducing the time step, you will create a smoother animation but at the
   same time you will increase the total size of the animation by including more frames. Sampling of 2-4 second appears to work reasonably well.
 - Sometimes it is desirealbel to limit the number of networks contributing to the animation. The _-n_ allows you to limit data requests only
   to the list you provide (separate multiple networks with comma). _-n all_ will request data from all available networks. When tuning the GMV,
   use a network with fewer number of stations to allow speedy processing. Once you are happy with the parameters, set the network list to _all_ or desired list.
 - You can select a reference station for display its trace below GMV panel. _-N_ is used to indicate the network of the reference station (see also the _-S_ option below).
 - You can provide the animation file name using the _-o_ option.
 - To avoid long initial video segment without major phase arrivals, use '_-p_' option to indicate how many seconds  after event origin time the video should start.
 - To avoid overprinting phase labels. use the _-P_ option to indicate the number of seconds between Phases marks on the trace. 
 - Quivers are used to depict horizontal motions when creating 3-component animations. The _-q_ option provides the number of data units per arrow 
   length unit, e.g., m/s per plot width; a smaller scale parameter makes the arrow longer.
 - To indicate the region for which the GMV should be created, use the _-r_ option. The selected region must be a key of the _proj\_regions_ dictionary in the param file. 
   The acceptable regions are: dict_keys(['af', 'ak', 'am', 'as', 'au', 'aua', 'eur', 'gl', 'np', 'na', 'sa', 'us']. You can define your own region by adding it to
   the _proj\_regions_ dictionary.
 - Marker size for stations should be selected proportional to the region. Selecting very small markers will make it difficult to correlate wave passage between stations and 
   very large marker size dominate animation and will distract the viewer. The _proj\_regions_ dictionary provides a reasonable marker size for each region. Howerev, you
   can override this value by using the _-s_ option that provides the marker size in points.
 - Reference station code to plot is provided using the _-S_ option. The station code along with the network code _-N_ option are used for selection of the reference station trace.
 - To override the default animation title that is based on the event magnitude (_-m_) and location (_-e_), use the _-T_ option. This option takes a single quoted GMV title as its argument.

To avoid  making one single large request for data to a data center, it is better to make multiple requests. The parameter _chunck_count_ in the parameter file determines the number of 
stations per request (chunk) that will be sent to each data center. This number should be adjusted based on the number of station-channels involved 
and the length of each request to avoid memory issues or data center timeouts.


CITATION:

To cite the use of this software please cite:
Trabant, C., A. R. Hutko, M. Bahavar, R. Karstens, T. Ahern, and R. Aster (2012), Data Products at the IRIS DMC: Stepping Stones for Research and Other Applications, Seismological Research Letters, 83(5), 846–854, https://doi.org/10.1785/0220120032.

Or cite the following DOI:
    10.17611/DP/USAGMV.1

CREDITS:

    - Chuck Ammon, Professor of Geosciences at Penn State’s original concept and visualizations.
    - Bob Woodward at IRIS – adapted the visualization code to MATLAB.
    - IRIS DMC Data products expanded and enhanced the MATLAB code.
    - IRIS DMC Data products Python code conversion.

 HISTORY
 - 2020-09-03 - Initial public release R.1.1
 - 2020-06-01 - Installed Generalized GMV 1.0
 - 2015-10-05 - Expanded to include Alaska
 - 2013-10-22 - Super GMVs online
 - 2011-03-24 - Customized GMV online
 - 2010-10-20 - 3-component GMVs online
 - 2010-02-25 - GMV online, automated production of GMV started

 
 COMMENTS/QUESTIONS:

    Please contact manoch@iris.washington.edu


