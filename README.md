 Incorporated Research Institutions for Seismology (IRIS)\
 Data Management Center (DMC)\
 Data Products Team\
 Ground Motion Visualization (GMV)

 2020-09-03

------------------------------------------------------------------------------------------------------------------------

 DESCRIPTION:

The Ground Motion Visualization (GMV, http://ds.iris.edu/ds/products/gmv/) is a video-based IRIS DMC data product that 
illustrates how seismic waves travel away from an earthquake location by animating the normalized recorded wave 
amplitudes at each seismometer location using colored markers. 
Color of each marker depicts amplitude of the vertical ground motion, as detected by the station’s seismometer and 
normalized to its peak amplitude. For seismic channels, either single-component (vertical or Z component) or three-component
(vertical plus two horizontal components) animations are possible (see http://ds.iris.edu/spud/gmv/18288940).

This Python bundle is the main code behind GMV production at IRIS DMC (http://ds.iris.edu/ds/newsletter/vol22/no2/522/generalized-gmvs-post-ta-ground-motion-visualizations/). The GMV production script (_gmv\_generalized.py_) can be configured via its parameter 
file (_gmv\_param.par_) or through command line arguments. Currently parameters are optimized for use with the Lambert 
conformal map projection and  seismic channels. However, with additional parameter tuning, it is possible to 
change the projection and/or the sensor technology. The code uses the FDSN Web Services 
(https://www.fdsn.org/webservices/) to retrieve waveform data from different FDSN data centers 
(https://service.iris.edu/irisws/fedcatalog/1/datacenters).


This bundle contains the following files:

     src/
       gmv_generalized.py
           - This is the main GMV production code. It is a Python version of the original MATLAB code used in GMV
             production. Calling the code with  -h  optoin displays a list of other options available to tune GMV 
             production. It also provides test examples to run.
          
       gmv_param.py
           - A Python file that contains all GMV parameters. You could modify to customize GMV production. All 
             parameter definitions in this file must follow Python rules. Each parameter group in this file is 
             commented for clarification.
     
       - gmv_utils.py—
           - A Python utility library used by the main script.


    CHANGES.txt
       - A text file containing the history of changes to this bundle.

    INSTALL.txt
       - The installation notes

    README.md
       - This file


 INSTALLATION:

    see the INSTALL.txt file


USAGE:
   
    gmv_generalized.py --band=LH,BH --comp=1 -n all -t 2020-07-22T06:12:42 -T "July 22, 2020, Alaska Peninsula, M 7.8" 
                       -m 7.8 -z 10.0 -e 55.1046,-158.4725 -r ak -d 1200 -s 6.0 -p -180 -q 3.5 -g 3 -D 0.05 
                       -G -o GMV_Example_Custom


    -h  --help		Output this message.
    -v, --verbose           [default: False] Turn on the verbose mode (no value needed).
    -b, --band		[default: LH,BH] Two-character channel bands to use in GMV production (separate by comma, 
                            for example to select from LHZ, BHZ, and HHZ channels, use : LH OR LH,BH OR  LH,BH,HH). 
                            See band and instrument code in https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/. 
    -c, --comp		[default: 1] Number of channel components (1: Z direction only, | 3: Z,N,E or Z,1,2 ).
    -d, --dur		[default: 2400] Duration of the GMV animation in seconds.
    -D, --std		[default: 0.05] Maximum acceptable standard deviation (SD) for trace vetting. Traces with SD for a 
                            window prior to the event's origin time larger than this value will not be included in 
                            animation ( will reduce blinking effect in GMV).
    -e, --eloc		[*required] Event location as lat,lon (example:  10.779,-62.907). For multiple events (super GMV, 
                            https://ds.iris.edu/ds/products/usarraygmv-super/), the value should be double-quoted and 
                            sets separated by a space (example: "24.760,-109.890 25.200,-109.870". Use of parentheses 
                            is optional: "(24.760,-109.890) (25.200,-109.870)".
    -g, --gain		[default: 3.0] Trace amplification to generate GMVs. For 3C GMVs, gain is only applied to the 
                            Z-component.
    -G, --gc		[default: False] Draw the great circle path between the event location and the reference station 
                            (no value needed).
    -l, --tstep             [default: 3.0] Time step in seconds to use for sampling traces and creating the video. 
                            Example: -l 2 samples traces every 2 seconds.
    -m, --emag		[*required] Event magnitude. For multiple events (Super GMV), double-quoted magnitude list and 
                            separate them with  a space (example: -m "7.3 6.2".
    -n, --net		[default: all] Network(s) to request data from. For multiple networks, separate them with comma. 
                            Use "all" to request from all networks (example: -n TA OR -n US,TA -n all).
    -N, --rnet		[default: IW] Network of the reference station. Example -N TA.
    -o, --output	        The output video file name (by default the video will be in MP4 format). 
                            Example: "-o Test_video" will output Test_video.mp4.
    -p, --delay	        [default: -20] Delay in seconds from the event origin time to start the video. Example: -p 120.
    -P, --phase	        [default: 200] Seconds between Phases marked on the trace. This is used to avoid overprinting 
                            phase labels. Example -P 30.
    -q, --qscale	        [default: 4.0] Quiver scale for 3C plots. Number of data units per arrow length unit, e.g., 
                            m/s per plot width; a smaller scale parameter makes the arrow longer.
    -r, --region	        [default: na] The region to create the GMV for. The selected region must be a key of the 
                            proj_regions. The acceptable regions are: ['af', 'ak', 'am', 'as', 'au', 'aua', 'eur', 'gl', 
                            'np', 'na', 'sa', 'us']).
    -s, --sizem	        [default: 4.0] Marker (https://matplotlib.org/api/markers_api.html) size in points.
    -S, --rsta		[default: MFID] Reference station code to plot.
    -t, --etime	        [*required] The event time as YYYY-MM-DDTHH:MM:SS.
    -T, --title	        [default: based on the event magnitude (-m) and location (-e)] Double-quoted GMV title.
    -z, --depth	        [*required] Event depth in km.
 
For additional configuration, please see the parameter file (gmv_param.py)


EXAMPLES:\
   _Note: GMV production is very involved. Depending on the duration, number of stations  and number of
          channels, a regular GMV production may take between one to two hours. Three-component GMVs may more than
          three hours to complete._

	1. Sample request with the least number of arguments:
		gmv_generalized.py -e 55.1046,-158.4725 -z 10.0 -m 7.8 -t 2020-07-22T06:12:42 -o GMV_Example_Default

	2. Sample complete request:
		gmv_generalized.py --band=LH,BH --comp=1 -n all -t 2020-07-22T06:12:42 
		-T "July 22, 2020, Alaska Peninsula, M 7.8" -m 7.8 -z 10.0 -e 55.1046,-158.4725 -r ak -d 1200 -s 6.0 -p -180 
		-q 3.5 -g 3 -D 0.05 -G -o GMV_Example_Custom

	3. Sample 3C request:
		gmv_generalized.py --band=LH,BH --comp=3 -n all -t 2020-06-18T12:49:53 
		-T "June 18, 2020, South Of Kermadec Islands, M 7.4" -m 7.4 -z 10.0 -e -33.29,-177.84 -r ak 
		-d 4129.0 -s 6.0 -p 466.0 -q 2.5 -g 5 -D 0.05 -G -o GMV_Example_3C
		
SELECTING PARAMETERS:

 GMVs are event-based animation and therefore you need to start with the event parameter. As a minimum, you need to 
 provide the fallowing event parameters:
 
 - location in the form of _latitude,longitude_ using the _-e_ option (-e 55.1046,-158.4725)
 - depth in kilometers using the _-z_ option (-z 10.0)
 - magnitude using the _-m_ option (-m 7.8). The magnitude lets the script select an optimum filter band to create GMVs
 - time using the _-t_ option (-t 2020-07-22T06:12:42). Time must follow the  yyyy-mm-ssThh24:mm:ss format and must be 
 in UTC time
 
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
 - Sometimes it is desirable to limit the number of networks contributing to the animation. The _-n_ allows you to limit data requests only
   to the list you provide (separate multiple networks with comma). _-n all_ will request data from all available networks. When tuning the GMV,
   use a network with fewer number of stations to allow speedy processing. Once you are happy with the parameters, set the network list to _all_ or desired list.
 - You can select a reference station for display its trace below GMV panel. _-N_ is used to indicate the network of the reference station (see also the _-S_ option below).
 - You can provide the animation file name using the _-o_ option.
 - To avoid long initial video segment without major phase arrivals, use '_-p_' option to indicate how many seconds  after event origin time the video should start.
 - To avoid overprinting phase labels. use the _-P_ option to indicate the number of seconds between Phases marks on the trace. 
 - Quivers are used to depict horizontal motions when creating 3-component animations. The _-q_ option provides the number of data units per arrow 
   length unit, e.g., m/s per plot width; a smaller scale parameter makes the arrow longer.
 - To indicate the region for which the GMV should be created, use the _-r_ option. The selected region must be a key of the _proj\_regions_ dictionary in the param file. 
   See the parameter file for definition of the existing regions.
   The acceptable regions are: dict_keys(['af', 'ak', 'am', 'as', 'au', 'aua', 'eur', 'gl', 'np', 'na', 'sa', 'us']. You can define your own region by adding it to
   the _proj\_regions_ dictionary.
 - Marker size for stations should be selected proportional to the region. Selecting very small markers will make it difficult to correlate station marker colors  and 
   very large marker size will dominate animation and  distract the viewer. The _proj\_regions_ dictionary provides a reasonable marker size for each region. Howerev, you
   can override this value by using the _-s_ option that provides the marker size in points.
 - Reference station code to plot is provided using the _-S_ option. The station code along with the network code _-N_ option are used for selection of the reference station trace that is plotted below the animation pane.
 - To override the default animation title that is based on the event magnitude (_-m_) and location (_-e_), use the _-T_ option. This option takes a double-quoted GMV title as its value.

NOTES:

- GMV generation, depending on the duration, number of stations  and number of channels, may take between one to two hours. Three-component GMVs may take as much as three hours to complete. Always start with a small duration for GMV and a small number of station to create GMV quickly. Tune the production parameters and then go for the actual production run.
          
- GMV production often requires many individual traces. To avoid  making one single large request for data, it is better to break data request to data centers to make multiple smaller requests. The parameter _chunck_count_ in the parameter file controls the number of 
stations per request (chunk) that will be sent to each data center. This number should be adjusted based on the number of station-channels involved 
and the length of each request to avoid memory issues or data center timeouts.


CITATION:

To cite the use of this software reference:\
Trabant, C., A. R. Hutko, M. Bahavar, R. Karstens, T. Ahern, and R. Aster (2012), Data Products at the IRIS DMC: Stepping Stones for Research and Other Applications, Seismological Research Letters, 83(5), 846–854, https://doi.org/10.1785/0220120032.

Or cite the following DOI:\
    doi:10.17611/dp/gmv.code.1

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


