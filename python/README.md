# OOI Data Explorations with Python

The python code provided here was developed primarily as a toolset for myself, as an OOI Data Team member, to 
facilitate accessing data from OOINet for the QC reviews, gap analyses, metadata checks, etc that I need to perform as 
part of my day-to-day work. I am providing this code to the larger community in the hopes that it will be of some use 
to you, since the steps I use for accessing and using the data are identical to those any user would employ. The code 
uses the [OOI M2M API](https://oceanobservatories.org/ooi-m2m-interface/) to access the data and either loads it into
the user workspace as an [xarray](http://xarray.pydata.org/en/stable/) dataset, or saves it to disk as a NetCDF file, 
depending on how you call it. There are examples below of how to setup and use the package, with more example notebooks 
and scripts available in the examples directory.

## Configuring System for Python (Install Bash, Git and Anaconda/Miniconda)

In order to use the python code in this repository, you will need to setup your computer with the proper tools. There
are several examples out there on how to this, so I'll avoid reinventing the wheel here. One of the best 
tutorials I've found has been developed by the folks at [Earth Lab](https://www.earthdatascience.org/). The 
[tutorial](https://www.earthdatascience.org/workshops/setup-earth-analytics-python/setup-git-bash-conda/) they have 
prepared will guide you through the process of setting up a system to use Python for Earth Science analysis from start
to finish, regardless of your computer's operating system. Experienced users can easily skip and skim through to the 
sections relevant to them.

One key difference between the tutorial and my personal preference is to use the full [Anaconda](https://www.anaconda.com) 
installation instead of [Miniconda](https://docs.conda.io/en/latest/miniconda.html) as called out in the tutorial. All 
the tools you would need, and then some, are installed with Anaconda. You can follow the tutorial exactly as written 
and all of the code in this project will work, but I recommend using Anaconda instead of Miniconda. There are several 
other things you may end up wanting to do with the tools available in Anaconda.

Note, for Windows users only, if you already have Anaconda/Miniconda installed on your machine, you do not need to 
uninstall/reinstall as described in the [tutorial](https://www.earthdatascience.org/workshops/setup-earth-analytics-python/setup-git-bash-conda/).
You can leave everything as-is. However, you do need to link Git Bash to Anaconda (or Miniconda); this happens
automagically if you follow the sequence in the tutorial by installing Git Bash before Anaconda. If you already have 
Anaconda installed, however, from the bash terminal add the following code to your `.bash_profile` file in your home 
directory (assuming you installed Anaconda in your home directory, which is the default):

```shell script
cd ~
echo ". ${HOME}/Anaconda3/etc/profile.d/conda.sh" >> ~/.bash_profile
source .bash_profile
```

### Obtaining the Code and Configuring the Environment

If you followed the tutorial above to configure Bash, Git and Anaconda on your system, the next few steps should 
proceed easily. With the basic processing tools in place, now you need to copy the python code to the local machine, 
setup a virtual environment, install a development copy of the code, and get access credentials in place.

From the bash terminal, clone the ooi-data-explorations code to your local machine:

```shell script
# download the ooi-data-explorations code
mkdir -p ~/Documents/GitHub
cd ~/Documents/GitHub
git clone https://github.com/oceanobservatories/ooi-data-explorations.git
cd ooi-data-explorations/python

# configure the OOI python environment
conda env create -f environment.yml
conda activate ooi

# install the package as a local development package
pip install -e . 
```

### Access Credentials

Access credentials are required to download the data from OOINet via the M2M interface. Directions on how to obtain 
these, in addition to details about the M2M system, are [available here](https://oceanobservatories.org/ooi-m2m-interface/).

* If you haven't already done so, either create a user account on the [OOI Data Portal](https://ooinet.oceanobservatories.org), 
or use the CILogon button with an academic or Google account (login button is towards the upper right corner of the 
web page) to login to the portal.
* Navigate to the drop down menu screen in the top-right corner of the menu bar 
* Click on the “User Profile” element of the drop down.
* Copy and save the following data from the user profile: API Username and API Token.

The python code uses the [netrc](https://docs.python.org/3.6/library/netrc.html) utility to obtain your access 
credentials. Users need to create a `.netrc` file in their home directory to store these access credentials. From the 
bash terminal (replacing the `<API Username>` and `<API Token>` below with your credentials from the 
[OOI Data Portal](https://ooinet.oceanobservatories.org)):

```shell script
cd ~
touch .netrc
chmod 600 .netrc
cat <<EOT >> .netrc
machine ooinet.oceanobservatories.org
    login <API Username>
    password <API Token>
EOT
```

## Using the Tools

The code is available in the [ooi_data_explorations](ooi_data_explorations/) directory with examples (both scripts
and notebooks) in the [examples](examples/) directory. The python code has been developed and used with both Python 3.6
and 3.7 on Windows and Linux machines. They are configured in a granular fashion to allow users to access the data in
a few different ways.

Before using these functions, it is important to understand how requests to the 
[OOI M2M API](https://oceanobservatories.org/ooi-m2m-interface/) are structured. A request is built around the 
reference designator (comprised of the site, node and sensor names), the data delivery method, and data stream (think 
of a stream as a data sets). Beginning and ending dates for the time period of interest are optional inputs. If omitted,
all of the data for a particular instrument of interest will be downloaded.

* Site -- 8 character uppercase string denoting the array and location within the array of the system. These are 
[defined](https://oceanobservatories.org/research-arrays/) on the OOI website.
* Node -- 5 character uppercase string (of which the first 2 characters are really the key) denoting the assembly the
instrument is connected to/mounted on. These can be thought of as physical locations within/under the higher level site 
designator.
* Sensor -- 12 character uppercase string that indicates, among other things, the instrument class and series. The class
and series [defined](https://oceanobservatories.org/instruments/) on the OOI website.
* Delivery Method -- Method of data delivery (lowercase).
    * `streamed` -- Real-time data delivery method for all cabled assets. Data is "streamed" to shore over the fiber
    optic network as it outputs from an instrument. 
    * `telemetered` -- Near real-time data delivery method for most uncabled assets. Data is recorded remotely
    by a data logger system and delivered in batches over a satellite or cellular network link on a recurring schedule 
    (e.g every 2 hours).
    * `recovered_host` -- Usually the same data set as the telemeterd for uncabled assets. Key difference is this data is 
	downloaded from the data logger system after the asset is recovered. In most cases, this is 1:1 with the 
    telemetered data unless there was an issue with telemetry during the deployment or the data was decimated
    (temporal and/or # of parameters) by the data logger system  prior to transmission.
    * `recovered_inst` -- Data recorded on and downloaded directly from an individual instrument after the instrument is
    recovered. Not all instruments internally record data, so this method will not be available for all instruments.
	* `recovered_wfp` -- Data recorded on and downloaded from the McLane Moored Profiler system used at several sites 
	in OOI. Telemetered data is decimated, this data set represents the full-resolution data.
	* `recovered_cspp` -- Data recorded on and downloaded from the Coastal Surface Piercing Profiler system used in
	the Endurance array. Telemetered data is decimated, this data set represents the full-resolution data.
* Stream -- A collection of parameters output by an instrument, or read from a file, and parsed into a named data set.
Stream names are all lowercase. Streams are mostly associated with the data delivery methods and there may be more than
one stream per method.

### Request Unmodified (Mostly, Kinda) Data

The core functions used to request and download data are `m2m_request` and `m2m_collect`, located in the 
[`common.py`](ooi_data_explorations/common.py) module. From those two functions, you can pretty much create your own 
library of functions to download and process whatever data you want from the system and either save it locally or
continue to work on it within your python environment. It is important to note, these functions require inputs that
map directly to those required by the [OOI M2M API](https://oceanobservatories.org/ooi-m2m-interface/).

The data requested and downloaded by the `m2m_request` and `m2m_collect` functions is somewhat modified from the
original: I switch the dimensions from `obs` to `time`, drop certain timing variables that were originally never meant 
to be exposed to the user, and clean up some basic metadata attributes. Beyond that, the data is provided as-is. No 
effort is made to select a subset of the variables, conduct QC, clean up the variable names, or in any other way alter 
the data obtained from the [OOI Data Portal](https://ooinet.oceanobservatories.org). For example:

```python
import os
from  ooi_data_explorations.common import m2m_request, m2m_collect

# Setup the needed information to request data from the pH sensor on the Oregon
# Shelf Surface Mooring near-surface (7 m depth) instrument frame (NSIF).
site = 'CE02SHSM'           # OOI Net site designator
node = 'RID26'              # OOI Net node designator
sensor = '06-PHSEND000'     # OOI Net sensor designator
method = 'telemetered'      # OOI Net data delivery method
stream = 'phsen_abcdef_dcl_instrument'  # OOI Net stream name
start = '2019-04-01T00:00:00.000Z'  # data for spring 2019 ...
stop = '2019-09-30T23:59:59.999Z'   # ... through the beginning of fall

# Request the data (this may take some time).
r = m2m_request(site, node, sensor, method, stream, start, stop)

# Use a regex tag to download only the pH sensor data from the THREDDS catalog
# created by our request.
data = m2m_collect(r, '.*PHSEN.*\\.nc$')

# Save the data to the users home directory under a folder called ooidata for
# further processing
out_path = os.path.join(os.path.expanduser('~'), 'ooidata')
out_path = os.path.abspath(out_path)
if not os.path.exists(out_path):
    os.makedirs(out_path)

# setup the output file
out_file = ('%s.%s.%s.%s.%s.nc' % (site.lower(), node.lower(), sensor.lower(), method, stream))
nc_out = os.path.join(out_path, out_file)

# save the data to disk
data.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf')
```

### Request Unmodified, Somewhat De-Obfuscated Data

Users really only need to use `m2m_request` and `m2m_collect` for the data requests. However, the user needs to
explicitly know all of the details and terms from above. Outside of OOI (and even inside of OOI), the terminology used
for sites, nodes, sensors, methods, and streams can be intimidating and confusing. In an attempt to clean up some of 
that terminology, limit the need for the user to learn all of the OOI lingo, and to align with some of the 
functionality in the Matlab and R utilities, I've organized a subset of all of the sources of OOI data into a 
[YAML structure](ooi_data_explorations/m2m_urls.yml) that users can query using a simpler set of terms as part of a 
data request. The YAML structure removes from consideration all so-called engineering sensors and instruments that 
cannot be accessed through the API (e.g. cameras or bio-accoustic sonar), as well as most of the non-science streams 
(engineering or metadata streams). The idea is to cover the most common needs, rather than all possible cases. Users can
still access any data set desired, but you will need to use the example from above to explicitly call any stream not
represented in the YAML structure.

The YAML structure I've created uses the OOI site codes as-is. It simplifies the node designation by taking the 100+
nodes and groups them according to an assembly type indicating where co-located sensors can be found. There are 6 
assembly types with differing numbers of subassemblies (see table below). Either the assembly or subassembly name can 
be used to request the data.

| Assembly | Subassembly | Description |
| --- | --- | --- |
| buoy | n/a | Surface buoys with meteorological, wave, and/or sea surface (~1 m depth) instrumentation |   
| midwater | nsif, riser, sphere, 200m_platform | Platforms located at various depths below the sea surface and above the seafloor |   
| seafloor |mfn, bep, low-pwr-jbox, medium-pwr-jbox| Platforms resting directly on the seafloor |   
| profiler | cspp, coastal-wfp, global-wfp, shallow-profiler, deep-profiler | Profiling systems with integrated instrumentation |   
| glider | coastal-glider, global-glider, profiling-glider | Autonomous, buoyancy driven underwater vehicles with integrated instrumentation |   
| auv | n/a | Autonomous underwater vehicles, currently only deployed in the Pioneer Array  |   

The shorter OOI instrument class name is used instead of the full sensor designator. The data delivery methods are as
defined above and are used to determine the stream(s) available. There is usually just one. In the few cases where 
there is more than one, the code defaults to selecting the first. I've curated the list to make sure this is the stream 
of interest to 99.9% of users. The key utility here is users do not have to know the stream name. You can still get at 
the other streams, if needed, but you have to explicitly know what they are and call them as shown in example above.

The last things to consider are date ranges to bound the request and an aggregation value. Date ranges are fairly 
self-explanatory. You need to select a starting and ending date for the data of interest, otherwise you will get all of 
the data for that instrument. That could potentially be a large request, so be careful. The dates entered need to be
recognizable as such. I'm using the [dateutil parser function](https://dateutil.readthedocs.io/en/stable/parser.html) to 
convert the dates you enter into the proper form for the M2M API. Alternatively, you can use the deployment number and
the M2M API will determine the dates. 

The aggregation value addresses a few cases where more than one instance of an instrument class is associated with
an assembly. For example the Global surface moorings have a midwater chain of 10 CTDs connected to an inductive modem 
line. If you do not specify the aggregation flag, you will only get the first of those 10 CTDs. If, however, you set the
aggregation flag to `0`, you will get all of them in a data set with a variable added called `sensor_count`, so you can
distinguish between them. Conversely, you can request a specific instrument by using its sequential number.

At the end of the day, users need to know the site, assembly, instrument class and data delivery method. Somewhat
simpler than the full site, node, sensor, method, and stream name, and hopefully more meaningful. Additionally, any
date/time string format, so long as it can be recognized by the 
[dateutil parser function](https://dateutil.readthedocs.io/en/stable/parser.html), will work for setting starting and
ending dates for the data requests as opposed to explicitly setting the date to a format of `YYYY-MM-DDThh:mm:ss.fffZ`.
Or, you can just use the deployment number and let the system figure out the dates. Take a look at the examples below: 

```python
from ooi_data_explorations.data_request import data_request

# Setup the needed information to request data from the pH sensor on the Oregon
# Shelf Surface Mooring near-surface (7 m depth) instrument frame (NSIF).
site = 'ce02shsm'           # OOI site designator
assembly = 'midwater'       # Assembly grouping name
instrument = 'phsen'        # OOI instrument class 
method = 'telemetered'      # data delivery method

# the first four inputs are required in the order given above, the following
# inputs are semi optional, You need to specify at least a start date and/or
# stop date, or use the deployment number
start = '2019-04-01'       # data for April 2019 ...
stop = '2019-05-30'        # ... through May 2019
deploy = 9                 # The Spring 2019 deployment number

# request and download the data using specific dates
data_01 = data_request(site, assembly, instrument, method, start=start, stop=stop)

# request and download the data using the deployment number
data_02 = data_request(site, assembly, instrument, method, deploy=deploy)

# Setup the needed information to request data from the CTDMO sensor on the 
# Global Irminger Surface Mooring inductive modem line.
site = 'gi01sumo'           # OOI site designator
assembly = 'riser'          # Subassembly grouping name
instrument = 'ctdmo'        # OOI instrument class 
method = 'recovered_inst'   # data delivery method

start = '2019-04-01'        # data for April 2019 ...
stop = '2019-05-30'         # ... through May 2019

# request and download the data using specific dates, this only returns the 
# first instance of the CTDMOs
data_03 = data_request(site, assembly, instrument, method, start=start, stop=stop)

# request and download the data for all 10 CTDMOs
data_04 = data_request(site, assembly, instrument, method, start=start, stop=stop, aggregate=0)

# request and download the data for the CTDMO 5 out of 10.
data_05 = data_request(site, assembly, instrument, method, start=start, stop=stop, aggregate=5)
```

### Request Processed and Modified Data

For most individuals, the above code should satisfy your needs. For some, like myself, the data simply is not
properly organized, named or fully processed enough to be really useful. The process_*.py modules in the cabled and
uncabled directories represent an attempt on my part to rework the data sets into more useful forms before conducting
any further work. Primarily, these re-works are for my own use, but they are available for others to use. The primary
steps are:

* Deleting certain variables that are of limited or no use to the actual data set.
* Renaming the alphabet soup parameter names to more meaningful, user friendly names (the original OOI names are
preserved as a variable level attribute termed `ooinet_variable_name`).
* Resetting the QC parameters to use the `flag_mask` and `flag_meaning` attributes from the 
[CF Metadata conventions](http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html#flags).
* Reseting incorrectly set units and other attributes for some of the variables.
* Reworking certain parameters by splitting or reshaping the data into more useful forms. 
* Update global attributes and otherwise clean-up the data sets.

Additionally, some of the instrument data is collected in bursts (e.g. every 15 minutes for 3 minutes at 1 Hz). This
can make some of the data sets fairly large. By applying a median average to each of the bursts, the size of the data
set can be reduced to a more workable form, and the point-to-point variability in each burst can be smoothed. Burst 
averaging is optional. Most of these functions are set to run from the command line. Examples of how these are run
can be found in the [examples directory](examples/). Bash scripts to automate downloading date using these processing
scripts are in the [utilities/harvesters](utilities/harvesters) directory.
