# OOI Data Explorations with Python

## Configuring System for Python

For the purposes of this README and the examples below, I will assume your
user name is `ooiuser`. Make sure you change it accordingly.

### Install Bash, Git and Anaconda/Miniconda

In order to use the python tools in this repository, you will need to setup
your computer with the proper tools. There are several examples out there on
how to this, so we'll avoid reinventing the wheel here. One of the best 
[tutorials](https://www.earthdatascience.org/workshops/setup-earth-analytics-python/setup-git-bash-conda/) 
we've found has been developed by the folks at [Earth Lab](https://www.earthdatascience.org/).
The [tutorial](https://www.earthdatascience.org/workshops/setup-earth-analytics-python/setup-git-bash-conda/)
they have prepared will guide you through the process of setting up a system 
to use Python for Earth Science analysis from start to finish, regardless of 
your computer's operating system. Experienced users can easily skip and skim 
to the sections relevant to them.

One key difference between the tutorial and our personal preference is to use 
the full Anaconda installation. All the tools you would need, and then some, 
are installed. From here on out, however, we'll assume you installed Miniconda.

_**Note, for Windows users only**_

If you already have Miniconda installed on your machine, you do not need to
uninstall/reinstall as described in the [tutorial](https://www.earthdatascience.org/workshops/setup-earth-analytics-python/setup-git-bash-conda/).
You can leave everything as-is. However, you do need to link Git Bash to 
Miniconda (or Anaconda); this automagically happens if you follow the tutorial. Add
the following code to your `.bashrc` or `.bash_profile` file in your home directory:

```shell script
cd /c/Users/ooiuser

```

### Obtaining the Code and Configuring the Environment

If you followed the tutorial above to configure Bash, Git and Miniconda on your
system, the next few steps should proceed easily. If not, I'll assume you are 
an advanced user and can piece most of this together yourself. With the basic 
processing tools in place, now we need to copy the code to the local machine, 
setup the local environment and run some basic tests to make sure everything 
is working.

From the Git Bash terminal, clone the ooi-data-explorations code to your local
machine:

```shell script
# download the ooi-data-explorations code
cd /c/Users/ooiuser/Documents/GitHub
git clone https://github.com/oceanobservatories/ooi-data-explorations.git
cd ooi-data-explorations/python

# configure the OOI python environment
conda create env -f environment.yml
conda activate ooi

# run the basic tests to make sure the key pieces are working
nosetests 
```

git clone

conda create env -f environment.yml

### Access Credentials

Access credentials are required to download the data from OOINet via the M2M 
interface. Directions on how to obtain these, in addition to details about the
M2M system, are [available](https://oceanobservatories.org/ooi-m2m-interface/).

* If you haven't already done so, either create a user account on the 
  [OOI Data Portal](https://ooinet.oceanobservatories.org), or use the CILogon
  button with an academic or Google account (login button is towards the upper 
  right corner of the web page).
* Login to your account.
* Navigate to the drop down menu screen in the top-right corner of the menu bar
* Click on the “User Profile” element of the drop down.
* Copy and save the following data from the user profile: API Username and API 
  Token.

The python code uses the netrc utility to obtain the users access credentials. 
Users need to create a `.netrc` or `_netrc` file in their home directory to store 
their access credentials. Linux and Mac users need to create a `.netrc` file. 
Windows users need to create a `_netrc` file, however I've found that recent 
versions of the netrc module are looking for a `.netrc` file. That's what works
on two of my Windows 10 machines (with Python and cURL). In all cases, the netrc 
file needs to be placed in the user's home directory (`/home/ooiuser` or 
`c:\Users\ooiuser`) and the permissions need to be changed to read/write 
only for the user.

Add the following to your `.netrc/_netrc` file:

```text
machine ooinet.oceanobservatories.org
    login <API Username>
    password <API Token>
```

We've added default credentials to the code, if you don't want to create your 
own credentials. The downside to using these is you two-fold:

* You won't get an email letting you know when the request is ready (OK,
maybe that isn't a downside).
* You won't get an email letting you know if anything was wrong with the
data that was later corrected. These emails are used to keep you, the user,
informed of any potential issues, corrections, updates to the data and the 
system. 

## Using the Tools

The python tools have been developed and used with both Python 3.6 and 3.7 on
Windows and Linux machines. They are configured in a granular fashion to allow 
users to access the data in a few different ways. The core functions needed are
all in the `common.py` module with configuration options available in a YAML 
file (rename `config.yaml.changeme` to `config.yaml`) that the user can 
customize to specify local system paths for saving the resulting data (if 
desired).

### Request Unmodified (Mostly, Kinda) Data

These modules request data for an instrument (from select streams), switch 
dimensions from `obs` to `time`, and drop certain timing variables that 
were originally never meant to be exposed to the user. Additional a few
metadata attributes are updated. Beyond that, the data is provided as-is.
No effort is made to select the variables or clean up the variable names.
These modules are essentially wrapper functions around three key functions 
that form the core of interacting with the M2M system: requesting the data, 
downloading (or collecting) the data, and then saving the data to disk. 
For example:

```python
import os
from  python.ooi_data_explorations.common import m2m_request, m2m_collect, CONFIG

# Setup the needed information to request data from the pH sensor on the Oregon
# Shelf Surface Mooring near-surface (7 m depth) instrument frame (NSIF).
site = 'CE02SHSM'           # OOI Net site designator
node = 'RID26'              # OOI Net node designator
sensor = '06-PHSEND000'     # OOI Net sensor designator
method = 'telemetered'      # OOI Net data delivery method
stream = 'phsen_abcdef_dcl_instrument'  # OOI Net stream name
level = 'nsif'              # local directory name, level below site
instrmt = 'phsen'           # local directory name, instrument below level
start = '2019-04-01T00:00:00.000Z'  # data for spring 2019 ...
stop = '2019-09-30T23:59:59.999Z'   # ... through the beginning of fall

# Request the data (this may take some time).
r = m2m_request(site, node, sensor, method, stream, start, stop)

# Use a regex tag to download only the pH sensor data from the THREDDS catalog
# created by our request.
data = m2m_collect(r, '.*PHSEN.*\\.nc$')

# Save the data to disk for further processing
out_path = os.path.join(CONFIG['base_dir']['m2m_base'], site.lower(), level, instrmt)
out_path = os.path.abspath(out_path)
if not os.path.exists(out_path):
    os.makedirs(out_path)

out_file = ('%s.%s.%s.%s.%s.nc' % (site.lower(), level, instrmt, method, stream))
nc_out = os.path.join(out_path, out_file)

data.to_netcdf(nc_out, mode='w', format='NETCDF4', engine='h5netcdf')
```

Of the above functions, users really only need `m2m_request` and `m2m_collect`. 
From those core functions, users can pretty much create their own library of 
functions to download whatever data they want from the system and either save 
it locally or continue to work on it within their python environment.

### Request Processed, Modified Data

request data, clean up variable names, remove variables that don't belong, fix attributes, apply burst averaging if so
desired

process_phsen

