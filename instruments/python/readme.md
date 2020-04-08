## Configuring System for Python

Install Anaconda

Setup the ooi environment

## Access Credentials

Access credentials are required to download the data from OOINet via the M2M 
interface. Directions on how to obtain those, in addition to details about the
M2M system are [available](https://oceanobservatories.org/ooi-m2m-interface/).

* If you haven't already done so, either create a user account on the 
  [OOI Data Portal](https://ooinet.oceanobservatories.org), or use the CILogon
  button with an academic or Google account (login button is towards the upper 
  right corner of the web page).
* Login to your account.
* Navigate to the drop down menu screen in the top-right corner of the menu bar
* Click on the “User Profile” element of the drop down.
* Copy and save the following data from the user profile: API Username and API 
  Token.

The python code uses the netrc utility to obtain the access credentials . Users
need to create a `.netrc` file in their home directory to store their access
credentials to download the data from the OOI system. Add the following to your .netrc file:

machine ooinet.oceanobservatories.org
    login <API Username>
    password <API Token>

## Using the Tools

The python tools have been developed using both 3.6 and 3.7 on windows and linux machines, and are configured in a
granular fashion to allow users to access the data in a few different ways. The core functions needed are all in the
common module with configuration options available in a YAML file (rename 'config.yaml.changeme' to 'config.yaml') that
the user can customize to specify local
system paths for saving the resulting data if desired.

The structure of the python code differs from the Matlab in order to take advantage of the core functionality provided
by python and more specifically xarray. Also, do not introduce new terminology. While the OOI terminology can be obtuse
at first, creating another layer simply furthers the confusion.

### Request Unmodified (Mostly) Data

These modules request data for the instrument (from select streams), switches the dimensions from `obs` to `time`, and
drops certain timing variables that were originally never meant to be exposed to the user. These modules are essentially
wrapper functions around three key functions that form the core of interacting with the M2M system: requesting the data,
downloading (or collecting the data) from the THREDDS catalog created to hold the results of the request, and saving the
data to the local hard drive for further work. For example:

```python
import os
from instruments.python.common import m2m_request, m2m_collect, CONFIG

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

Of the above functions, users really only need `m2m_request` and `m2m_collect`. From those core functions, users can
pretty much create their own library of functions to download whatever data they want from the system.m2m_collect.

### Request Processed, Modified Data

request data, clean up variable names, remove variables that don't belong, fix attributes, apply burst averaging if so
desired

process_phsen

