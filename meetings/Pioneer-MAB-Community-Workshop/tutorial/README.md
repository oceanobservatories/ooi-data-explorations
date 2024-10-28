# pyOOI
#### Author: Andrew Reed

---
The modules and tools within this repo are designed to assist in requesting, importing, downloading, and processing data from the Ocean Observatories Initiative API by M2M requests. This is a streamlined package based on https://github.com/reedan88/OOINet while fusing important functionality from https://github.com/oceanobservatories/ooi-data-explorations/tree/master/python.

The Ocean Observatories Initiative (OOI) is a science-driven ocean observing network that delivers real-time data from more than 800 instruments to address critical science questions regarding the world’s ocean. OOI data are freely available online to anyone with an Internet connection. OOI features several methods accessing data: via the web portal at ooinet.oceanobservatories.org or via OOINet's API.

The pyOOI package is designed to simplify and streamline the searching, requesting, and downloading of datasets from OOINet via the API. It is designed for interactive computing in jupyter notebooks for the science end-user. It allows for searching of specific datasets based on the location and/or specific instrument and built-in functions reprocess datasets to allow for multiple overlapping netCDF files to be loaded remotely into an xarray dataset. It will also retrieve associated metadata and vocabulary for specific datasets.

---
## Setup
Please follow the setup outlined in the repo's [parent directory README](../).

---
## Modules

### M2M
This module consists of a number of functions which simplify interacting with OOINet's API. These include functions to search for datasets, retrieve metadata, deployment information, and instrument vocab. It also includes tools for downloading netCDF datasets from the OOI THREDDS and Gold Copy THREDDS servers.

### Bottles
This module consists of a number of functions to process, interpret, and clean up the discrete sample summary spreadsheets from OOI for use as validation against deployed instrumentation.