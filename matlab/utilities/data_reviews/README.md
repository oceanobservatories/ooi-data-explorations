# OOI Data Review Tool
Uses a Matlab App to request and load data downloaded from the OOI system for 
review, with the ability to either download associated annotations (after the 
site, node and stream are specified), or load locally saved annotations. Basic
workflow is to load locally saved data, plot selected variables, review, create,
or update annotations, save those annotations locally and repeat. Can be used
to formally review the data for a particular instrument (e.g. download and save 
all the PHSEN data for the CE01ISSM NSIF and create, update, or correct 
annotations associated with the different data types such as recovered_host or 
recovered_inst) or just to quickly peruse datasets to get a sense of overall
performance.

## Usage
`TODO`

## Assumptions
Developed and tested with Matlab 2019b and 2020a on a PC (Windows 10) and a Mac
(Mojave, 10.14.6). Earlier versions of Matlab will not be supported.

### Access Credentials for the OOI M2M API
In order to access the OOI M2M API, you need to setup your OOINet credentials as
a `weboptions` object that Matlab can then use as part of a `webread` request
to the API.  To set this up, in Matlab (replacing `<API Username>` and 
`<API Token>` with your OOINet credentials):

``` matlab
>> username = '<API Username>';
>> password = '<API Token>';
>> options = weboptions('Timeout', 180, 'HeaderFields', {'Authorization', ...
                 ['Basic ' matlab.net.base64encode([username ':' password])]});
>> save('ooinet.credentials.mat', 'options');
```

Put the `ooinet.credentials.mat` file in your Matlab path (best option is your
Matlab home directory as it is on the Matlab path by default) so it is available
for this application.  The GUI uses these credentials to put together the 
hierarchical data request structure and to pull annotations.

### Setup for M2M Data Requests
To request data via the M2M system (after structuring the request using the
hierarchical dropdowns to bound the request) you will need to have the python
code available in the [ooi-data-explorations](https://github.com/oceanobservatories/ooi-data-explorations/tree/master/python) 
repository installed on your system per the directions in the 
[README](https://github.com/oceanobservatories/ooi-data-explorations/blob/master/python/README.md). 

Matlab does not support Anaconda or the use of conda to change the python
environment, so it needs to be explicitly setup to use the correct code. To 
determine the correct path elements to add, from a separte terminal, run the 
following commands to snapshot the system path before and after activating the
ooi python environment:

``` bash
# If using Windows, use %PATH% for the path variable name rather than $PATH
echo $PATH
conda activate ooi
echo $PATH
```
The elements that are prepended to the path after activating the ooi environment
are the ones to note/copy. I use a `startup.m` file in my Matlab home directory 
to then setup the path and the python environment. 

``` matlab
% startup.m
username = getenv('username');
setenv('PYTHONUNBUFFERED', '1');
setenv('PATH', ['C:\Users\' username '\Anaconda3\envs\ooi;', ...
    'C:\Users\' username '\Anaconda3\envs\ooi\Library\mingw-64\bin;', ...
    'C:\Users\' username '\Anaconda3\envs\ooi\Library\usr\bin;', ...
    'C:\Users\' username '\Anaconda3\envs\ooi\Library\bin;', ...
    'C:\Users\' username '\Anaconda3\envs\ooi\Scripts;', ...
    'C:\Users\' username '\Anaconda3\envs\ooi\bin;', ...
    getenv('PATH')])
[~] = pyenv('Version',['C:\Users\' username '\Anaconda3\envs\ooi\pythonw.exe'], ...
            'ExecutionMode','InProcess');
clear username
```

The above works on my PC. Mac/Linux users will need to adjust accordingly. Note,
Windows used a `;` to separate paths, while the path separator is a `:` for
Mac/Linux. For example, the above `startup.m` becomes the following on a 
co-worker's Mac:

``` matlab
% startup.m
setenv('PYTHONUNBUFFERED', '1');
setenv('PATH', ['/Applications/anaconda3/envs/ooi/bin:'  getenv('PATH')])
[~] = pyenv('Version', '/Applications/anaconda3/envs/ooi/bin/python',...
  'ExecutionMode', 'InProcess');
```

Everyone's system will be somewhat different. You may need to play around with
the above to make sure you have all the path elements set correctly.

### Testing
To confirm that you have setup Matlab correctly, you can run the following set
of commands:
``` matlab
>> % test that the system paths match what was set in startup.m via setenv
>> [~, ~] = system('python -m site', '-echo');
sys.path = [ 
    'C:\\Users\\cwingard\\Documents\\GitHub\\ooicgsn-data-tools\\data_reviews', 
    'C:\\Users\\cwingard\\Anaconda3\\envs\\ooi\\python37.zip', 
    'C:\\Users\\cwingard\\Anaconda3\\envs\\ooi\\DLLs', 
    'C:\\Users\\cwingard\\Anaconda3\\envs\\ooi\\lib', 
    'C:\\Users\\cwingard\\Anaconda3\\envs\\ooi', 
    'C:\\Users\\cwingard\\Anaconda3\\envs\\ooi\\lib\\site-packages', 
    'c:\\users\\cwingard\\documents\\github\\ooi-data-explorations\\python', 
    'C:\\Users\\cwingard\\Anaconda3\\envs\\ooi\\lib\\site-packages\\win32', 
    'C:\\Users\\cwingard\\Anaconda3\\envs\\ooi\\lib\\site-packages\\win32\\lib', 
    'C:\\Users\\cwingard\\Anaconda3\\envs\\ooi\\lib\\site-packages\\Pythonwin', 
] 
USER_BASE: 'C:\\Users\\cwingard\\AppData\\Roaming\\Python' (doesn't exist) 
USER_SITE: 'C:\\Users\\cwingard\\AppData\\Roaming\\Python\\Python37\\site-packages' (doesn't exist) 
ENABLE_USER_SITE: True 
>> 
>> % test that we can run python
>> py.math.exp(1)

ans =

    2.7183

>> % test that we can load a module from our python environment and run it
>> np = py.importlib.import_module('numpy');
>> np.exp(1)

ans =

    2.7183

>> % test that we can load and execute a command from ooi_data_explorations
>> od = py.importlib.import_module('ooi_data_explorations.common');
>> od.list_nodes('CE01ISSM')

ans = 

  Python list with no properties.

    ['MFC31', 'MFD35', 'MFD37', 'RID16', 'SBC11', 'SBD17']
```
 
If all of the above works and the paths match what was set in your `startup.m`,
then you should be good to use `data_reviews`. Please note, however, that this
is still very much a beta product. There will be issues. Don't hesitate to let
me know by raising an [issue](https://github.com/oceanobservatories/ooicgsn-data-tools/issues).
