% OOINET_DEFAULTS create global variables used by the different functions

% Set the base URL for the M2M interface
BASE_URL = 'https://ooinet.oceanobservatories.org/api/m2m/';  % Base M2M URL
try
    % Check to see if we are on the internal network, if so change the BASE_URL
    intra = webread('https://ooinet.intra.oceanobservatories.org', weboptions('Timeout', 1));
    BASE_URL = 'https://ooinet.intra.oceanobservatories.org/api/m2m/';  % Base intranet M2M URL
    INTRANET = true;
catch
    INTRANET = false;
end %try

% set the sub-URLs for the different API endpoints
ANNO_URL = '12580/anno/';                     % Annotation Information
ASSET_URL = '12587/asset/';                   % Asset and Calibration Information
DEPLOY_URL = '12587/events/deployment/inv/';  % Deployment Information
SENSOR_URL = '12576/sensor/inv/';             % Sensor Information
VOCAB_URL = '12586/vocab/inv/';               % Vocabulary Information
STREAM_URL = '12575/stream/byname/';          % Stream Information
PARAMETER_URL = '12575/parameter/';           % Parameter Information

% other constants used in the code
FILL_INT = -9999999;
FILL_FLOAT = NaN;

% load the access credentials
try
    load('ooinet.credentials.mat')  % returns a variable called options
catch
    error(['Unable to load access credentials. Users need to create a ' ...
           'weboptions object with their personal OOINet API keys.'])
end %try
