#Tested in R 3.6.3 ("Holding the Windsock") in RStudio 1.2.1335
#This example plots a time-series of wind data from the Pioneer Offshore surface buoy for the day of June 01, 2019.

#Installation
#install.packages("httr")
#install.packages("ncdf4")
#install.packages("devtools")
require(devtools)
install_github("oceanobservatories/ooim2mr")


#Required packages.
require(ooim2mr)
require(ncdf4)
require(jsonlite)
require(httr)
require(stringr)


###----------CHANGE ME---------------###
user = "YOUR-OOI-API-USER-HERE"   #OOI API Username
token = "YOUR-OOI-API-TOKEN-HERE"  #OOI API Token

site = 'CP04OSSM'
node = 'SURFACE'
instrument = 'MET'
method = 'recovered_host'
stream = "metbk_a_dcl_instrument"
start_date = '2019-06-01'
stop_date = '2019-06-01'
drop_paired = TRUE
simplify_data = TRUE
###---------------------------------###


#Generate a URL, keeping start and stop times as defaults.
url = ooi_create_url(site = site,
                     node = node,
                     instrument = instrument,
                     method = method,
                     stream = stream,
                     start_date = start_date,
                     stop_date = stop_date) 

#Issue the request.
response = ooi_submit_request(url,user = user, token = token)  

#Get the OpenDAP URLs.
remote = ooi_get_location(response,drop_paired = drop_paired)  

#Download the data.
local = ooi_download_data(remote) 

#Read the data in.
lol = ooi_get_data(local,simplify_data = simplify_data) 

data = data.frame(lol[['data']]) #Data is always in the first sublist.
varunits = data.frame(lol[['variables_units']])  #Variables and associated units are always the second sublist.
data$time = as.POSIXct(data$time,tz='UTC',origin='1900-01-01')  #Convert time to something understandable. Origin can be found in the units of time.

#Quick plot of chl-a data.
plot(data$time,data$met_relwind_speed,main = "CP04OSSM Wind Speed",ylab = 'Rel. Win Speed (m/s)',xlab = 'Time (UTC)')