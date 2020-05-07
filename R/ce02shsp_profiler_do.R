#Tested in R 3.6.3 ("Holding the Windsock") in RStudio 1.2.1335
#This example plots a time-series of the profiler oxygen optode at CE02SHSP for the year 2019.

#Installation
#install.packages("httr")
#install.packages("ncdf4")
#install.packages("devtools")
require(devtools)
install_github("IanTBlack/ooim2mr")

#Required packages.
require(ooim2mr)
require(ncdf4)
require(jsonlite)
require(httr)
require(stringr)

#Packages for plotting.
require(ggplot2)
require(RColorBrewer)

###----------CHANGE ME---------------###
user = "OOIAPI-BCJPAYP2KUVXFX"   #OOI API Username
token = "D3HV2X0XH1O"  #OOI API Token
site = 'CE02SHSP'
node = 'PROFILER'
instrument = 'DOSTA'
method = 'recovered_cspp'
start_date = '2019-01-01'
stop_date = '2019-12-31'
drop_paired = TRUE
simplify_data = TRUE
####---------------------------------###


info = ooi_site_info(site,user,token)
max_depth = info[["Surface Piercing Profiler: SP001"]][["Dissolved Oxygen: 01-DOSTAJ000"]][["maxdepth"]]
min_depth = info[["Surface Piercing Profiler: SP001"]][["Dissolved Oxygen: 01-DOSTAJ000"]][["mindepth"]]


url = ooi_create_url(site = site,node = node,instrument = instrument,method = method,start_date = start_date,stop_date = stop_date) #Generate a URL.
response = ooi_submit_request(url,user = user, token = token)  #Issue the request.
remote = ooi_get_location(response,drop_paired = drop_paired)  #Check data status and return remote locations of data on the Thredds server.
local = ooi_download_data(remote)  #Uncomment this option if you want to download the data.
lol = ooi_get_data(local,simplify_data = simplify_data)  #Read in the remote NetCDFs.
data = data.frame(lol[['data']]) #Organize the data so it is easier on the eyes. Offers a table of the data if data is 1D. Data stays as a list of lists if 2D.
varunits = data.frame(lol[['variables_units']])  #Offers a table of the variables and associated units.

data$time = as.POSIXct(data$profiler_timestamp,tz='UTC',origin='1970-01-01')  #Convert time to something understandable. profiler_timestamp is a common value among CSPP sensors.
data = data[(data$pressure_depth<max_depth & data$pressure_depth>min_depth),]  #Only keep values that are between the max and min depth determined through the ooi_site_depth() function.


#-------------------------------------#



myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 300))


oxy_plot = ggplot(data,aes(time,pressure_depth,colour = dissolved_oxygen)) +
  geom_point() +
  scale_y_reverse() +
  sc +
  labs(x = 'Time',y = 'Pressure (dbars)',colour = 'Dissolved Oxygen (umol/kg)', title = 'CE02SHSP Dissolved Oxygen',subtitle = 'All Deployments For 2019')
oxy_plot
