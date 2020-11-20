#Description: An example script that will download CE07SHSM VELPT (Aquadopp) recovered instrument and host data for the year 2019.
#Tested in R 3.6.2 on a Windows 10 machine.
#Author: Ian Black (blackia@oregonstate.edu)
#Date: 2020-09-26

install.packages("ggplot2")

install.packages("devtools")  #Install the devtools library from CRAN.
library(devtools) #The devtools package allows the use of non-CRAN packages (e.g. GitHub).
install_github("oceanobservatories/ooim2mr") #Install the ooim2mr package from GitHub.
detach("package:devtools",unload=TRUE)

require(ooim2mr)  #This allows use of custom OOINet functions.

#These packages are required for the ooim2mr package.
library(jsonlite) #Allows for parsing of responses from OOINet.
library(httr)  #Allows for issuing requests to OOINet.
library(stringr)  #Allows for string manipulation.
library(ncdf4) #Allows for manipulation of NetCDFs.

#This package lets you plot things in a more attractive format compared to base R.
library(ggplot2)


save_directory = 'C:/Users/Ian/Desktop' #Where data will be saved.
user = 'OOIAPI-BCJPAYP2KUVXFX' #Your own OOINet credentials should go here.
token = 'ARZNUYHENKG'

site = 'CE07SHSM'  #WA Shelf mooring.
node = 'SURFACE'  #Buoy
instrument = 'VELPT'  #Aquadopp
start_date = '2019-01-01' #Format in YYYY-MM-DD.
stop_date = '2019-12-31'

url_host = ooi_create_url(site = site, #Generate URL for data that is recorded by the mooring controller.
                     node = node, 
                     instrument = instrument, 
                     method = 'recovered_host',
                     start_date = start_date,
                     stop_date = stop_date)

url_inst = ooi_create_url(site = site, #Generage URL for data that is offloaded from the instrument after recovery.
                          node = node, 
                          instrument = instrument, 
                          method = 'recovered_inst',
                          start_date = start_date,
                          stop_date = stop_date)

response_host = ooi_submit_request(url = url_host,user = user,token = token)  #Issue data requests.
response_inst = ooi_submit_request(url = url_inst,user = user,token = token)

remote_host = ooi_get_location(response = response_host,drop_paired = TRUE) #Get the remote locations of the data.
remote_inst = ooi_get_location(response = response_inst,drop_paired = TRUE)

local_host = ooi_download_data(remote_host,save_directory) #Download data to a specified directory.
local_inst = ooi_download_data(remote_inst,save_directory)

dv_host = ooi_get_data(local_host,simplify_data = FALSE) #Download the data. Return is a list of lists. The first list is the data. The second list are the variables in the dataset.
dv_inst = ooi_get_data(local_inst,simplify_data = FALSE)


data_host = data.frame(dv_host[['data']])  #Separate data from variables.
host_vars = data.frame(dv_host[['variables_units']])

data_inst = data.frame(dv_inst[['data']])
inst_vars = data.frame(dv_inst[['variables_units']])
print(inst_vars)

host_pres <- data_host$pressure_mbar  #Identify the data we care about and put it in its own dataframe.
host_time <- as.POSIXlt(data_host$time,origin="1900-01-01",tz = "GMT")  #The time variable is seconds since 1900-01-01.
host <- data.frame(host_time,host_pres) 

inst_pres <- data_inst$pressure_mbar
inst_time <- as.POSIXlt(data_inst$time,origin="1900-01-01",tz = "GMT")
inst <- data.frame(inst_time,inst_pres)

#Show the length of the datasets.
print(length(host_time))  
print(length(inst_time))
#It would appear the recovered instrument dataset has one more data point than the recovered host.
#Not bad. It probably doesn't matter what dataset you use.

#Plot of host data.
h <- ggplot(host,aes(host_time,host_pres)) +
  geom_line() +
  xlab("time") + 
  ylab("millibars")+
  ggtitle("CE07SHSM VELPT BUOY Host")



#Plot of instrument data.
i <- ggplot(inst,aes(inst_time,inst_pres)) +
  geom_line() +
  xlab("time") + 
  ylab("millibars")+
  ggtitle("CE07SHSM VELPT BUOY Inst")


h #Show plots
i