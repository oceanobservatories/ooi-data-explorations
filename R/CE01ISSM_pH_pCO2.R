# Description
# An example script that requests pH, pCO2, and CTD data from the 7m package located at the OOI Oregon Inshore site (CE01ISSM) for summer 2019.
# Data is then run through the seacarb package.
# Seacarb Manual: https://cran.r-project.org/web/packages/seacarb/seacarb.pdf

# List of Coastal OOI sites that have packages with co-located pH, pCO2, and CTD sensors.
# site = 'CE01ISSM', node = 'MIDWATER' (Oregon Inshore Surface Mooring - 7m)
# site = 'CE01ISSM', node = 'SEAFLOOR' (Oregon Inshore Surface Mooring - 25m)
# site = 'CE02SHSM', node = 'MIDWATER' (Oregon Shelf Surface Mooring - 7m)
# site = 'CE02SHBP', node = 'SEAFLOOR' (Oregon Shelf Benthic Package - 80m)
# site = 'CE04OSSM', node = 'MIDWATER' (Oregon Offshore Surface Mooring - 7m)
# site = 'CE04OSBP', node = 'SEAFLOOR' (Oregon Offshore Benthic Package - 580m)
# site = 'CE04OSPS', node = 'MIDWATER' (Oregon Offshore 200m Platform - 200m)
# site = 'CE04OSPS', node = 'PROFILER' (Oregon Offshore Profiler - 20 to 200m)
# site = 'CE06ISSM', node = 'MIDWATER' (Washington Inshore Surface Mooring - 7m)
# site = 'CE06ISSM', node = 'SEAFLOOR' (Washington Inshore Surface Mooring - 29m)
# site = 'CE07SHSM', node = 'MIDWATER' (Washington Shelf Surface Mooring - 7m)
# site = 'CE07SHSM', node = 'SEAFLOOR' (Washington Shelf Surface Mooring - 87m)
# site = 'CE09OSSM', node = 'MIDWATER' (Washington Offshore Surface Mooring - 7m)
# site = 'CE09OSSM', node = 'SEAFLOOR' (Washington Offshore Surface Mooring - 600m)
# site = 'CP01CNSM', node = 'MIDWATER' (Pioneer Central Surface Mooring - 7m)
# site = 'CP01CNSM', node = 'SEAFLOOR' (Pioneer Central Surface Mooring - 133m)
# site = 'CP03ISSM', node = 'MIDWATER' (Pioneer Inshore Surface Mooring - 7m)
# site = 'CP03ISSM', node = 'SEAFLOOR' (Pioneer Inshore Surface Mooring - 92m)
# site = 'CP04OSSM', node = 'MIDWATER' (Pioneer Offshore Surface Mooring - 7m)
# site = 'CP04OSSM', node = 'SEAFLOOR' (Pioneer Offshore Surface Mooring - 450m)
# site = 'RS01SBPS', node = 'MIDWATER' (Oregon Slope 200m Platform - 200m)
# site = 'RS01SBPS', node = 'PROFILER' (Oregon Slope Profiler - 5 to 200m)
# site = 'RS03AXPS', node = 'MIDWATER' (Axial Base 200m Platform - 200m)
# site = 'RS03AXPS', node = 'PROFILER' (Axial Base Profiler - 5 to 200m)


# Install required packages. You can comment out this section if you've already installed them.
try({remove.packages("ooirtest")})  #Remove the old OOI R repo. May throw an error.
pkgs <- c('devtools','httr','ncdf4','lubridate','seacarb') #Packages to install.
for (pkg in pkgs){
  install.packages(pkg)
}


#Install the ooim2mr package.
require(devtools)  #Load the devtools package so that we can install the OOI R repo from GitHub.
install_github("IanTBlack/ooim2mr")  #Install the new build of the OOI R repo.
detach("package:devtools",unload=TRUE)  #Unload the devtools package because it might interfere with the seacarb package.


# Load the necessary packages.
require(ncdf4)  #Used for accessing OpenDAP NetCDFs.
require(httr)  #Used for issuing requests.
require(ooim2mr)  #Used for getting data from OOINet.
require(jsonlite)  #Used for parsing request responses.
require(stringr)  #Used for navigating the ooim2mr lookup table.
require(data.table)  #Makes it easier to merge data by nearest time.
require(lubridate)  #For rounding dates.
require(seacarb)  #For people good at chemistry.


#Set the directory to save data to.
path = "C:/Users/Ian/Desktop/test"  #Change this to match a format for your operating system.
setwd(path) #Set the working directory. Data folders will be added to this path further down the script.


#Define OOINet credentials.
user = 'OOIAPI-BCJPAYP2KUVXFX'  #OOI.CSPP@gmail.com username
token = 'D3HV2X0XH1O'  #OOI.CSPP@gmail.com token


#Set default information for requests.
site = "CE01ISSM"  #The site we want data from. https://oceanobservatories.org/site/ce01issm/
node = "MIDWATER"  #The node or package we want data from.
method = "recovered_inst"  #How we want to get data. For pH and pCO2, this will either be recovered_inst, recovered_host, telemetered, or streamed.
start_date = '2019-05-01'  #Start date of data request.
stop_date = '2019-09-30'  #Stop date of data request.


#Create request urls.
pH_url = ooi_create_url(site,node,'pH',method,start_date = start_date,stop_date = stop_date) #This URL will request pH data from the CE01ISSM NSIF between May 01 and Sept 30, 2019. 
pCO2_url = ooi_create_url(site,node,'pCO2',method,start_date = start_date,stop_date = stop_date)
CTD_url = ooi_create_url(site,node,'CTD',method,start_date = start_date,stop_date = stop_date)


#Submit requests.
pH_r = ooi_submit_request(pH_url,user,token)
pCO2_r = ooi_submit_request(pCO2_url,user,token)
CTD_r = ooi_submit_request(CTD_url,user,token)


#Get the OpenDAP urls.
pH_opendap = ooi_get_location(pH_r,drop_paired = TRUE)  #Don't use paired CTD data because it is filled with NaNs. Download the CTD data separately.
pCO2_opendap = ooi_get_location(pCO2_r,drop_paired = TRUE) 
CTD_opendap = ooi_get_location(CTD_r,drop_paired = TRUE)


#Download data. If running on MacOS or Linux, you can bypass this section and input OpenDAP urls directly into the ooi_get_data function.
pH_path = sprintf('%s%s',path,'/pH') 
dir.create(pH_path) #Create a subfolder for pH data.
pH_files = ooi_download_data(pH_opendap,directory = pH_path) #Download the data.

pCO2_path = sprintf('%s%s',path,'/pCO2')
dir.create(pCO2_path)
pCO2_files = ooi_download_data(pCO2_opendap,directory = pCO2_path)

CTD_path = sprintf('%s%s',path,'/CTD')
dir.create(CTD_path)
CTD_files = ooi_download_data(CTD_opendap,directory = CTD_path)


#Bring data into the workspace.
pH_lol = ooi_get_data(pH_files,simplify_data= TRUE)  #Merge data from multiple NetCDFs and drop data products that are confusing or generally useless.
pH_data = data.frame(pH_lol[['data']])  #The first list of the ooi_get_data return is always the data.
pH_vars = data.frame(pH_lol[['variables_units']]) #The second list of the ooi_get_data return is always a list of variables and units.

pCO2_lol = ooi_get_data(pCO2_files,simplify_data = TRUE)
pCO2_data = data.frame(pCO2_lol[['data']])
pCO2_vars = data.frame(pCO2_lol[['variables_units']])

CTD_lol = ooi_get_data(CTD_files,simplify_data = TRUE)
CTD_data = data.frame(CTD_lol[['data']])
CTD_vars = data.frame(CTD_lol[['variables_units']])


#Mooring CTDs sample ~ every 15 minutes. pH and pCO2 sensors sample ~ every 2 hours.
#Round dates to the nearest two hours (0000,0200,0400, etc.) and average grouped datetimes for each dataset.
pH_data$dt = round_date(as.POSIXct(pH_data$time, tz = 'UTC', origin = '1900-01-01'),"2 hours")  #Round datetimes to the nearest 2 hour mark.
pH_data= aggregate(pH_data,by=list(hour = pH_data$dt),FUN=mean)  #Average rows that have the same datetime.

pCO2_data$dt = round_date(as.POSIXct(pCO2_data$time, tz = 'UTC', origin = '1900-01-01'),"2 hours")
pCO2_data= aggregate(pCO2_data,by=list(hour = pCO2_data$dt),FUN=mean)

CTD_data$dt = round_date(as.POSIXct(CTD_data$time, tz = 'UTC', origin = '1900-01-01'),"2 hours")
CTD_data = aggregate(CTD_data,by=list(hour = CTD_data$dt),FUN=mean)


#Merge data by nearest neighbor. 
CTD_table = data.table(CTD_data)  #Convert dataframes to datatables so they can be merged easily.
pH_table = data.table(pH_data)
pCO2_table = data.table(pCO2_data)
CTD_pH_merge = CTD_table[pH_table, on = "hour"]  #Merge based on the "hour" column.
CTD_pCO2_merge = CTD_table[pCO2_table,on = "hour"]
data = CTD_pH_merge[CTD_pCO2_merge, on = "hour"]


#Only keep columns we care about so that we can run them through the carb function.
#Order: datetime, hydrostatic pressure (dbars), temperature (degC), salinity (PSU), pH, pCO2(uatm)
data = data[,c('hour','ctdbp_seawater_pressure','ctdbp_seawater_temperature','practical_salinity','phsen_abcdef_ph_seawater','pco2_seawater'),with = FALSE]


#Drop rows that contain nulls and missing values. Missing values can cause the carb function to fail.
data <- na.omit(data) 
data <- data[complete.cases(data),]

#Run data through the carb function.
#Flag 21 for carb function is for PCO2 and pH respectively. See pg. 24 of seacarb manual.
sw_carb_sys = carb(flag = 21, 
                   var1 = data$pco2_seawater,
                   var2 = data$phsen_abcdef_ph_seawater,  #phsen_abcdef_ph_seawater is named this way because the API was written by computer scientists and not oceanographers.
                   S = data$practical_salinity,
                   T = data$ctdbp_seawater_temperature, 
                   P = data$ctdbp_seawater_pressure/10) #Pressure must be in bars. OOI pressure is generally in decibars.

datetime = data$hour  #Pull the datetime from the data.
sw_carb_sys = cbind(datetime,sw_carb_sys) #Combine time with each row of the sw_carb_sys dataframe.

#Look at the data.
View(sw_carb_sys)