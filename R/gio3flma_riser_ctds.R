#Desciption
#This script will request data from all of the riser CTD on the Global Irminger Flanking A Mooring between June 01, 2018 and June 02,2018.
#It does not consider efficiency when issuing data requests to OOINet.
#It will make a request and then check to see if the data has been made available before moving on to request the next dataset.

require(devtools)
install_github("oceanobservatories/ooim2mr")

require(jsonlite)
require(rio)
require(stringr)
require(httr)
require(ncdf4)
require(ooim2mr)

user = "YOUR-OOI-API-USER-HERE"   #OOI API Username
token = "YOUR-OOI-API-TOKEN-HERE"  #OOI API Token


#Request mooring information. The following function returns a list of lists that contains the node, instrument, and stream metatdata associated with the requested mooring.
mooring = ooi_site_info('GI03FLMA', user, token )

#We can search through the "mooring" object for information about each instrument, such as the reference designator and the depth it is deployed at.
#If you are using RStudio, you can type View(mooring) in the console to visualize the list of lists.
refdes = c()
inst_depth = c()
for (node in names(mooring)){  #For each node on the mooring...
  for (inst in names(mooring[[node]])){  #For each instrument on each node...
    rd = mooring[[node]][[inst]][['refdes']]  #Identify the reference designator for that combination.
    depth = mooring[[node]][[inst]][['maxdepth']]  #Find the depth of that instrument.
    refdes = c(refdes,rd)  #Concatenate previous loop with this loop to get a list of reference designators.
    inst_depth = c(inst_depth,depth)  #Concatenate previous loop with this loop to get a list of instrument depths.
  }
}
inst_df = cbind(refdes,inst_depth)  #Combine the list of reference designators and instrument depths.
inst_df = inst_df[order(inst_depth),]  #Order each instrumetn by depth so our figure is in order later.


#We only want CTD data, so we can drop any reference designators that don't contain "CTD".
ctd_df <- inst_df[grep('CTD', inst_df[, "refdes"]),]  #Drop instruments taht aren't CTDs.
banana = str_split(ctd_df[,'refdes'], "-")  #Split the reference designator so we can pull out the ctd id.
ctd_id = c()
for (id in banana){  #For each split reference designator...
  ctd = id[4]  #The ctd id is in the 4th position.
  ctd_id = c(ctd_id,ctd)  #Concatenate previous loops to get a list of ctd ids.
}

#Next, we can generate request URLs for Flanking Mooring A using a for loop that loops through each CTD id. Data between 06/01/2018 and 06/02/2018.
urls = c()
for (ctd in ctd_id){  #For each ctd in our list of ctds...
  site = 'GI03FLMA'
  node = 'MIDWATER'
  instrument = ctd
  method = 'recovered_inst'
  start_date = '2018-06-01'
  stop_date = '2018-06-02'
  url = ooi_create_url(site,node,instrument,method,start_date=start_date,stop_date=stop_date) #Make a url for requesting data using the above information.
  urls = c(urls,url)  #Concatenate each loop so that we have a list of request urls.
}

#This next part issues a request for each url in the list of urls and then checks the data status before issuing the next request.
#To speed things up, you could issue all of your data requests before you check the data status.
remotes = list()
for (url in urls){  #For each url in the list of urls.
  response = ooi_submit_request(url,user,token)  #Submit the request.
  remotes[url] = ooi_get_location(response,drop_paired=TRUE)  #Get the remote location of the data. Add it to the list of remote locations.
}

#Next, we'll download the fileServer NetCDFs that are paired with the remote dodsC files.
files = c()
for (remote in remotes){  #For each remote dataset...
  for (link in remote){   #For each link in each remote dataset...
    file = ooi_download_data(remote)  #Download the data.
    files = c(files,file)  #Concatenate filenames.
  }
}
files = data.frame(files)  #Put the filenames into a data.frame.

#Then we will bring the data into the workspace by the ctd id.
data = list()
for (id in ctd_id){  #For each CTD id.
  ctd <- toString(files[grep(id, files[, "files"]),])  #Identify which files have data for that CTD.
  hold =  ooi_get_data(ctd,simplify_data = TRUE)
  data[id] = list(hold)  #Bring that data into the workspace so that we get a list of lists of data.
}


#Finally, we will quickly create a plot in a new window, where the number of subplots is equal to the number of CTDs.
#dev.new(width = 6, height = 9, unit = "in")  #Create a new window for the plot
par(mfrow = c(length(ctd_id),1),mar = c(1,1,1,1))  #Create a number of subplots that is equal to the number of CTD ids.
for (i in 1:length(data)){  #For each list in data.
  d = data[[i]][['data']]  #Move further down the list to make it easier to access the data.
  time = as.POSIXct(d$time,tz = 'GMT',origin = '1900-01-01')  #Get the time.
  sal = d$practical_salinity  #Get the salinity.
  plot(time,sal,main = sprintf('%s: %sm',ctd_id[i],ctd_df[,'inst_depth'][i]),ylim = c(32,35),ylab = 'Salinity')   #Plot the data.
}
