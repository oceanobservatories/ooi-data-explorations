#Comment out the next 5 lines if you already have these packages installed.
#using Pkg; #Using the Julia Pkg package....
#packages = ["HTTP","JSON","DataFrames","NCDatasets"] #Add packages to be downloaded here.
#for package in packages
#    Pkg.add(package)  #Download package.
#end

using HTTP, JSON, DataFrames, NCDatasets  #Bring the packages in so that we can use them.

#------CHANGE ME-----#
user = "OOI-API-USER-HERE"
token = "OOI-API-TOKEN-HERE"
directory = "C:/Users/Ian/Desktop/test"
#--------------------#

site = "CE01ISSP"
node = "SP001"
instrument = "09-CTDPFJ000"
method = "recovered_cspp"
stream = "ctdpf_j_cspp_instrument_recovered"
start_date = "2018-06-01"
stop_date = "2019-06-02"

url = ooi_create_url(site,node,instrument,method,stream,start_date = start_date, stop_date = stop_date)
response = ooi_submit_request(url,user,token)
dodsC = ooi_get_location(response)
fileServer = ooi_download_data(dodsC,directory)
