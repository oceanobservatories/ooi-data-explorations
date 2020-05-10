function ooi_get_location(response)
    response2json =JSON.parse(String(response.body))  #Convert the response to a JSON object.
    check = string(response2json["allURLs"][2],"/status.json")  #Build the checker URL.
    for i = 1:1800  #For roughly 30 minutes...
        try  #Try to get the status from the check URL...
            status = HTTP.get(check)  #Attempt to make a request to status.json and get the response.
            status = String(status.body)
            if occursin("complete",status) == true  #If it says complete...
                println("Request complete.")
                break  #Exit the for loop once complete.
            else  #Habit. No real reason to have an else here.
            end
        catch  #If status returns a 404 error....
            println("Checking data status...")
            sleep(1) #Sleep for 1 second before checking again.
        end #End of try/catch.
    end #End of for loop.
    data_url = response2json["allURLs"][1]  #Data catalog is located here.
    data_response = HTTP.get(data_url,cookies = true)  #Make a request for html text of the catalog.
    data = String(data_response.body)
    nc_pattern = r"dataset=ooi.+?.nc'>"  #Use this string to parse out NetCDF urls.
    url_parts = findall(nc_pattern,data)  #Identify the indexes of the NetCDF urls.
    url_df = DataFrame(num = [],url = [])  #Create a holder dataframe.
    for i = 1:length(url_parts)  #For each index range...
        url_part = data[url_parts[i]]  #Pull the url out based on the index range.
        nc_url = replace(url_part,"dataset=" => "https://opendap.oceanobservatories.org/thredds/dodsC/")  #Add a dodsC url by replacing dataset=.
        nc_url = replace(nc_url,"'>"=>"#fillmismatch") #Apply a fillmismatch flag by replacing the trailing '>
        nc_url = string(nc_url)  #Make sure it is a string.
        push!(url_df,(i,nc_url))  #Append a new row to the dataframe.
    end #End of for loop.
    location = Vector{String}(url_df.url)  #Create an array of NetCDF urls from the dataframe url column.
    println(location)  #Print the remote locations so we know that they are being formatted properly.
    return location
end #End of function.
