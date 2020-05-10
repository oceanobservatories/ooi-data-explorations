# Submit a data request to OOINet.
#
# @return A queryable response object if the request is successful.
#
# @param url A request URL string generated from ooi_create_url or by hand.
# @param user Your OOINet username.
# @param token Your OOINet token.



function ooi_submit_request(url,user,token)
    response = HTTP.get(url)  #Get the response from OOINet.
    if response.status == 200 #If the response says successful (200)
        println("Request successful.")
        return response #Return the response object.
    elseif response.status == 401
        println("Unauthorized request. Recommend checking user and token.")
    elseif response.status == 404
        println("No data available for the requested window.")
    else
        println(string("Unanticipated error code (",response.status,")."))
    end  #End of elif ladder.
end #End of function.
