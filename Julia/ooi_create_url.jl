# Create an OOI data request URL from input parameters.
#
# @return A URL string that can be input into ooi_submit_request.
#
# @param site An eight character OOI site.
# @param node A five character OOI node.
# @param instrument A twelve character OOI instrument.
#   The third character must be a dash (-)
# @param method Method for data delivery.
# @param stream  Stream for data.
# @param start_date UTC in format YYYY-mm-dd.
# @param start_time UTC in format HH:MM:SS.
# @param stop_date UTC in format YYYY-mm-dd.
# @param stop_time UTC in format HH:MM:SS.



function ooi_create_url(site,node,instrument,method,stream;
                        start_date ="2010-01-01",start_time = "00:00:00",
                        stop_date= "2040-12-31",stop_time="23:59:59")

    url = string("https://", user,":",token,
                 "@ooinet.oceanobservatories.org/api/m2m/12576/sensor/inv/",
                 site,"/",node,"/",instrument,"/",method,"/",stream,
                 "?beginDT=",start_date,"T",start_time,
                 ".000Z&endDT=",stop_date,"T",stop_time,".999Z")

    return url
end  #End of function.
