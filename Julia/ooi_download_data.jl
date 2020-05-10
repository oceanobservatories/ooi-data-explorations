# Download OOI data to a given directory.
#
# @return A list of filepaths where data is located.
#   Files are renamed based on site, node, instrument, method, etc.
#
# @param dodsC A list of dodsC OpenDAP urls.
# @param directory An operating specific directory to save data to.



function ooi_download_data(dodsC,directory)
    file_df = DataFrame(num = [],url = [])  #Create a holder dataframe.
    for i = 1:length(dodsC)
        drop_mismatch = replace(dodsC[i],"#fillmismatch" => "")
        fileServer = replace(drop_mismatch,"dodsC"=>"fileServer")
        banana = split(fileServer,'_')
        window_nc = banana[end]
        banana = split(banana[6],'/')
        deployment = banana[end]
        banana = split(fileServer,'-')
        site = banana[2]
        node = banana[3]
        instrument = banana[5]
        method  = banana[6]
        banana = split(banana[7],"/")
        stream = banana[1]
        filename = string(site,'_',node,'_',instrument,'_',method,'_',stream,'_',deployment,'_'window_nc)
        filepath = string(directory,'/',filename)
        download(fileServer,filepath)
        push!(file_df,(i,filepath))
    end #End of for loop.
    file_df = Vector{String}(file_df.url)
    println("Data are now downloaded.")
    return file_df
end #End of function.
