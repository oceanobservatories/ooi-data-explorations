import re
import pandas as pd

class QualityFlags():
    """QARTOD QC-flag definitions"""
    GOOD = 1
    UNKNOWN = 2
    SUSPECT = 3
    BAD = 4
    MISSING = 9

FLAGS = QualityFlags

def check_fill(flag):
    """Check if an OOI discrete sample flag is a fill value (-9999999)"""
    if pd.isna(flag):
        return True
    elif str(flag) == "-9999999":
        return True
    elif "1" not in str(flag):
        return True
    else:
        return False


def parse_flag(flag):
    """Parse an OOI discrete sample flag. Returns fill or nan when appropriate."""
    locs=[]
    for match in re.finditer("1", flag[::-1], re.S):
        locs.append(match.span()[0])
    return locs


def interp_ctd_flag(flag):
    """Interpret OOI discrete sample CTD flags to QARTOD QC-flags."""
    # First filter for fill
    if check_fill(flag):
        return flag
    else:
        parsed_flag = parse_flag(flag)
        max_bit = max(parsed_flag)
        if max_bit == 1:
            return QualityFlags.MISSING
        elif max_bit == 2:
            return QualityFlags.GOOD
        elif max_bit == 3:
            return QualityFlags.SUSPECT
        elif max_bit == 4:
            return QualityFlags.BAD
        else:
            return QualityFlags.UNKNOWN


def interp_discrete_flag(flag):
    """Interpret OOI discrete sample Bottle flags to QARTOD QC-flags."""
    # First filter for fill values
    if check_fill(flag):
        return flag
    else:
        parsed_flag = parse_flag(flag)
        max_bit = max(parsed_flag)
        if max_bit == 1:
            return QualityFlags.MISSING
        elif max_bit == 2:
            return QualityFlags.GOOD
        elif max_bit == 3:
            return QualityFlags.SUSPECT
        elif max_bit == 4:
            return QualityFlags.BAD
        else:
            return QualityFlags.UNKNOWN

def interp_replicate_flag(flag):
    """
    Interpret OOI discrete sample Bottle replicate flags. Returns a boolean
    if a sample is a duplicate/replicate sample.
    """
    # First filter for fill values
    if check_fill(flag):
        return flag
    else:
        parsed_flag = parse_flag(flag)
        max_bit = max(parsed_flag)
        if max_bit == 3 or max_bit == 4:
            return True
        else:
            return False


def interp_niskin_flag(flag):
    """Interpret OOI discrete Niskin Bottle flags to QARTOD QC-flags."""
    if check_fill(flag):
        return flag
    else:
        parsed_flag = parse_flag(flag)
        max_bit = max(parsed_flag)
        if max_bit == 1:
            return QualityFlags.MISSING
        elif max_bit == 2:
            return QualityFlags.GOOD
        elif max_bit == 3 or max_bit == 4 or max_bit == 5:
            return QualityFlags.SUSPECT
        else:
            return QualityFlags.UNKNOWN


def convert_times(x):
    """Parse times in OOI discrete summary spreadsheets to python datetimes."""
    if type(x) is str:
        x = x.replace(" ","")
        x = pd.to_datetime(x, utc=False)
    else:
        pass
    return x


def not_statistically_sigificant(x):
    """
    Replace OOI discrete nutrient sample values that are not statistically
    significant with zero.
    """
    if type(x) is str:
        if "<" in x:
            x = 0
    return x


def findNearest(bottleData, buoyLoc, maxDist):
    """Find the bottle sample values within a maximum distance from the buoy
    
    Parameters
    ----------
    bottleData: (pd.DataFrame -> strings or floats)
        A tuple of (latitude, longitude) values in decimal degrees of the bottle sample location
    buoyLoc: (tuple -> floats)
        A tuple of (latitude, longitude) values in decimal degrees of the buoy location
    maxDist: (float)
        Maximum distance in km away for a sample location from the buoy location
    
    Returns
    -------
    mask: (boolean)
        Returns True or False boolean if sampleLoc < maxDist from buoyLoc
    """
    # Get the startLat/startLon as floats
    startLat = bottleData["Start Latitude [degrees]"].apply(lambda x: float(x))
    startLon = bottleData["Start Longitude [degrees]"].apply(lambda x: float(x))
    
    # Calculate the distance
    distance = []
    for lat, lon in zip(startLat, startLon):
        sampleLoc = (lat, lon)
        distance.append(hs.haversine(sampleLoc, buoyLoc))
    
    # Filter the results
    return [d <= maxDist for d in distance]


def findSamples(bottleData, buoyLoc, buoyDepth, maxDist, depthTol):
    """Find the bottle sample values within a maximum distance from the buoy
    
    Parameters
    ----------
    bottleData: (pd.DataFrame -> strings or floats)
        A tuple of (latitude, longitude) values in decimal degrees of the bottle sample location
    buoyLoc: (tuple -> floats)
        A tuple of (latitude, longitude) values in decimal degrees of the buoy location
    buoyDepth: (float)
        Deployment depth of the instrument
    maxDist: (float)
        Maximum distance in km away for a sample location from the buoy location
    depthTol: (float)
        Maximum depth difference to select samples from the buoyDepth
    
    Returns
    -------
    mask: (boolean)
        Returns True or False boolean if sampleLoc < maxDist from buoyLoc
    """
    # Filter for the nearest samples
    nearest = findNearest(bottleData, buoyLoc, maxDist)
    bottleData = bottleData[nearest]
    
    # Filter based on depth
    depthMin = buoyDepth - depthTol
    depthMax = buoyDepth + depthTol
    bottleData = bottleData[(bottleData["CTD Depth [m]"] >= depthMin) & (bottleData["CTD Depth [m]"] <= depthMax)]
    
    return bottleData