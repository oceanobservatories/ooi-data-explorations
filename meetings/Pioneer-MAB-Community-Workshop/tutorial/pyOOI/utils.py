import datetime
import pandas as pd


def ntp_seconds_to_datetime(ntp_seconds):
    """Convert OOINet timestamps to unix-convertable timestamps."""
    # Specify some constant needed for timestamp conversions
    ntp_epoch = datetime.datetime(1900, 1, 1)
    unix_epoch = datetime.datetime(1970, 1, 1)
    ntp_delta = (unix_epoch - ntp_epoch).total_seconds()

    return datetime.datetime.utcfromtimestamp(ntp_seconds - ntp_delta)


def convert_time(ms):
    """Calculate UTC timestamp from OOI milliseconds"""
    if ms is None:
        return None
    else:
        return datetime.datetime.utcfromtimestamp(ms/1000)


def unix_epoch_time(date_time):
    """Convert a datetime to unix epoch microseconds."""
    # Convert the date time to a string
    date_time = int(pd.to_datetime(date_time).strftime("%s"))*1000
    return date_time