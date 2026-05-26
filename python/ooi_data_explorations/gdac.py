#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Provides access to the IOOS GliderDAC (NGDAC) via ERDDAP for use in
    creating QARTOD test values for OOI glider instrument classes across the
    Coastal Endurance, Coastal Pioneer, and Global arrays.
"""
import numpy as np
import pandas as pd
import requests
import sys
import xarray as xr

from concurrent.futures import ThreadPoolExecutor
from erddapy import ERDDAP
from functools import partial
from tqdm import tqdm

# ERDDAP variable names always included regardless of user request
_REQUIRED_ERDDAP_VARS: list[str] = [
    'precise_time', 'precise_lon', 'precise_lat', 'depth', 'pressure'
]

# Maps internal (post-rename) variable names to their ERDDAP names
_SCIENCE_VAR_MAP: dict[str, str] = {
    'temperature': 'temperature',
    'conductivity': 'conductivity',
    'salinity': 'salinity',
    'density': 'density',
    'bback': 'backscatter',
    'fluorometric_cdom': 'CDOM',
    'estimated_chlorophyll': 'chlorophyll',
    'dissolved_oxygen': 'dissolved_oxygen',
    'par': 'PAR',
}

# Rename map applied to downloaded datasets (ERDDAP name -> internal name)
_RENAME_MAP: dict[str, str] = {
    'precise_time': 'time',
    'precise_lon': 'longitude',
    'precise_lat': 'latitude',
    'backscatter': 'bback',
    'CDOM': 'fluorometric_cdom',
    'chlorophyll': 'estimated_chlorophyll',
    'PAR': 'par',
}

# Public list of valid science variable names (internal/post-rename names)
GLIDER_SCIENCE_VARIABLES: list[str] = list(_SCIENCE_VAR_MAP.keys())

# Maps OOI site prefix to GliderDAC institution name.
# TODO: Confirm the correct institution strings for CP and Global arrays
#     by querying the NGDAC ERDDAP search endpoint.
GLIDER_INSTITUTIONS: dict[str, str] = {
    'CE': 'ooi_coastal_endurance',
    'CP': 'ooi_coastal_pioneer',       # TODO: confirm
    'GA': 'ooi_global_argentine_basin', # TODO: confirm
    'GI': 'ooi_global_irminger_sea',    # TODO: confirm
    'GP': 'ooi_global_station_papa',    # TODO: confirm
    'GS': 'ooi_global_southern_ocean',  # TODO: confirm
}


def glider_variables() -> list[str]:
    """
    Return the science variables available from the GliderDAC.

    These are the valid values for the ``variables`` parameter of
    :func:`collect_glider`. Navigation and coordinate variables
    (time, latitude, longitude, depth, pressure) are always included
    and are not listed here.

    :return: list of valid internal variable names
    """
    return GLIDER_SCIENCE_VARIABLES.copy()


def glider_institution(site: str) -> str:
    """
    Return the GliderDAC institution name for an OOI site.

    Uses the two-character site prefix to look up the corresponding
    institution string used in GliderDAC ERDDAP searches.

    :param site: OOI site designator (e.g. 'CE05MOAS', 'CP05MOAS')
    :return: GliderDAC institution name
    :raises ValueError: If the site prefix has no known institution mapping
    """
    prefix = site[:2].upper()
    if prefix not in GLIDER_INSTITUTIONS:
        raise ValueError(
            f"No GliderDAC institution known for site prefix '{prefix}'. "
            f"Known prefixes: {sorted(GLIDER_INSTITUTIONS.keys())}. "
            f"Update GLIDER_INSTITUTIONS in gdac.py to add it."
        )
    return GLIDER_INSTITUTIONS[prefix]


def _create_box(lat: float, lon: float, extent: float) -> list[float]:
    """
    Create a lat/lon bounding box around a point.

    :param lat: Center latitude in decimal degrees
    :param lon: Center longitude in decimal degrees
    :param extent: Half-width of the box in nautical miles
    :return: [min_lat, max_lat, min_lon, max_lon]
    """
    return [
        lat - (extent / 60.),
        lat + (extent / 60.),
        lon - (extent * np.cos(np.radians(lat)) / 60.),
        lon + (extent * np.cos(np.radians(lat)) / 60.)
    ]


def _list_gliders(bounding_box: list[float], institution: str) -> np.ndarray:
    """
    List all delayed-mode glider datasets for an institution in the
    GliderDAC that fall within the bounding box.

    :param bounding_box: [min_lat, max_lat, min_lon, max_lon]
    :param institution: GliderDAC institution name (see
        :data:`GLIDER_INSTITUTIONS` or :func:`glider_institution`)
    :return: array of GliderDAC dataset IDs
    """
    server = ERDDAP(server='ngdac')
    search_url = server.get_search_url(
        response='csv',
        search_for='delayed',
        institution=institution,
        min_lat=bounding_box[0],
        max_lat=bounding_box[1],
        min_lon=bounding_box[2],
        max_lon=bounding_box[3],
    )
    search = pd.read_csv(search_url)
    return search['Dataset ID'].values


def _download_glider(
    dataset_id: str,
    bounding_box: list[float],
    erddap_vars: list[str]
) -> xr.Dataset | None:
    """
    Download a single glider dataset from the GliderDAC.

    Creates a fresh ERDDAP connection per call so that concurrent
    downloads do not share mutable server state.

    :param dataset_id: GliderDAC dataset ID
    :param bounding_box: [min_lat, max_lat, min_lon, max_lon]
    :param erddap_vars: ERDDAP variable names to request
    :return: xarray Dataset or None if no data found in the bounding box
    """
    server = ERDDAP(server='ngdac')
    server.constraints = {
        'latitude>=': bounding_box[0],
        'latitude<=': bounding_box[1],
        'longitude>=': bounding_box[2],
        'longitude<=': bounding_box[3],
    }
    server.protocol = 'tabledap'
    server.variables = erddap_vars
    server.dataset_id = dataset_id

    try:
        ds = server.to_xarray()
    except requests.exceptions.HTTPError:
        return None

    ds = ds.swap_dims({'obs': 'precise_time'})
    ds = ds.reset_coords()
    for key in ['profile_id', 'time', 'longitude', 'latitude', 'trajectoryIndex', 'rowSize']:
        if key in ds.variables:
            ds = ds.drop_vars(key)

    ds = ds.squeeze(drop=True)
    rename = {k: v for k, v in _RENAME_MAP.items() if k in ds}
    ds = ds.rename(rename)
    ds = ds.sortby('time')
    return ds


def collect_glider(
    latitude: float,
    longitude: float,
    institution: str,
    variables: list[str] | None = None,
    extent: float = 2.7
) -> xr.Dataset:
    """
    Collect glider data from the GliderDAC near a given location.

    Navigation variables (time, latitude, longitude, depth, pressure)
    are always included. Science variables are controlled by the
    ``variables`` parameter. Call :func:`glider_variables` to see all
    valid names.

    The ``extent`` parameter controls the size of the bounding box. The
    default (2.7 nm) is appropriate for selecting observations near a
    specific mooring. Use a larger value (e.g., 200 nm) to capture all
    glider data across a whole array deployment area.

    :param latitude: Center latitude in decimal degrees
    :param longitude: Center longitude in decimal degrees
    :param institution: GliderDAC institution name. Use
        :func:`glider_institution` to look up the correct string from
        an OOI site designator.
    :param variables: Internal names of science variables to include.
        Defaults to all available science variables when None.
    :param extent: Half-width of the search bounding box in nautical
        miles (default 2.7)
    :return: Combined xarray Dataset sorted by time
    :raises ValueError: If any name in ``variables`` is not recognized
    """
    if variables is None:
        variables = GLIDER_SCIENCE_VARIABLES
    else:
        invalid = [v for v in variables if v not in _SCIENCE_VAR_MAP]
        if invalid:
            raise ValueError(
                f"Unknown variable(s): {invalid}. "
                f"Call glider_variables() to see valid names."
            )

    erddap_vars = _REQUIRED_ERDDAP_VARS + [_SCIENCE_VAR_MAP[v] for v in variables]
    bounding_box = _create_box(latitude, longitude, extent)
    gliders = _list_gliders(bounding_box, institution)

    dl = partial(_download_glider, bounding_box=bounding_box, erddap_vars=erddap_vars)
    with ThreadPoolExecutor(max_workers=5) as executor:
        frames = list(tqdm(
            executor.map(dl, gliders),
            total=len(gliders),
            desc='Downloading GliderDAC data',
            file=sys.stdout
        ))

    data = xr.concat([f for f in frames if f is not None], dim='time')
    data = data.sortby('time')
    return data
