#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Christopher Wingard
@brief Provides access to the IOOS GliderDAC (GDAC) via ERDDAP for use in
    creating QARTOD test values for OOI glider instrument classes across the
    Coastal Endurance, Coastal Pioneer, and Global arrays.
"""
import dask
import io
import numpy as np
import pandas as pd
import re
import requests
import time
import xarray as xr

from dask.diagnostics import ProgressBar
from erddapy import ERDDAP

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

# Maps internal science variable names to their GliderDAC QARTOD primary flag
# names. Only variables for which QC flags are currently available in the
# GliderDAC are listed.
_QC_FLAG_MAP: dict[str, str] = {
    'conductivity': 'qartod_conductivity_primary_flag',
    'temperature': 'qartod_temperature_primary_flag',
    'pressure': 'qartod_pressure_primary_flag',
    'salinity': 'qartod_salinity_primary_flag',
    'density': 'qartod_density_primary_flag',
}

# Public list of valid science variable names (internal/post-rename names)
GLIDER_SCIENCE_VARIABLES: list[str] = list(_SCIENCE_VAR_MAP.keys())

# Maps OOI site prefix to GliderDAC institution name.
GLIDER_INSTITUTIONS: dict[str, str] = {
    'CE': 'ooi_coastal_endurance',
    'CP': 'ooi_coastal_global_scale_nodes_cgsn_',
    'GA': 'ooi_coastal_global_scale_nodes_cgsn_',
    'GI': 'ooi_coastal_global_scale_nodes_cgsn_',
    'GP': 'ooi_coastal_global_scale_nodes_cgsn_',
    'GS': 'ooi_coastal_global_scale_nodes_cgsn_',
}

# Maps OOI site prefix to the lowercase dataset ID prefix used in GliderDAC.
# Pioneer and all Global sites share an institution string, so the dataset ID
# prefix is required to disambiguate them in ERDDAP searches.
GLIDER_DATASET_PREFIXES: dict[str, str] = {
    'CE': 'ce',
    'CP': 'cp',
    'GA': 'ga',
    'GI': 'gi',
    'GP': 'gp',
    'GS': 'gs',
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


def glider_search_for(site: str) -> str:
    """
    Return the ERDDAP full-text search string for delayed-mode glider
    datasets from an OOI site.

    Combines ``'delayed'`` with a dataset ID prefix so that Pioneer and
    Global sites -- which share an institution string -- are correctly
    disambiguated.

    :param site: OOI site designator (e.g. 'CE05MOAS', 'CP05MOAS')
    :return: ERDDAP search_for string (e.g. 'delayed datasetID=ce')
    :raises ValueError: If the site prefix has no known dataset prefix
        mapping
    """
    prefix = site[:2].upper()
    if prefix not in GLIDER_DATASET_PREFIXES:
        raise ValueError(
            f"No GliderDAC dataset prefix known for site prefix '{prefix}'. "
            f"Known prefixes: {sorted(GLIDER_DATASET_PREFIXES.keys())}. "
            f"Update GLIDER_DATASET_PREFIXES in gdac.py to add it."
        )
    return f"delayed datasetID={GLIDER_DATASET_PREFIXES[prefix]}"


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


def _list_gliders(
    institution: str,
    search_for: str,
    bounding_box: list[float] | None = None
) -> np.ndarray:
    """
    List all delayed-mode glider datasets for an institution in the
    GliderDAC, optionally filtered to a bounding box.

    :param institution: GliderDAC institution name (see
        :data:`GLIDER_INSTITUTIONS` or :func:`glider_institution`)
    :param search_for: ERDDAP full-text search string (see
        :func:`glider_search_for`)
    :param bounding_box: Optional [min_lat, max_lat, min_lon, max_lon].
        When None, all datasets for the institution are returned.
    :return: array of GliderDAC dataset IDs
    """
    server = ERDDAP(server='ngdac')
    kwargs: dict = dict(response='csv', search_for=search_for, institution=institution)
    if bounding_box is not None:
        kwargs.update(
            min_lat=bounding_box[0],
            max_lat=bounding_box[1],
            min_lon=bounding_box[2],
            max_lon=bounding_box[3],
        )
    search_url = server.get_search_url(**kwargs)
    search = pd.read_csv(search_url)
    return search['Dataset ID'].values


def _stratified_subsample(gliders: np.ndarray, subsample: int) -> np.ndarray:
    """
    Subsample glider dataset IDs with stratification by deployment date
    to ensure full temporal coverage.

    Deployment dates are parsed from the dataset ID (e.g.
    ``ce_388-20160718T1600``). The full date range is divided into
    ``n_target`` equal-sized buckets and one dataset is drawn at random
    from each bucket, guaranteeing uniform coverage across the record.
    Glider diversity follows naturally from temporal coverage since
    different gliders operated at different times.

    :param gliders: Array of GliderDAC dataset IDs.
    :param subsample: Percentage of datasets to retain [1, 100].
    :return: Subsampled array, length <= ``n_target``.
    """
    n_target = max(1, round(len(gliders) * subsample / 100))
    if n_target >= len(gliders):
        return gliders

    date_pat = re.compile(r'-(\d{8}T\d{4})')
    records = []
    for gid in gliders:
        m = date_pat.search(gid)
        date = pd.to_datetime(m.group(1), format='%Y%m%dT%H%M') if m else pd.NaT
        records.append({'dataset_id': gid, 'date': date})

    df = pd.DataFrame(records).sort_values('date', na_position='last').reset_index(drop=True)

    # Always anchor to the earliest and latest deployments, then fill
    # n_target - 2 slots by stratified random sampling from the middle.
    first = df.iloc[0]['dataset_id']
    last = df.iloc[-1]['dataset_id']

    if n_target <= 2:
        anchors = [first] if n_target == 1 else [first, last]
        return np.array(anchors)

    middle = df.iloc[1:-1]
    n_middle = n_target - 2
    bucket_size = len(middle) / n_middle
    rng = np.random.default_rng()
    selected: list[str] = [first]

    for i in range(n_middle):
        lo = int(i * bucket_size)
        hi = max(lo + 1, int((i + 1) * bucket_size))
        pick = middle.iloc[lo:hi].sample(1, random_state=int(rng.integers(1_000_000)))
        selected.append(pick['dataset_id'].iloc[0])

    selected.append(last)
    return np.array(selected)


def _download_glider(
    dataset_id: str,
    bounding_box: list[float] | None,
    erddap_vars: list[str],
    max_retries: int = 3,
    timeout: int = 300
) -> xr.Dataset | None:
    """
    Download a single glider dataset from the GliderDAC.

    Creates a fresh ERDDAP connection per call so that concurrent
    downloads do not share mutable server state. Retries on timeout or
    HTTP errors with linear backoff between attempts.

    :param dataset_id: GliderDAC dataset ID
    :param bounding_box: Optional [min_lat, max_lat, min_lon, max_lon]
        to constrain the downloaded data. When None, all data for the
        dataset are downloaded.
    :param erddap_vars: ERDDAP variable names to request
    :param max_retries: Number of download attempts before giving up
        (default 3)
    :param timeout: Per-request read timeout in seconds (default 300)
    :return: xarray Dataset or None if all attempts fail
    """
    # set up the erddapy server object and create the dataset request URL
    server = ERDDAP(server='ngdac')
    if bounding_box is not None:
        server.constraints = {
            'latitude>=': bounding_box[0],
            'latitude<=': bounding_box[1],
            'longitude>=': bounding_box[2],
            'longitude<=': bounding_box[3],
        }
    server.protocol = 'tabledap'
    server.variables = erddap_vars
    server.dataset_id = dataset_id
    url = server.get_download_url(response='csv')

    # download the data into a pandas dataframe (downloading CSV content as it is quicker than NetCDF)
    df = None
    units: dict[str, str] = {}
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=timeout)
            response.raise_for_status()
            buf = io.StringIO(response.text)
            columns = pd.read_csv(buf, nrows=0).columns.tolist()  # pull the column names
            buf.seek(0)
            units_vals = pd.read_csv(buf, header=None, skiprows=1, nrows=1).iloc[0].tolist()  # pull the unit strings
            units = dict(zip(columns, units_vals))
            buf.seek(0)
            df = pd.read_csv(buf, skiprows=[1])  # skip units row, keep header, pull the data
            time_col = next(
                (c for c in ['precise_time', 'time'] if c in df.columns), None
            )
            if time_col is None:
                return None
            df[time_col] = pd.to_datetime(df[time_col], utc=True).dt.tz_convert(None)
            df = df.set_index(time_col)
            break
        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                time.sleep(10 * (attempt + 1))

    if df is None or df.empty:
        return None

    # convert to a xarray dataset mapping the unit strings to the variables
    ds = df.to_xarray()
    for key in ['trajectoryIndex', 'rowSize']:
        if key in ds:
            ds = ds.drop_vars(key)
    units_renamed = {_RENAME_MAP.get(k, k): v for k, v in units.items()}
    for name in list(ds.data_vars) + list(ds.coords):
        if name in units_renamed:
            ds[name].attrs['units'] = units_renamed[name]
    rename = {k: v for k, v in _RENAME_MAP.items() if k in ds}
    ds = ds.rename(rename)

    # add the glider serial number as a variable
    glider_num = int(dataset_id.split('-')[0].rsplit('_', 1)[-1])
    ds['glider'] = xr.DataArray(
        np.full(ds.sizes['time'], glider_num, dtype=int), dims=['time']
    )

    # remove the fill values (helps to minimize dataset size in memory
    ds = ds.set_coords(['longitude', 'latitude', 'glider'])
    ds = ds.dropna('time', how='all')
    if ds.sizes['time'] == 0:
        return None
    ds = ds.sortby('time')

    return ds


def collect_glider(
    institution: str,
    search_for: str,
    latitude: float | list[float] | None = None,
    longitude: float | list[float] | None = None,
    variables: list[str] | None = None,
    extent: float = 2.7,
    qc_flags: bool = False,
    subsample: int | None = None
) -> xr.Dataset:
    """
    Collect glider data from the GliderDAC for an institution, with
    optional spatial filtering.

    Navigation variables (time, latitude, longitude, depth, pressure)
    are always included. Science variables are controlled by the
    ``variables`` parameter. Call :func:`glider_variables` to see all
    valid names.

    Spatial filtering is controlled by ``latitude`` and ``longitude``:

    - Pass scalar values with an ``extent`` to build a bounding box
      centered on that point. The default extent (2.7 nm) is suited to
      selecting data near a specific mooring.
    - Pass two-element lists ``[min, max]`` to use an explicit bounding
      box.
    - Pass ``None`` (the default) for no spatial filter; all delayed-mode
      datasets for the institution are returned.

    :param institution: GliderDAC institution name. Use
        :func:`glider_institution` to look up the correct string from
        an OOI site designator.
    :param search_for: ERDDAP full-text search string. Use
        :func:`glider_search_for` to derive the correct value from an
        OOI site designator.
    :param latitude: Center latitude (scalar), explicit [min, max]
        bounds, or None for no spatial filter.
    :param longitude: Center longitude (scalar), explicit [min, max]
        bounds, or None for no spatial filter.
    :param variables: Internal names of science variables to include.
        Defaults to all available science variables when None.
    :param extent: Half-width of the bounding box in nautical miles
        when ``latitude`` and ``longitude`` are scalars (default 2.7).
        Ignored when they are lists or None.
    :param qc_flags: When True, also download the QARTOD primary summary
        flag for each requested variable that has one available in the
        GliderDAC. Flag variables are retained in the returned dataset
        under their GliderDAC names (e.g.
        ``qartod_temperature_primary_flag``). Currently only CTD
        variables (temperature, conductivity, salinity, density) have
        corresponding flags; other variables are silently skipped.
    :param subsample: If set, randomly select this percentage of the
        available datasets rather than downloading all of them. Must be
        an integer in [1, 100]. Useful for large arrays where a
        representative sample is sufficient. Defaults to None (all
        datasets downloaded).
    :return: Combined xarray Dataset sorted by time
    :raises ValueError: If any name in ``variables`` is not recognized,
        or if ``subsample`` is not in the range [1, 100]
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
    if qc_flags:
        erddap_vars += [_QC_FLAG_MAP[v] for v in variables if v in _QC_FLAG_MAP]

    if latitude is None or longitude is None:
        bounding_box = None
    elif isinstance(latitude, list) and isinstance(longitude, list):
        bounding_box = [latitude[0], latitude[1], longitude[0], longitude[1]]
    else:
        bounding_box = _create_box(float(latitude), float(longitude), extent)

    gliders = _list_gliders(institution, search_for, bounding_box)

    if subsample is not None:
        if not 1 <= subsample <= 100:
            raise ValueError(
                f"subsample must be an integer in [1, 100], got {subsample}."
            )
        gliders = _stratified_subsample(gliders, subsample)

    downloads = [
        dask.delayed(_download_glider)(gid, bounding_box, erddap_vars)
        for gid in gliders
    ]
    with ProgressBar():
        frames = dask.compute(*downloads, scheduler='threads', num_workers=3)

    data = xr.concat([f for f in frames if f is not None], dim='time')
    data = data.sortby(['glider', 'time'])
    return data
