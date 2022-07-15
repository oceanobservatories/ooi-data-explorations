#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author Andrew Reed
@brief Check if a particular reference designator - data stream - parameter are
    in a given QC-lookup table
"""
import glob
import ast
import pandas as pd
import os


def parse_gross_range_table(table):
    """Parse the gross range qc-lookup table.
    
    Parameters
    ----------
    table: (pd.DataFrame)
        The gross range qc-lookup table loaded into python as a pandas dataframe.
        
    Returns
    -------
    table: (pd.DataFrame)
         The gross range qc-lookup table with the entries parsed into functional entries.
    """
    
    # Parse the table entries which are dictionaries
    table["parameters"] = table["parameters"].apply(lambda x: ast.literal_eval(x))
    table["qcConfig"] = table["qcConfig"].apply(lambda x: ast.literal_eval(x))

    # Next, parse out the inputs, fail_span, and suspect_span into their own columns
    table["inp"] = table["parameters"].apply(lambda x: x.get("inp"))
    table["fail_span"] = table["qcConfig"].apply(lambda x: x.get("qartod").get("gross_range_test").get("fail_span"))
    table["suspect_span"] = table["qcConfig"].apply(lambda x: x.get("qartod").get("gross_range_test").get("suspect_span"))

    return table


def parse_climatology_table(table):
    """Parse the climatology qc-lookup table.
    
    Parameters
    ----------
    table: (pd.DataFrame)
        The climatology qc-lookup table loaded into python as a pandas dataframe.
        
    Returns
    -------
    table: (pd.DataFrame)
         The climatology qc-lookup table with the entries parsed into functional entries.
    """
    table["parameters"] = table["parameters"].apply(lambda x: ast.literal_eval(x))
    table["inp"] = table["parameters"].apply(lambda x: x.get("inp"))
    return table


# +
def check_entry_gross_range(refdes, stream, variables, gross_range_test_values):
    """
    Check if the reference designator, stream, and variable(s) are in the gross range qc-lookup table.
    
    Parameters
    ----------
    refdes: (str)
        The reference designator of the instrument
    stream: (str)
        The particular data stream for the given reference designator
    variables: (str or list)
        The variable(s) for the reference designator-stream to check are in the qc-lookup table
    gross_range_test_values: (pd.DataFrame)
        The qc-lookup table for the gross range test which should contain the reference designator,
        stream, and variable(s) if they have been added
        
    Returns
    -------
    results: (dict)
        A dictionary with a boolean in "missing" key if the given variable for the reference
        designator - stream is in the qc-lookup table. Meant to be parsed into a pandas
        dataframe.
    """
    # Get the subsite - node - sensor
    subsite, node, sensor = refdes.split("-", 2)
    
    # Initialize a dictionary to store the results
    results = {
        "refdes": [],
        "stream": [],
        "variable": [],
        "missing": []
    }
    
    for var in variables:
        # Now find match in the gross range table
        match = gross_range_test_values[(gross_range_test_values["subsite"] == subsite) &
                                        (gross_range_test_values["node"] == node) &
                                        (gross_range_test_values["sensor"] == sensor) &
                                        (gross_range_test_values["stream"] == stream) &
                                        (gross_range_test_values["inp"] == var)]
        # Check if there is no match
        if len(match) == 0:
            results["refdes"].append(refdes)
            results["stream"].append(stream)
            results["variable"].append(var)
            results["missing"].append(True)
        else:
            results["refdes"].append(refdes)
            results["stream"].append(stream)
            results["variable"].append(var)
            results["missing"].append(False)
            
    return results


def check_entry_climatology(refdes, stream, variables, climatology_test_values, qartod_dir):
    """
    Check if the reference designator, data stream, and variables are in the climatology table.
    
    Parameters
    ----------
    refdes: (str)
        The reference designator of the instrument
    stream: (str)
        The particular data stream for the given reference designator
    variables: (str or list)
        The variable(s) for the reference designator-stream to check are in the qc-lookup table
    climatology_test_values: (pd.DataFrame)
        The qc-lookup table for the gross range test which should contain the reference designator,
        stream, and variable(s) if they have been added
    qartod_dir: (str)
        The path to where the climatology
        
    Returns
    -------
    results: (dict)
        A dictionary with the following keys:
        * refdes: reference designator
        * stream: stream
        * variable: variable
        * missing: boolean value indicating if the refdes-stream-variable is in the qc-lookup table
        * tableExists: if "missing" is False, a boolean which checks if the climatologyTable for the
                       variable is in the climatologyTable
        * tableKeyMatch: if "tableExists" is True, a boolean which checks if the key in the name of
                         the climatologyTable matches what is in the "inp" value from the qc-lookup
                         table
    """
    # Get the subsite - node - sensor
    subsite, node, sensor = refdes.split("-", 2)
    
    # Initialize a dictionary to store the results
    results = {
        "refdes": [],
        "stream": [],
        "variable": [],
        "missing": [], 
        "tableExists": [],
        "tableKeyMatch": []
    }
    
    for var in variables:
        # Now find match in the gross range table
        match = climatology_test_values[(climatology_test_values["subsite"] == subsite) &
                                        (climatology_test_values["node"] == node) &
                                        (climatology_test_values["sensor"] == sensor) &
                                        (climatology_test_values["stream"] == stream) &
                                        (climatology_test_values["inp"] == var)]

        # If the match is not missing, we then want to check that the climatology table exists
        if len(match) == 0:
            results["refdes"].append(refdes)
            results["stream"].append(stream)
            results["variable"].append(var)
            results["missing"].append(True)
            results["tableExists"].append(False)
            results["tableKeyMatch"].append(False)
        else:
            # Now check that the table is there
            for cind in match.index:
                tablePath = match.loc[cind, "climatologyTable"]
                # Check that the table exists
                tableExists = os.path.exists(f"{qartod_dir}/{tablePath}")
                if tableExists:
                    # Check that the tableKey matches the input key
                    tableName = tablePath.split("/")[1].split("-")[-1].split(".")[0]
                    tableKeyMatch = (tableName == var)
                else:
                    tableKeyMatch = False
                # Update the check table
                results["refdes"].append(refdes)
                results["stream"].append(stream)
                results["variable"].append(var)
                results["missing"].append(False)
                results["tableExists"].append(tableExists)
                results["tableKeyMatch"].append(tableKeyMatch)
    return results


def check_entries(refdes, stream, variables, qartod_dir=None, tests=["gross_range", "climatology"]):
    """
    Check if the reference designator, data stream, and variables are in the qc-lookup tables.
       
    Current iteration is set for "gross_range" and "climatology" tests. Designed to be flexible
    enough to accommodate future tests so long as they follow the qc-lookup format set by the
    existing test tables.
    
    Parameters
    ----------
    refdes: (str)
        The reference designator of the instrument
    stream: (str)
        The particular data stream for the given reference designator
    variables: (str or list)
        The variable(s) for the reference designator-stream to check are in the qc-lookup table
    qartod_dir: (str)
        The path to where qc-lookup tables for the qartod tests are located
    tests: (list -> str: ["gross_range", "climatology"])
        A list of qartod tests which to check the qc-lookup value tables for the relevant 
        reference designator - stream - variable entries.
        
    Returns
    -------
    gross_range_results: (dict)
        A dictionary with a boolean in "missing" key if the given variable for the reference
        designator - stream is in the qc-lookup table. Meant to be parsed into a pandas
        dataframe.

    climatology_results: (dict)
        A dictionary with the following keys:
        * refdes: reference designator
        * stream: stream
        * variable: variable
        * missing: boolean value indicating if the refdes-stream-variable is in the qc-lookup table
        * tableExists: if "missing" is False, a boolean which checks if the climatologyTable for the
                       variable is in the climatologyTable
        * tableKeyMatch: if "tableExists" is True, a boolean which checks if the key in the name of
                         the climatologyTable matches what is in the "inp" value from the qc-lookup
                         table
    """
    # Get the subsite, node, sensor
    subsite, node, sensor = refdes.split("-", 2)
    inst_class = sensor.split("-")[-1][0:5].lower()

    # Find the directory where the qartod tables are stored
    if qartod_dir is None:
        qartod_dir = "/home/"
        qartod_dir = glob.glob(qartod_dir + "**/" + "qartod/" + inst_class, recursive=True)[0]
        
    # Load the test tables
    for test in tests:
        if "gross_range" in test:
            gross_range_test_table = pd.read_csv(f"{qartod_dir}/{inst_class}_qartod_gross_range_test_values.csv")
            gross_range_test_table = parse_gross_range_table(gross_range_test_table)
        elif "climatology" in test:
            climatology_test_table = pd.read_csv(f"{qartod_dir}/{inst_class}_qartod_climatology_test_values.csv")
            climatology_test_table = parse_climatology_table(climatology_test_table)
        else:
            raise FileExistsError(f"{test} not yet available.")
    
    # Check if there is more than one variable
    if type(variables) is not list and type(variables) == str:
        variables = variables.split(",")
        
    for test in tests:
        if "gross_range" in test:
            gross_range_results = check_entry_gross_range(refdes, stream, variables, gross_range_test_table)
        elif "climatology" in test:
            climatology_results = check_entry_climatology(refdes, stream, variables, climatology_test_table, qartod_dir)

    # Return the results
    return gross_range_results, climatology_results
