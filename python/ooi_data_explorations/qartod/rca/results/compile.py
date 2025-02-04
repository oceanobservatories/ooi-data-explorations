#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob

read_files = glob.glob("*gross_range.csv.fixed")

with open("grossRange.csv", "wb") as outfile:
    for f in read_files:
        with open(f, "rb") as infile:
            outfile.write(infile.read())
