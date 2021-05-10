# vim: set fileencoding=utf-8 :

"""
This module will parse an xsdir file.
"""

import io
import re
import pathlib
import os

import pandas as pd
import numpy as np

_date_pattern = re.compile(r'\s*\d{2}\/\d{2}\/\d{2,4}\s*')

def _addExtras(XSDIR):
    """
    addLibType will append extra information to the XSDIR DataFrame. This extra
    information is extracted/calculated from the content already in the
    DataFrame.
    """
    libTypes = []
    ZAs = []
    libraries = []

    for entry in XSDIR.itertuples():
        zaid = entry.ZAID
        za = zaid.split('.')[0]
        try:
            za = int(za)
        except ValueError:
            za = za
        ZAs.append(za)
        lType = zaid[-1]
        if lType == 'c':
            if zaid[-2] == 'n':
                lType = 'nc'
        libTypes.append(lType)

        path = entry.path
        l = len(path.parents)
        if l < 2:
            libraries.append(entry.path)
        else:
            libraries.append(entry.path.parents[l-2])

    XSDIR['library'] = libraries
    XSDIR['lib_type'] = libTypes
    XSDIR['ZA'] = ZAs
    XSDIR['T(K)'] = round(XSDIR['temperature']/8.6173E-11, 1)

def readXSDIR(filename=pathlib.Path(os.environ['DATAPATH'], 'xsdir')):
    """
    readXSDIR will read the XSDIR and return two pandas.DataFrame objects;

        1. atomic weight ratios
        2. xsdir entries
    """
    AWRs = []
    entries = []

    columnType = {
        "ZAID":"U",
        "AWR":float,
        "path":"U",
        "access":int,
        "file_type":int,
        "address":int,
        "table_length":int,
        "record_length":int,
        "num_entries":int,
        "temperature":float,
        "ptable":bool,
    }

    with filename.open('r') as xsdirFile:
        for line in xsdirFile:
            if line.strip() == 'atomic weight ratios':
                break

        # Parse atomic weight ratios
        for line in xsdirFile:
            if _date_pattern.match(line):
                break
            AWRs.extend(line.split())

        # Make sure we are in the directory listing
        for line in xsdirFile:
            if line.strip() == 'directory':
                break

        # Parse entries
        lines = ""
        for line in xsdirFile:
            # Check if entry extends to next line
            if line.strip().endswith('+'):
                continuationIndex = line.rfind("+")
                line = line[:continuationIndex] + xsdirFile.readline()

            lines += line
    AWRs = pd.DataFrame(np.reshape(AWRs, (-1, 2)), columns=["ZA", "AWR"]) \
             .astype({"ZA":int, "AWR":
                                                            float})
    entries = pd.read_csv(io.StringIO(lines), sep='\s+', names=list(columnType)) \
                .fillna(0) \
                .astype(columnType)
    entries.path = entries.path.apply(pathlib.Path)
    _addExtras(entries)
    return (AWRs, entries)

if __name__ == "__main__":
    print("\nI'm parsing XSDIR file(s).\n")

    AWRs, entries = readXSDIR()

    # xs = xsdir()
