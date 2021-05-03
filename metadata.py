# vim: set fileencoding=utf-8 :

import pathlib
import os
import collections

import pandas as pd

import xsdir
import ace

def fastNeutron(XSDIR, index):
    """
    Add metadata from all fast neutron ACE files.

    datapath: Where the data files are relative to
    XSDIR: XSDIR of DataFrame
    """
    path = pathlib.Path(datapath, XSDIR.loc[index].path)
    print(path)

    address = XSDIR.loc[index].address
    ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

    meta = {}
    XSDIR.loc[index, 'length'] = int(ACE.NXS[1])
    NE = int(ACE.NXS[3])
    XSDIR.loc[index, 'NE'] = NE
    XSDIR.loc[index, 'Emax'] = round(ACE.XSS[NE], 1)
    if (ACE.JXS[12] != 0) or (ACE.JXS[13] != 0):
        XSDIR.loc[index, 'GPD'] = True
    else:
        XSDIR.loc[index, 'GPD'] = False

    if ACE.JXS[2] != 0:
        if ACE._XSS[(ACE.JXS[2] - 1)] > 0:
            XSDIR.loc[index, 'nubar'] = 'nubar'
        else:
            XSDIR.loc[index, 'nubar'] = 'both'
    else:
        XSDIR.loc[index, 'nubar'] = 'no'

    # Charged particle  see XTM:96-200
    if ACE.NXS[7] > 0:
        XSDIR.loc[index, 'CP'] = True
    else:
        XSDIR.loc[index, 'CP'] = False
    # Delayed neutron
    if ACE.JXS[24] > 0:
        XSDIR.loc[index, 'DN'] = True
    else:
        XSDIR.loc[index, 'DN'] = False
    # Unresolvedresonance
    if ACE.JXS[23] > 0:
        XSDIR.loc[index, 'UR'] = True
    else:
        XSDIR.loc[index, 'UR'] = False

    # XSDIR = XSDIR.append(pd.Series(meta), ignore_index=True)
    # return meta

metadata = {
    'c': fastNeutron,
    'nc': fastNeutron
}

if __name__ == "__main__":
    print("Adding metadata to XSDIR DataFrame")

    global datapath
    datapath = pathlib.Path(os.environ['DATAPATH'])
    AWRs, XSDIR = xsdir.readXSDIR(pathlib.Path(datapath, 'xsdir'))

    # Additional columns for metadata
    metaColumns = {
        # FastNeutron
        'NE': int,
        'length': int,
        'Emax': float,
        'GPD': bool,
        'nubar': str,
        'CP': bool,
        'DN': bool,
        'UR': bool,
        # Thermal Scattering
    }

    # Add columns to DataFrame for metadata
    for name, dtype in metaColumns.items():
        if name not in XSDIR.columns:
            XSDIR[name] = pd.Series(dtype=dtype)

    lib_type = 'nc'
    libIndices = (XSDIR['lib_type'] == lib_type) & (XSDIR['ZA'] == 92235)
    entries = XSDIR.loc[libIndices]

    for index in entries.index:
        metadata.get(lib_type)(XSDIR, index)
