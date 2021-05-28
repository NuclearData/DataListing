# vim: set fileencoding=utf-8 :
#! /usr/bin/env python3

import pathlib
import os
import collections
import multiprocessing as mp
import argparse
import textwrap

import pandas as pd
import IPython.display
from tqdm.auto import tqdm

import xsdir
import ace

xsdirName = "xsdir.json"

class DataDirectory:
    def __init__(self, xsdirPath):
        self.problems = []
        self.datapath = xsdirPath.parent

        self.metadataFunctions = {
            'c': self._fastNeutron,
            'nc': self._fastNeutron,
            't': self._thermalScattering,
            'p': self._photon,
        }

        AWRs, self.XSDIR = xsdir.readXSDIR(xsdirPath)

    def _fastNeutron(self, entry):
        """
        Add metadata from all fast neutron ACE files.

        index: Location in self.XSDIR of row where we are adding metadata
        """
        path = pathlib.Path(self.datapath, entry.path)

        address = entry.address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {'ZAID': entry.ZAID}
        meta[ 'length'] = int(ACE.NXS[1])
        NE = int(ACE.NXS[3])
        meta[ 'NE'] = NE
        meta[ 'Emax'] = round(ACE.XSS[NE], 1)
        if (ACE.JXS[12] != 0) or (ACE.JXS[13] != 0):
            meta[ 'GPD'] = True
        else:
            meta[ 'GPD'] = False

        if ACE.JXS[2] != 0:
            if ACE._XSS[(ACE.JXS[2] - 1)] > 0:
                meta[ 'nubar'] = 'nubar'
            else:
                meta[ 'nubar'] = 'both'
        else:
            meta[ 'nubar'] = 'no'

        # Charged particle  see XTM:96-200
        if ACE.NXS[7] > 0:
            meta[ 'CP'] = True
        else:
            meta[ 'CP'] = False
        # Delayed neutron
        if ACE.JXS[24] > 0:
            meta[ 'DN'] = True
        else:
            meta[ 'DN'] = False
        # Unresolvedresonance
        if ACE.JXS[23] > 0:
            meta[ 'UR'] = True
        else:
            meta[ 'UR'] = False

        return meta

    def _thermalScattering(self, entry):
        """
        Add metadata from thermal scattering ACE files.

        index: Location in self.XSDIR of row where we are adding metadata
        """
        path = pathlib.Path(self.datapath, entry.path)

        address = entry.address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {'ZAID': entry.ZAID}
        meta['NA'] = int(ACE.NXS[3] + 1)
        meta['NE'] = int(ACE.NXS[4])

        if ACE.NXS[7] == 2:
            meta['representation'] = 'continuous'
        else:
            meta['representation'] = 'discrete'

        return meta

    def _photon(self, entry):
        """
        Add metadata from photon ACE files.

        index: Location in self.XSDIR of row where we are adding metadata
        """
        path = pathlib.Path(self.datapath, entry.path)

        address = entry.address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {'ZAID': entry.ZAID}

        meta['length'] = int(ACE.NXS[1])
        meta['NE'] = int(ACE.NXS[3])
        return meta

    def _default(self, entry):
        """
        _default does nothing, but prevents Python from crashing when
        """
        path = pathlib.Path(self.datapath, entry.path)

        address = entry.address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {'ZAID': entry.ZAID}

        meta['length'] = int(ACE.NXS[1])
        return meta

    def extend(self, index):
        """
        extend will add metadata to a row of self.XSDIR given the row's index
        """
        entry = self.XSDIR.loc[index]
        return self.metadataFunctions.get(
            entry.lib_type, self._default)(entry)


class DisplayData:
    def __init__(self, XSDIR, lib_type=None):

        if lib_type:
            self.XSDIR = XSDIR.query('lib_type == @lib_type')
        else:
            self.XSDIR = XSDIR
        self.lib_type = lib_type

        self.displayColumns = {
            'c':  ['ZAID', 'AWR', 'library', 'ZA', 'T(K)', 'NE',
                       'Emax', 'GPD', 'nubar', 'CP', 'DN', 'UR'],
            'nc': ['ZAID', 'AWR', 'library', 'ZA', 'T(K)', 'NE',
                       'Emax', 'GPD', 'nubar', 'CP', 'DN', 'UR'],
            't': ['ZAID', 'library', 'ZA', 'T(K)', 'NE', 'NA', 'representation'],
            'h': ['ZAID', 'AWR', 'library', 'ZA', 'T(K)'],
            'p': ['ZAID', 'AWR', 'library', 'ZA', 'T(K)'],
            'm': ['ZAID', 'AWR', 'library', 'ZA', 'T(K)'],
            'y': ['ZAID', 'AWR', 'library', 'ZA', 'T(K)'],
            None: self.XSDIR.columns
        }

    def __call__(self, ZA=None, columns=[]):
        """
        """
        if isinstance(self.lib_type, list):
            lt = self.lib_type[0]
        else:
            lt = self.lib_type

        if not columns:
            columns = self.displayColumns.get(lt, self.XSDIR.columns)

        if ZA:
            XSDIR = self.XSDIR.query('ZA == @ZA')
        else:
            XSDIR = self.XSDIR

        IPython.display.display(XSDIR[[*columns]])

    def _default(self, ZA=None, columns=[]):
        """
        The display function when lib_type doesn't exist
        """
        if columns:
            IPython.display.display(self.XSDIR[[*columns]])
        else:
            IPython.display.display(self.XSDIR)


def loadXSDIR(filename=xsdirName):
    """
    loadXSDIR will create a pandas DataFrame from a json file on disk. It will
    first check to see if the filename exists
    """
    path = pathlib.Path(filename)

    if path.exists():
        return pd.read_json(path)
    else:
        print(textwrap.dedent("""
            Can't read XSDIR information from file:
            \t{} 
            as it doens't exist""").format(path))

        print("Please first run: python listing.py to generate {}"
              .format(xsdirName))

def generateJSON(xsdirPath, N=max(1, mp.cpu_count()-1)):
    """
    generateJSON will generate the JSON version of the XSDIR pandas DataFrame.
    It saves the fil
    """
    ddir = DataDirectory(xsdirPath)

    # results = process_map(ddir.extend, ddir.XSDIR.index, max_workers=N,
    #                       chunksize=1)
    with mp.Pool(N) as pool:
        results = list(
            tqdm(pool.imap(ddir.extend, ddir.XSDIR.index), total=len(ddir.XSDIR)))

    dtype = {
        # Continuous-energy neutron
        'length': int,
        'NE': int,
        'Emax': float,
        'GPD': bool,
        'nubar': str,
        'CP': bool,
        'DN': bool,
        'UR': bool,

        # Thermal Scattering
        'NA': int,
        'representation': str,
    }
    results = pd.DataFrame([r for r in results if r]).fillna(0).astype(dtype)

    ddir.XSDIR = pd.merge(ddir.XSDIR, results, on='ZAID')

    with open(xsdirName, 'w') as jsonFile:
        json = ddir.XSDIR.to_json(orient='records', default_handler=str, indent=2)
        jsonFile.write(json)

    if ddir.problems:
        print("There were problems reading data from these ZAIDs:")
        for ZAID in ddir.problems:
            print("\t{}".format(ZAID))

    print("Finished generating JSON from XSDIR.")


def processInput():
    description= "Preparing to list available ACE data"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--xsdir', type=pathlib.Path,
        default = pathlib.Path(os.environ['DATAPATH'], 'xsdir'),
        help="Path to xsdir file. Defaults to $DATAPATH/xsdir")

    parser.add_argument('-N', type=int,
        default=max(1, mp.cpu_count()-1),
        help="Number of parallel threads.")

    parser.add_argument('--dont-generate', action='store_true', default=False,
        help="Don't generate {}".format(xsdirName))

    return parser.parse_args()


if __name__ == "__main__":

    args = processInput()
    if not args.dont_generate:
        generateJSON(args.xsdir, args.N)

    XSDIR = loadXSDIR()
