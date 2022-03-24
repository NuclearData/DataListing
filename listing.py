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
            'h': self._chargedParticle, # hydrogen/proton
            'o': self._chargedParticle, # deuteron
            'r': self._chargedParticle, # triton
            's': self._chargedParticle, # helion/He-3
            'a': self._chargedParticle, # alpha/He-4
        }

        AWRs, self.XSDIR = xsdir.readXSDIR(xsdirPath)

    def _fastNeutron(self, entry, ACE):
        """
        Add metadata from all fast neutron ACE files.

        entry: Row of pandas data frame object
        """
        meta = {}
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

    def _thermalScattering(self, entry, ACE):
        """
        Add metadata from thermal scattering ACE files.

        entry: Row of pandas data frame object
        """
        meta = {}
        meta['NA'] = int(ACE.NXS[3] + 1)
        meta['NE'] = int(ACE.NXS[4])
        meta['IZ'] = ACE._IZ

        if ACE.NXS[7] == 2:
            meta['representation'] = 'continuous'
        else:
            meta['representation'] = 'discrete'

        return meta

    def _photon(self, entry, ACE):
        """
        Add metadata from photon ACE files.

        entry: Row of pandas data frame object
        """
        meta = {}
        meta['NE'] = int(ACE.NXS[3])
        return meta

    def _chargedParticle(self, entry, ACE):
        """
        Add metadata from the charged-particle ACE files.

        entry: Row of pandas data frame object
        """

        meta = {}
        NE = int(ACE.NXS[3])
        meta['target'] = entry.ZA
        meta[ 'NE'] = NE
        meta[ 'Emax'] = round(ACE.XSS[NE], 1)
        return meta

    def _default(self, entry, ACE):
        """
        _default does nothing, but prevents Python from crashing when
        """
        return {}

    def extend(self, index):
        """
        extend will add metadata to a row of self.XSDIR given the row's index
        """
        entry = self.XSDIR.loc[index]
        path = pathlib.Path(self.datapath, entry.path)

        address = entry.address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {'ZAID': entry.ZAID}

        meta['length'] = int(ACE.NXS[1])
        meta['Date'] = ACE.process_date

        meta2 = self.metadataFunctions.get(
            entry.lib_type, self._default)(entry, ACE)

        meta.update(meta2)
        return meta


class DisplayData:
    def __init__(self, XSDIR, lib_type=None):

        if lib_type:
            self.XSDIR = XSDIR.query('lib_type == @lib_type')
        else:
            self.XSDIR = XSDIR
        self.lib_type = lib_type

        self.displayColumns = {
            'c':  ['ZAID', 'AWR', 'library', 'path','ZA', 'T(K)', 'Date', 'NE',
                       'Emax', 'GPD', 'nubar', 'CP', 'DN', 'UR'],
            'nc': ['ZAID', 'AWR', 'library',  'path','ZA', 'T(K)', 'Date', 'NE',
                       'Emax', 'GPD', 'nubar', 'CP', 'DN', 'UR'],
            't': ['ZAID', 'library', 'path','ZA', 'T(K)', 'Date', 'NE', 'NA',
                  'representation'],
            'h': ['ZAID', 'AWR', 'library', 'path', 'ZA', 'T(K)', 'Date', 'NE',
                  'Emax'],
            'o': ['ZAID', 'AWR', 'library', 'path', 'ZA', 'T(K)', 'Date', 'NE',
                  'Emax'],
            'r': ['ZAID', 'AWR', 'library', 'path', 'ZA', 'T(K)', 'Date', 'NE',
                  'Emax'],
            's': ['ZAID', 'AWR', 'library', 'path', 'ZA', 'T(K)', 'Date', 'NE',
                  'Emax'],
            'a': ['ZAID', 'AWR', 'library', 'path', 'ZA', 'T(K)', 'Date', 'NE',
                  'Emax'],
            None: ['ZAID', 'AWR', 'library', 'path', 'ZA', 'T(K)', 'Date'],
        }

    def __call__(self, ZA=None, columns=[]):
        """
        """
        if isinstance(self.lib_type, list):
            lt = self.lib_type[0]
        else:
            lt = self.lib_type

        defaultColumns = ['ZAID', 'AWR', 'library', 'path','ZA', 'T(K)', 'Date']
        if not columns:
            columns = self.displayColumns.get(lt, defaultColumns)

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

    # ddir.XSDIR = ddir.XSDIR.query('ZA == 1001 or ZA == "lwtr"')
    with mp.Pool(N) as pool:
        results = list(
            tqdm(pool.imap(ddir.extend, ddir.XSDIR.index), total=len(ddir.XSDIR)))

    dtype = {
        'Date': str,

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

        # Charged-particle
        'target': int,
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

    # Get default XSDIR
    defaultXSDIR = max(
        pathlib.Path(os.environ['DATAPATH']).glob("xsdir*"),
        key = lambda p: p.stat().st_ctime
    )

    parser.add_argument('--xsdir', type=pathlib.Path,
        default = defaultXSDIR,
        help="Path to xsdir file. Defaults to $DATAPATH/xsdir")

    parser.add_argument('-N', type=int,
        default=max(1, mp.cpu_count()-1),
        help="Number of parallel threads.")

    parser.add_argument('--dont-generate', action='store_true', default=False,
        help="Don't generate {}".format(xsdirName))

    return parser.parse_args()


if __name__ == "__main__":
    __spec__ = None

    args = processInput()
    if not args.dont_generate:
        generateJSON(args.xsdir, args.N)

    XSDIR = loadXSDIR()

