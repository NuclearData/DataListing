# vim: set fileencoding=utf-8 :
#! /usr/bin/env python3

import pathlib
import os
import collections
import multiprocessing
import argparse
import textwrap

import pandas as pd
import IPython.display

import xsdir
import ace

xsdirName = "xsdir.json"

class DataDirectory:
    def __init__(self, xsdirPath):
        self.datapath = xsdirPath.parent

        self.metadataFunctions = {
            'c': self._fastNeutron,
            'nc': self._fastNeutron
        }

        AWRs, self.XSDIR = xsdir.readXSDIR(xsdirPath)

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
            if name not in self.XSDIR.columns:
                self.XSDIR[name] = pd.Series(dtype=dtype)

        # libIndices = (self.XSDIR['lib_type'] == 'nc') & \
        #              (self.XSDIR['ZA'] == 1001)
        # self.XSDIR = self.XSDIR.loc[libIndices]

    def _fastNeutron(self, index):
        """
        Add metadata from all fast neutron ACE files.

        index: Location in self.XSDIR of row where we are adding metadata
        """
        path = pathlib.Path(self.datapath, self.XSDIR.loc[index].path)
        print("{}\t{}".format(self.XSDIR.loc[index].ZAID, path))

        address = self.XSDIR.loc[index].address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {}
        self.XSDIR.loc[index, 'length'] = int(ACE.NXS[1])
        NE = int(ACE.NXS[3])
        self.XSDIR.loc[index, 'NE'] = NE
        self.XSDIR.loc[index, 'Emax'] = round(ACE.XSS[NE], 1)
        if (ACE.JXS[12] != 0) or (ACE.JXS[13] != 0):
            self.XSDIR.loc[index, 'GPD'] = True
        else:
            self.XSDIR.loc[index, 'GPD'] = False

        if ACE.JXS[2] != 0:
            if ACE._XSS[(ACE.JXS[2] - 1)] > 0:
                self.XSDIR.loc[index, 'nubar'] = 'nubar'
            else:
                self.XSDIR.loc[index, 'nubar'] = 'both'
        else:
            self.XSDIR.loc[index, 'nubar'] = 'no'

        # Charged particle  see XTM:96-200
        if ACE.NXS[7] > 0:
            self.XSDIR.loc[index, 'CP'] = True
        else:
            self.XSDIR.loc[index, 'CP'] = False
        # Delayed neutron
        if ACE.JXS[24] > 0:
            self.XSDIR.loc[index, 'DN'] = True
        else:
            self.XSDIR.loc[index, 'DN'] = False
        # Unresolvedresonance
        if ACE.JXS[23] > 0:
            self.XSDIR.loc[index, 'UR'] = True
        else:
            self.XSDIR.loc[index, 'UR'] = False

    def _default(self, index):
        """
        _default does nothing, but prevents Python from crashing when 
        """
        pass

    def extend(self, index):
        """
        extend will add metadata to a row of self.XSDIR given the row's index
        """
        lib_type = self.XSDIR.loc[index].lib_type

        self.problems = []
        try:
            self.metadataFunctions.get(lib_type, self._default)(index)
        except Exception as e:
            print("Problem with {}".format(self.XSDIR.loc[index].ZAID))
            self.problems.append(index)


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


class DisplayData:
    def __init__(self, XSDIR, lib_type=None):

        if lib_type:
            self.XSDIR = XSDIR.query('lib_type == @lib_type')
        else:
            self.XSDIR = XSDIR
        self.lib_type = lib_type

        self._display = {
            'c': self._fastNeutron,
            'nc': self._fastNeutron,
            None: self._default
        }

    def __call__(self, ZA=None, columns=[]):
        """
        """
        if isinstance(self.lib_type, list):
            lt = self.lib_type[0]
        else:
            lt = self.lib_type

        self._display.get(lt, self._default)(ZA, columns)

    def _default(self, ZA=None, columns=[]):
        """
        The display function when lib_type doesn't exist
        """
        if columns:
            IPython.display.display(self.XSDIR[[*columns]])
        else:
            IPython.display.display(self.XSDIR)

    def _fastNeutron(self, ZA=None, columns=[]):
        """
        The display function when asking for a continuous-energy neutron
        material
        """
        if not columns:
            columns = ['ZAID', 'AWR', 'library', 'ZA', 'T(K)', 'NE',
                       'Emax', 'GPD', 'nubar', 'CP', 'DN', 'UR']

        if ZA:
            XSDIR = self.XSDIR.query('ZA == @ZA')
        else:
            XSDIR = self.XSDIR

        IPython.display.display(XSDIR[[*columns]])



def processInput():
    description= "Preparing to list available ACE data"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--xsdir', type=pathlib.Path,
        default = pathlib.Path(os.environ['DATAPATH'], 'xsdir'),
        help="Path to xsdir file. Defaults to $DATAPATH/xsdir")

    parser.add_argument('-N', type=int,
        default=multiprocessing.cpu_count()-1,
        help="Number of parallel threads.")

    parser.add_argument('--dont-generate', action='store_true', default=False,
        help="Don't generate {}".format(xsdirName))

    return parser.parse_args()


def generateJSON(xsdirPath):
    """
    generateJSON will generate the JSON version of the XSDIR pandas DataFrame.
    It saves the fil
    """
    ddir = DataDirectory(xsdirPath)
    for index in ddir.XSDIR.index:
        ddir.extend(index)

    # with multiprocessing.Pool(args.N) as pool:
    #     pool.map(ddir.extend, ddir.XSDIR.index)

    with open(xsdirName, 'w') as jsonFile:
        json = ddir.XSDIR.to_json(orient='records', default_handler=str, indent=2)
        jsonFile.write(json)


if __name__ == "__main__":

    args = processInput()
    if not args.dont_generate:
        generateJSON(args.xsdir)

    XSDIR = loadXSDIR()
    display = DisplayData(XSDIR)
    display()

