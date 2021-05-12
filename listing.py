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
        self.problems = []
        self.datapath = xsdirPath.parent

        self.metadataFunctions = {
            'c': self._fastNeutron,
            'nc': self._fastNeutron,
            't': self._thermalScattering,
            'h': self._proton,
            'p': self._photon,
            'm': self._multigroup,
            'm': self._dosimetry,
        }

        AWRs, self.XSDIR = xsdir.readXSDIR(xsdirPath)

    def _fastNeutron(self, index):
        """
        Add metadata from all fast neutron ACE files.

        index: Location in self.XSDIR of row where we are adding metadata
        """
        path = pathlib.Path(self.datapath, self.XSDIR.loc[index].path)

        address = self.XSDIR.loc[index].address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {'ZAID': self.XSDIR.loc[index].ZAID}
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

    def _thermalScattering(self, index):
        """
        Add metadata from thermal scattering ACE files.

        index: Location in self.XSDIR of row where we are adding metadata
        """
        path = pathlib.Path(self.datapath, self.XSDIR.loc[index].path)

        address = self.XSDIR.loc[index].address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {'ZAID': self.XSDIR.loc[index].ZAID}
        meta['NA'] = int(ACE.NXS[3] + 1)
        meta['NE'] = int(ACE.NXS[4])

        if ACE.NXS[7] == 2:
            meta['representation'] = 'continuous'
        else:
            meta['representation'] = 'discrete'

        return meta

    def _photon(self, index):
        """
        Add metadata from photon ACE files.

        index: Location in self.XSDIR of row where we are adding metadata
        """
        path = pathlib.Path(self.datapath, self.XSDIR.loc[index].path)

        address = self.XSDIR.loc[index].address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {'ZAID': self.XSDIR.loc[index].ZAID}

        meta['length'] = int(ACE.NXS[1])
        meta['NE'] = int(ACE.NXS[3])
        return meta

    def _multigroup(self, index):
        """
        Add metadata from multigroup ACE files.

        index: Location in self.XSDIR of row where we are adding metadata
        """
        path = pathlib.Path(self.datapath, self.XSDIR.loc[index].path)

        address = self.XSDIR.loc[index].address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {'ZAID': self.XSDIR.loc[index].ZAID}

        meta['length'] = int(ACE.NXS[1])
        return meta

    def _proton(self, index):
        """
        Add metadata from proton ACE files.

        index: Location in self.XSDIR of row where we are adding metadata
        """
        path = pathlib.Path(self.datapath, self.XSDIR.loc[index].path)

        address = self.XSDIR.loc[index].address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {'ZAID': self.XSDIR.loc[index].ZAID}

        meta['length'] = int(ACE.NXS[1])
        return meta

    def _dosimetry(self, index):
        """
        Add metadata from dosimetry ACE files.

        index: Location in self.XSDIR of row where we are adding metadata
        """
        path = pathlib.Path(self.datapath, self.XSDIR.loc[index].path)

        address = self.XSDIR.loc[index].address
        ACE = ace.ace(filename=path, headerOnly=False, start_line=address)

        meta = {'ZAID': self.XSDIR.loc[index].ZAID}

        meta['length'] = int(ACE.NXS[1])
        return meta

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
        print("{}\t{}".format(self.XSDIR.loc[index].ZAID,
                              self.XSDIR.loc[index].path))

        try:
            return self.metadataFunctions.get(lib_type, self._default)(index)
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
            None: self._default
        }

    def __call__(self, ZA=None, columns=[]):
        """
        """
        if isinstance(self.lib_type, list):
            lt = self.lib_type[0]
        else:
            lt = self.lib_type

        if not columns:
            columns = self.displayColumns[lt]

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
    # for index in ddir.XSDIR.index:
    #     ddir.extend(index)

    with multiprocessing.Pool(args.N) as pool:
        results = pool.map(ddir.extend, ddir.XSDIR.index)

    results = pd.DataFrame([r for r in results if r])

    ddir.XSDIR = pd.merge(ddir.XSDIR, results, on='ZAID')

    with open(xsdirName, 'w') as jsonFile:
        json = ddir.XSDIR.to_json(orient='records', default_handler=str, indent=2)
        jsonFile.write(json)

    if ddir.problems:
        print("There were problems reading data from these ZAIDs:")
        for ZAID in ddir.problems:
            print("\t{}".format(ZAID))


if __name__ == "__main__":

    args = processInput()
    if not args.dont_generate:
        generateJSON(args.xsdir)

    XSDIR = loadXSDIR()
    display = DisplayData(XSDIR)
    display()
