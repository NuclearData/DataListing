# vim: set fileencoding=utf-8 :

"""
This module will parse an xsdir file.
"""

import pathlib
import re
import os
import collections

_float_pattern = r'^|\s(\d*\.\d*)(?=$|\s)'
_int_pattern   = r'^|\s(\d+)(?=$|\s)'
_card_pattern  = re.compile(r'^\s{,3}(?=\d)', re.MULTILINE)
_isotope_pattern = re.compile(r'\s*(\d+)\s+(\d+\.\d+)')
_zaid_pattern = re.compile('''
        (?P<ZA>.+?)\.               # ZA
        (?P<lib>[0-9]{2,3})         # Extension
        (?P<type>[a-zA-Z]{1,2})     # Library type
        ''', re.VERBOSE)

class xsdir(collections.defaultdict):
    """
    xsdir will parse an xsdir file and organize all the iformation for "easy"
    retrieval.  The attribute 'zaids' is a dictionary with the full zaids (e.g.,
    '1001.70c') as the keys and the values are the xsdir_entry objects.

    The xsdir object itself is a dictionary with (not full) zaids (e.g., 1001)
    as the keys and the values are lists of every entry that matches the
    incomplete zaid; e.g.,
        [1001.70c, 1001.71c, 1001.72c, 1001.73c, 1001.74c, 1001.62c, 1001.66c,
        1001.60c, 1001.50c, 1001.42c, 1001.53c, 1001.24c, 1001.50d, 1001.30y,
        1001.50m, 1001.70h, 1001.24h]
    """
    def __init__(self, filename=pathlib.Path(os.environ['DATAPATH'], 'xsdir')):
        super(xsdir, self).__init__(list)

        if isinstance(filename, pathlib.PosixPath):
            self.filename = filename.absolute()
        else:
            self.filename = pathlib.Path(filename)
        self.atomic_weight_ratios = collections.OrderedDict()
        self.atomic_mass = collections.OrderedDict()

        self.zaids = {}
        self.entries = []

        self._parse()

    def _parse(self):
        with self.filename.open() as xsdirFile:
            line = '\n'.join(xsdirFile.readlines())
            xsdirFile.close()

            self.sections = line.split('directory')

            # Parse elements
            elements = _card_pattern.split(self.sections[0])
            for E in elements:
                Isotopes = _isotope_pattern.findall(E)
                if Isotopes:
                    self.atomic_mass[int(Isotopes[0][0])] = Isotopes[0][1]

                    for ZA, amr in Isotopes[1:]:
                        self.atomic_weight_ratios[int(ZA)] = amr

            # Parse directory
            lines = self.sections[1].split('\n')
            multipleLine = False
            for entry in lines:

                if entry:
                    if multipleLine:
                        Line += ' '.format(entry.lstrip())
                        multipleLine = False

                    elif entry.rstrip().endswith(' +'):
                        # Multiple line xsdir entry
                        Line = entry.rstrip()[:-1]
                        multipleLine = True
                        continue

                    else:
                        Line = entry
                        multipleLine = False


                    # Create xsdir_entry object
                    E = xsdir_entry(Line)

                    # Check to assure atomic weight ratio is provided
                    # try:
                    #     iZAID = int(E.za)
                    #     if not iZAID in self.atomic_weight_ratios:
                    #         if not iZAID in self.atomic_mass:
                    #             print("No atomic weight ratio for {}  {}".format(
                    #                 iZAID, E.atomic_weight_ratio))
                    # except ValueError:
                    #     pass

                    # Store entry in self
                    self.entries.append(E)
                    za = E.za
                    self[za].append(E)

                    # Store data individually instead of in a list
                    self.zaids[E.zaid] = E

class xsdir_entry(object):
    """
    xsdir_entry contains the information for each entry in an xsdir file.  For
    more information on what the data holds, look at Appendix F.II XSDIR---Data
    Directory File.
    """
    def __init__(self, entry):
        """
        entry: String---line from xsdir file
        """
        self.entry = entry
        words = entry.split()

        # Get zaid and library
        self.zaid = words[0]
        zaid_found = _zaid_pattern.search(self.zaid)
        if zaid_found:
            zaid_dict = zaid_found.groupdict()
            try:
                self.za = self.ZA = int(zaid_dict['ZA'])
            except ValueError:
                self.za = self.ZA = zaid_dict['ZA']

            self.lib = zaid_dict['lib']
            self.lib_type = zaid_dict['type']

        # Get remaining entries
        self.atomic_weight_ratio = float(words[1])
        self.filename = pathlib.Path(words[2])
        self.access = words[3]
        self.file_type = int(words[4])
        self.address = int(words[5])
        self.table_length = int(words[6])
        try: self.record_length = int(words[7])
        except IndexError: self.record_length=None
        try: self.num_entries = int(words[8])
        except IndexError: self.num_entries=None

        # Get temperature
        try: self.temperature = float(words[9])
        except IndexError: self.temperature=None

        # Probability table flag
        try: self.ptable = words[10]
        except IndexError: self.ptable = None

    def __repr__(self): return self.zaid

    def __str__(self):
        if self.temperature:
            entry = '{zaid} {atomic_weight_ratio} {filename} {access} ' \
                    '{file_type} {address} {table_length} {record_length} ' \
                    '{record_length} {temperature:.6E} {ptable}'.format(
                            **self.__dict__)
        else:
            entry = '{zaid} {atomic_weight_ratio} {filename} {access} ' \
                    '{file_type} {address} {table_length} {record_length} ' \
                    '{record_length} {temperature} {ptable}'.format(
                            **self.__dict__)
        return entry

if __name__ == "__main__":
    print("\nI'm parsing XSDIR file(s).\n")

    xs = xsdir()

