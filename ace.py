# vim: set fileencoding=utf-8 :

"""
This module provides abilities for parsing ace files which come out of NJOY.
"""

import os
import collections
import warnings
import copy
import textwrap

import numpy


class ace(object):
    """
    ace is an object which parses and stores the data from an ace file.  For
    every variable described in Appendix F of the MCNP manual, is duplicated
    here prefixed with '_'. The data is usually also stored elsewhere, but will
    at least be stored according to the documentation.

    Note: much of the indexing here has unecessary '+1' or '-1'. This was
    done intentionally to make the indices similar to what is done in the
    MCNP manual and---presumably---in the code.
    """

    def __init__(self, filename=None, headerOnly=False, zaid=None, xsdir=None,
                 start_line=1):
        """
        headerOnly: If True, then only the header (i.e., NXS, JXS arrays) are
            read and parsed.

        start_line: Start reading the file at this line.
        """

        assert filename or (zaid and xsdir), \
            "Must have either a filename or a zaid and an xsdir object."

        if filename:
            self.filename = filename
            self.start_line = start_line
        elif zaid and xsdir:
            try:
                xsdir_entry = xsdir.zaids[zaid]
            except KeyError:    # zaid is not a full or accurate zaid

                if isinstance(zaid, int):
                    xsdir_entry = xsdir[zaid][0]
                elif isinstance(zaid, str) or isinstance(zaid, unicode):
                    xsdir_entry = xsdir[zaid][0]

            # Get filename of data file
            parent, tail = os.path.split(xsdir.filename)
            self.filename = os.path.join(os.path.abspath(parent),
                    xsdir_entry.filename)

            # Where to start reading in the file
            self.start_line = xsdir_entry.address

        # Don't know what to do
        else: raise SyntaxError("I can't determine what data to process.")

        super(ace, self).__init__()
        self.headerOnly = headerOnly

        # Set aside space for storing cross sections
        self.xs = collections.OrderedDict()

        # Open file and move to starting line number
        self._file = open(self.filename, 'r')
        for i in range(self.start_line-1):
            line = self._file.readline()

        # Get header (first 12 lines)
        self._processHeader()

        if not self.headerOnly:
            # Read XSS array
            self.XSS = self._XSS = numpy.fromfile(self._file, dtype=float,
                sep=' ', count = self._NXS[0]) # Number of entries on XSS array

            # Process XSS array
            # self._processXSS()

        self._file.close()

    def _processHeader(self):
        """
        _processHeader is called to process the header
        """
        # Determine if we are using old- or new-style header
        line = self._file.readline().strip()
        words = line.split()

        if len(words) > 3:
            # Old-style header
            self.isNewStyle = False

            self._processOldStyleHeader(words)
        else:
            # New-style header
            self.isNewStyle = True

            self._processNewStyleHeader(words)

        header = []
        for i in range(10): header.append(self._file.readline().strip())

        # IZ, AW
        izaw = ' '.join(header[:4]).split()
        self._IZ = numpy.array( [ izaw[i] for i in range(len(izaw)) if not i%2],
                dtype=numpy.int)
        self._AW = numpy.array( [ izaw[i] for i in range(len(izaw)) if i%2],
                dtype=numpy.float)

        # NXS
        NXS = ' '.join(header[4:6])
        self._NXS = numpy.fromstring(NXS, dtype='i8', sep=' ')
        self.NXS = collections.OrderedDict()
        for i, n in enumerate(self._NXS): self.NXS[i+1] = n

        # JXS
        JXS = ' '.join(header[6:])
        self._JXS = numpy.fromstring(JXS, dtype='i8', sep=' ')
        self.JXS = collections.OrderedDict()
        for i, n in enumerate(self._JXS): self.JXS[i+1] = n

    def _processNewStyleHeader(self, firstWords):
        """
        _processNewStyleHeader will process the header according to the
        new-style.

        firstWords: The first line of the data table split by white space
        """
        # Process first line
        self.Version = firstWords[0]

        self._HZ = self.ZAID = self.full_zaid = firstWords[1]
        self.zaid, self.zaid_suffix = self.ZAID.split('.')

        self.Source = firstWords[2]

        # Process second line
        words = self._file.readline().strip().split()
        # Atomic weight ratio
        self._AW0 = self.atomic_weight_ratio = float(words[0])
        # Temperature
        self._TZ = self.temperature = float(words[1])
        # Processing date
        self._HD = self.process_date = words[2]
        # Number of comment lines
        self._NComments = int(words[3])

        # Read the comment lines
        comments = []
        for n in range(self._NComments):
            comments.append(self._file.readline().strip())

        self.CommentLines = '\n'.join(comments)

    def _processOldStyleHeader(self, firstWords):
        """
        _processOldStyleHeader will process the header according to the
        old-style.

        firstWords: The first line of the data table split by white space
        """
        # ZAID
        self._HZ = self.ZAID = firstWords[0]
        self.zaid, self.suffix = self.ZAID.split('.')
        self.full_zaid = self.ZAID
        self.zaid = self.zaid
        try:
            self.zaid = int(self.zaid)
            self.Z = int(self.zaid/1000)
            self.A = int(self.zaid-self.Z*1000)

        except ValueError:  # Probably S(a,B)
            self.Z = None
            self.A = None

        # Atomic weight ratio
        self._AW0 = self.atomic_weight_ratio = float(firstWords[1])
        # Temperature
        self._TZ = self.temperature = float(firstWords[2])
        # Processing date
        self._HD = self.process_date = firstWords[3]

        line = self._file.readline().rstrip()
        # Comment
        self._HK = self.comment = line[:70].rstrip()
        # Material ID
        self._HM = self.mat_ID = line[70:].strip()

    def _processNXS(self):
        """
        _processNXS will process the array self._NXS and set other class
        attributes according to what is found
        """
        raise NotImplementedError

    def _processJXS(self):
        """
        _processJXS will process the array self._JXS and set other class
        attributes according to what is found
        """
        raise NotImplementedError

    def _processXSS(self):
        """
        _processXSS will process the array self._XSS and set other class
        attributes according to what is found
        """
        raise NotImplementedError

    def __repr__(self): return self.filename

if __name__ == "__main__":
    print("\n\nI'm an ACE!\n")

    path = '/Users/jlconlin/Documents/Data/type1/endf71x/Fm/100255.716nc'

    ACE = ace(filename=path, headerOnly=True)
