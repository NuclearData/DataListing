# vim: set fileencoding=utf-8 :

from __future__ import print_function, division


"""
This module provides abilities for parsing ace files which come out of NJOY.
"""

import os
import collections
import warnings
import copy
import textwrap

import numpy


# Laws describing the energy distribution of outgoing particles
class _Law_base(object):
    """
    _Law_base is the base object for all the outgoing particle distributions

    LDAT: Location of data in xss array
    LDIS: Relative location of data (JXS[x])
    xss: Array where data is located
    law: What law is being processed (may not be helpful in base case).
    """
    def __init__(self, LDAT, LDIS, xss, law, **kwargs):
        super(_Law_base, self).__init__()
        self.LDAT = LDAT
        self.LDIS = LDIS
        self.xss = xss
        self.law = law

        self._process()

        del self.xss

    def _process(self):
        warnings.warn(
            "No function has been defined for law: {}.".format(self.law))
#       raise NotImplementedError

    def __repr__(self): return "Law: {}".format(self.law)


class _tabular_distribution(_Law_base):
    """
    _tabular_distribution is a base class which provides common functionality
    for multiple Law classes. In particular, it provides functions for Laws 4,
    44, 61 and 67.

    The data that _tabular_distribution will extract is
        self.nr: Number of interpolation regions
        self.NBT: ENDF interpolation parameters
        self.INT: ENDF interpolation parameters
        self.ne: Number of energies at which distributions are tabulated
        self.Ein: Incident neutron energies
        self.loc: Locations of distributions
    """

    def _process(self):

        # Number of interpolation regions
        self.nr = _NR = self.xss[int((self.LDAT+(1-1))-1)]
        # ENDF interpolation parameters
        self.NBT = self.xss[(int(self.LDAT+(2-1)+0*_NR)-1):int((self.LDAT+(2-1)+1*_NR)-1)]
        self.INT = self.xss[int((self.LDAT+(2-1)+1*_NR)-1):int((self.LDAT+(2-1)+2*_NR)-1)]

        # Incident energies
        K = int(self.LDAT+(2-1)+2*_NR)
        self.ne = _NE = self.xss[K-1]
        K += 1
        self.Ein = self.xss[int((K+0*_NE)-1):int((K+1*_NE)-1)]
        # Locations of distributions (relative to JXS[11] or JXS[19])
        self.loc = _L = self.xss[int((K+1*_NE)-1):int((K+2*_NE)-1)]

        # Get outgoing distribution
        self._process_outgoing_distribution()

    def _process_outgoing_distribution(self): raise NotImplementedError

class _spectrum_distribution(_Law_base):
    """
    _spectrum_distribution is a base class which provides common functionality
    for multiple Law classes.  In particular, it provides functions for Laws 5,
    7, and 9.

    The data that _spectrum_distribution will extract is:
        self.nr: Number of interpolation regions
        self.NBT: ENDF interpolation parameters
        self.INT: ENDF interpolation parameters
        self.ne: Number of incident energies tabulated
        self.Ein: Incident energy table
        self.T: Tabulated function of incident energies

    """
    def _process(self):
        # Number of interpolation regions
        self.nr = _NR = int(self.xss[ (self.LDAT+(1-1))-1 ])
        # ENDF interpolation parameters
        self.NBT = self.xss[(self.LDAT+(2-1)+0*_NR)-1:(self.LDAT+(2-1)+1*_NR)-1]
        self.INT = self.xss[(self.LDAT+(2-1)+1*_NR)-1:(self.LDAT+(2-1)+2*_NR)-1]

        # Incident energies
        self.K = K = self.LDAT+(2-1)+2*_NR
        self.ne = _NE = int(self.xss[ K-1 ])
        K += 1
        self.Ein = self.xss[ (K+0*_NE)-1:(K+1*_NE)-1 ]   # Incident energy
        self.T   = self.xss[ (K+1*_NE)-1:(K+2*_NE)-1 ]   # Tabulated T's

        # Get extra stuff
        self._process_extras()

    def _process_extras(self): raise NotImplementedError

class Law1(_Law_base):
    """
    Tabular Equiprobable Energy Bins (From ENDF Law 1)
    Table F.14a
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 1
        super(Law1, self).__init__(LDAT, LDIS, xss, self.law)
#       print(" Law: {}".format(self.law))

    def _process(self):
        # Number of interpolation regions
        self.nr = _NR = self.xss[ (self.LDAT+(1-1))-1 ]
        # ENDF interpolation parameters
        self.NBT = self.xss[(self.LDAT+(2-1)+0*_NR)-1:(self.LDAT+(2-1)+1*_NR)-1]
        self.INT = self.xss[(self.LDAT+(2-1)+1*_NR)-1:(self.LDAT+(2-1)+2*_NR)-1]

        # Incident energies
        K = self.LDAT+(2-1)+2*_NR
        self.ne = _NE = self.xss[ K-1 ]
        K += 1
        self.Ein  = self.xss[ (K+0*_NE)-1:(K+1*_NE)-1 ]

        # Outgoing energies
        K = self.LDAT+(3-1)+2*_NR+_NE
        self.net = _NET = self.xss[ K-1 ] # Number of outgoing energies

        # Return array
        self.outgoing = numpy.zeros(
            int(_NE), dtype=[('Ein', 'f8'),('Eout','f8', int(_NET))])
        self.outgoing['Ein'] = self.Ein

        # Store outgoing energies
        K = self.LDAT+4+2*_NR+_NE-1
        for i, E in enumerate(self.Ein):
            mask = self.outgoing['Ein'] == E
            eout = self.xss[ (K+i*_NET)-1:(K+(i+1)*_NET)-1 ]
            self.outgoing['Eout'][mask] = eout

class Law3(_Law_base):
    """
    Level Scattering (From ENDF Law 3)
    Table F.14c

    The data stored here is simply a length-2 array, self.L. See the MCNP
    documentation for more information.
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 3
        super(Law3, self).__init__(LDAT, LDIS, xss, self.law)
#       print(" Law: {}".format(self.law))

    def _process(self): self.L = self.xss[int(self.LDAT-1):int((self.LDAT+2)-1)]

class Law4(_tabular_distribution):
    """
    Continuous Tabular Distribution (From ENDF Law 1)
    Table F.14d

    In addition to the data extracted from the base class, the following is
    extracted for each incident energy and stored in a namedtuple:
        np: Number of interpolation points
        Eout: Outgoing energies
        pdf:
        cdf:
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 4
        super(Law4, self).__init__(LDAT, LDIS, xss, self.law)
#       print(" Law: {}".format(self.law))

    def _process_outgoing_distribution(self):
        self.outgoing = collections.OrderedDict()
        K = int(self.LDAT + (3+2*self.nr+2*self.ne) - 1)

        out_data = collections.namedtuple('outgoing', 'np Eout pdf cdf')
        # Iterate through all the incident energies
        for e in self.Ein:
            # Interpolation scheme
            inttp = self.xss[int(K-1)]
            if inttp > 10:
                _ND = int(inttp)
                _INTT = inttp - 10*_ND
            else:
                _ND = None
                _INTT = inttp

            _NP   = self.xss[int((K+1)-1)]    # Number of points
            _Eout = self.xss[int((K+2+0*_NP)-1):int((K+2+1*_NP)-1)] # Outgoing Energies
            _PDF  = self.xss[int((K+2+1*_NP)-1):int((K+2+2*_NP)-1)] # PDF
            _CDF  = self.xss[int((K+2+2*_NP)-1):int((K+2+3*_NP)-1)] # CDF

            # Store outgoing data
            self.outgoing[e] = out_data(_INTT,  _Eout, _PDF, _CDF)

            # Update K for next iteration
            K += 2+3*_NP

class Law7(_spectrum_distribution):
    """
    Simple Maxwell Fission Spectrum (From ENDF–6 File 5 LF=7)
    Table F.14f

    In addition to the data provided by the _spectrum_distribution base class,
    Law7 also extracts
        self.U: Restriction energy
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 7
        super(Law7, self).__init__(LDAT, LDIS, xss, self.law)
#       print(" Law: {}".format(self.law))

    def _process_extras(self):

        K = self.K
        self.U = self.xss[ (K+2*self.ne)-1 ]

class Law9(_spectrum_distribution):
    """
    Evaporation Spectrum (From ENDF–6 File 5 LF=9)
    Table F.14g
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 9
        super(Law9, self).__init__(LDAT, LDIS, xss, self.law)
#       print(" Law: {}".format(self.law))

    def _process_extras(self):
        K = self.K
        self.U = self.xss[ (K+2*self.ne)-1 ]

class Law44(_tabular_distribution):
    """
    Kalbach-87 Formalism (LDAT, self.xssFrom ENDF File 6 Law 1, LANG=2)
    Table F.14k

    In addition to the data extracted from the base class, the following is
    extracted for each incident energy and stored in a namedtuple:
        np: Number of interpolation points
        Eout: Outgoing energies
        pdf:
        cdf:
        R: Precompound fraction
        A: Angular distribution slope value
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 44
        super(Law44, self).__init__(LDAT, LDIS, xss, self.law)
#       print(" Law: {}".format(self.law))

    def _process_outgoing_distribution(self):
        self.outgoing = collections.OrderedDict()
        K = self.LDAT + (3+2*self.nr+2*self.ne) -1

        out_data = collections.namedtuple('outgoing', 'np Eout pdf cdf R A')
        # Iterate through all the incident energies
        for e in self.Ein:
            # Interpolation scheme
            inttp = self.xss[int(K-1)]
            if inttp > 10:
                _ND = int(inttp)
                _INTT = inttp - 10*self._ND
            else:
                _ND = None
                _INTT = inttp

            _NP   = self.xss[int((K+1)-1)]    # Number of points
            _Eout = self.xss[int((K+2+0*_NP)-1):int((K+2+1*_NP)-1)] # Outgoing Energies
            _PDF  = self.xss[int((K+2+1*_NP)-1):int((K+2+2*_NP)-1)] # PDF
            _CDF  = self.xss[int((K+2+2*_NP)-1):int((K+2+3*_NP)-1)] # CDF
            _R    = self.xss[int((K+2+3*_NP)-1):int((K+2+4*_NP)-1)] # Precompound fraction r
            _A    = self.xss[int((K+2+4*_NP)-1):int((K+2+5*_NP)-1)] # Angular dist. slope

            # Store outgoing data
            self.outgoing[e] = out_data(_INTT,  _Eout, _PDF, _CDF, _R, _A)

            # Update K for next iteration
            K += 2+5*_NP

class Law61(_tabular_distribution):
    """
    Like LAW 44 but tabular angular distribution instead of Kalbach-87
    Table F.14l

    In addition to the data extracted from the base class, the following is
    extracted for each incident energy and stored in a namedtuple:
        np: Number of interpolation points
        Eout: Outgoing energies
        pdf:
        cdf:
        angle: Angular distribution (stored in a namedtuple)
            np: Number of points
            cos: Cosine scattering angular grid
            pdf:
            cdf:
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 61
        super(Law61, self).__init__(LDAT, LDIS, xss, self.law)
#       print(" Law: {}".format(self.law))

    def _process_outgoing_distribution(self):

        self.outgoing = collections.OrderedDict()

        out_data = collections.namedtuple('outgoing',
            field_names='interp np Eout pdf cdf, angle')

        # Iterate through al the incident energies
        K = int(self.LDAT + (3+2*self.nr+2*self.ne) - 1)
        for e in self.Ein:
            # Interpolation scheme
            inttp = self.xss[int(K-1)]
            if inttp > 10:
                _ND = int(inttp)
                _INTT = inttp - 10*self._ND
            else:
                _ND = None
                _INTT = inttp

            _NP   = self.xss[int((K+1)-1)]    # Number of points
            _Eout = self.xss[int((K+2+0*_NP)-1):int((K+2+1*_NP)-1)] # Outgoing Energies
            _PDF  = self.xss[int((K+2+1*_NP)-1):int((K+2+2*_NP)-1)] # PDF
            _CDF  = self.xss[int((K+2+2*_NP)-1):int((K+2+3*_NP)-1)] # CDF

            _LC = self.xss[int((K+2+3*_NP)-1):int((K+2+4*_NP)-1)]
            K += 2+4*_NP

            # Angular distribution for outgoing energies
            angle = collections.namedtuple('angle', 'interp np cos pdf cdf')
            angle_dist = []
            for loc in _LC:
                if loc > 0:
                    dist = collections.OrderedDict()       # Distribution for an outgoing energy

                    L = self.LDIS+abs(loc)-1
                    interp = _JJ  = self.xss[int((L+(1-1))-1)]
                    np     = _aNP = self.xss[int((L+(2-1))-1)]
                    cos = self.xss[int((L+(3-1)+0*_aNP)-1):int((L+(3-1)+1*_aNP)-1)]
                    pdf = self.xss[int((L+(3-1)+1*_aNP)-1):int((L+(3-1)+2*_aNP)-1)]
                    cdf = self.xss[int((L+(3-1)+2*_aNP)-1):int((L+(3-1)+3*_aNP)-1)]

                    dist = angle(interp, np, cos, pdf, cdf)

                # Isotropic outgoing
                elif loc == 0: dist = 'isotropic'

                # Save angular distribution
                angle_dist.append(dist)

            # Update K so it can find the next bit of data
            K = L+(3-1)+3*_aNP

            # Store angular data
            self.outgoing[e] = out_data(_INTT, _NP, _Eout, _PDF, _CDF,
                    angle_dist)

class Law66(_Law_base):
    """
    N-body phase space distribution (LDAT, self.xssFrom ENDF File 6 Law 6)
    Table F.14m

    The data stored here is simply a length-2 array, self.L where
        self.L[0]: NPSX---Number of bodies in the phase space.
        self.L[1]: A_p---Total mass ratio for the NPSX particles
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 66
        super(Law66, self).__init__(LDAT, LDIS, xss, self.law)
#       print(" Law: {}".format(self.law))

    def _process(self): self.L = self.xss[ self.LDAT-1:(self.LDAT+2)-1 ]

class Law67(_tabular_distribution):
    """
    Laboratory Angle–Energy Law (LDAT, self.xssFrom ENDF File 6 Law 7)
    Table F.14n

    In addition to the data extracted from the base class, the following is
    extracted for each incident energy (stored in self.outgoing), stored in a
    named tuple:
        interp: Interpolation scheme
        cos: Secondary cosines
        E_dist: Energy distribution of secondary particles (namedtuple)
            interp: Interpolation parameter
            E: Secondary energies
            pdf:
            cdf:
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 67
        super(Law67, self).__init__(LDAT, LDIS, xss, self.law)
#       print(" Law: {}".format(self.law))

    def _process_outgoing_distribution(self):

        self.outgoing = collections.OrderedDict()

        # Get outgoing distributions
        K = self.LDAT+3+2*self.nr+2*self.ne-1

        out_data = collections.namedtuple('outgoing', 'interp cos cos_E_dist ')
        sec_energies = collections.namedtuple('sec_energies', 'interp E pdf cdf')

        # Iterate through incident energies
        for e in self.Ein:

            _INTMU = self.xss[ (K+0)-1 ] # Interpolation scheme
            _NMU   = self.xss[ (K+1)-1 ] # Number of secondary cosines
            _XMU = self.xss[ (K+2+0*_NMU)-1:(K+2+1*_NMU)-1 ] # Secondary cosines
            _LMU = self.xss[ (K+2+1*_NMU)-1:(K+2+2*_NMU)-1 ] # Loc. of secondary cosines

            E_dist = collections.OrderedDict()  # Energy distribution
            J = K+2+2*_NMU
            for x in _XMU:

                _INTEP = self.xss[ (J+0)-1 ] # Inerpolation parameter
                _NPEP  = self.xss[ (J+1)-1 ] # Number of secondary energies
                _EP    = self.xss[ (J+2+0*_NPEP)-1:(J+2+1*_NPEP)-1 ] # Sec. ener.
                _PDF   = self.xss[ (J+2+1*_NPEP)-1:(J+2+2*_NPEP)-1 ]
                _CDF   = self.xss[ (J+2+2*_NPEP)-1:(J+2+3*_NPEP)-1 ]
                J += 2+3*_NPEP

                # Store data for each outgoing energy (for each outgoing cosine)
                SE = sec_energies(_INTEP, _EP, _PDF, _CDF)
                E_dist[x] = SE

            # Store data for each outgoing cosine
            self.outgoing[e] = out_data(_INTMU, _XMU, E_dist)

            # Update K for next iteration
            K = J

class Law2(_Law_base):
    """
    Discrete Photon Energy
    Table F.14b
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 2
        super(Law2, self).__init__(LDAT, LDIS, xss, self.law)
        print(" Law: {}".format(self.law))

#   def _process(self):
#       pass

class Law5(_Law_base):
    """
    General Evaporation Spectrum (LDAT, self.xssFrom ENDF–6 File 5 LF=5)
    Table F.14e
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 5
        super(Law5, self).__init__(LDAT, LDIS, xss, self.law)
        print(" Law: {}".format(self.law))

#   def _process(self):
#       pass

class Law11(_Law_base):
    """
    Energy Dependent Watt Spectrum (LDAT, self.xssFrom ENDF–6 File 5 LF=11)
    Table F.14h

    Look at 92235.12c for example.
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 11
        super(Law11, self).__init__(LDAT, LDIS, xss, self.law)

#   def _process(self):
#       pass

class Law22(_Law_base):
    """
    Tabular Linear Functions (LDAT, self.xssfrom UK Law 2)
    Table F.14i
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 22
        super(Law22, self).__init__(LDAT, LDIS, xss, self.law)
        print(" Law: {}".format(self.law))

#   def _process(self):
#       pass

class Law24(_Law_base):
    """
    (LDAT, self.xssFrom UK Law 6)
    Table F.14j
    """
    def __init__(self, LDAT, LDIS, xss, **kwargs):
        self.law = 24
        super(Law24, self).__init__(LDAT, LDIS, xss, self.law)
        print(" Law: {}".format(self.law))

#   def _process(self):
#       pass

def _noLaw(LDAT, xss, law, **kwargs):
    """
    This function will be executed when a function is expected to process a
    particular energy distribution but has not yet been written.  This
    function simply raises a warning.

    self.LDAT:
    self.xss: Array where all data is contained
    law: Law name or number that is expected.
    """
    warnings.warn("No function has been defined for law: {}.".format(law))
    return None

_laws = {
#   -1: _noLaw, # What to execute when no function defining law is found
    1: Law1,   # Tabular Equiprobable Energy Bins (From ENDF Law 1)
    2: Law2,   # Discrete Photon Energy
    3: Law3,   # Level Scattering (From ENDF Law 3)
    4: Law4,   # Continuous Tabular Distribution (From ENDF Law 1)
    5: Law5,   # General Evaporation Spectrum (From ENDF–6 File 5 LF=5)
    7: Law7,   # Simple Maxwell Fission Spectrum (From ENDF–6 File 5 LF=7)
    9: Law9,   # Evaporation Spectrum (From ENDF–6 File 5 LF=9)
    11:Law11,  # Energy Dependent Watt Spectrum (From ENDF–6 File 5 LF=11)
    22:Law22,  # Tabular Linear Functions (from UK Law 2)
    24:Law24,  # (From UK Law 6)
    44:Law44,  # Kalbach-87 Formalism (From ENDF File 6 Law 1, LANG=2)
    61:Law61,  # Like LAW 44 but tabular angular distribution instead of Kalbach-87
    66:Law66,  # N-body phase space distribution (From ENDF File 6 Law 6)
    67:Law67,  # Laboratory Angle–Energy Law (From ENDF File 6 Law 7)
}

def _common_spectrum_law(LDAT, xss, D, **kwargs):
    """
    _common_spectrum_law provides common functionality between different laws
    that extract data of a spectrum. It returns K, an index used to find the
    data.

    This function is used by _law5, _law7, and _law9.
    """
    D['nr'] = _NR = int(self.xss[ (self.LDAT+(1-1))-1 ])  # Number of interpolation regions
    # ENDF interpolation parameters
    D['NBT'] = self.xss[ (self.LDAT+(2-1)+0*_NR)-1:(self.LDAT+(2-1)+1*_NR)-1]
    D['INT'] = self.xss[ (self.LDAT+(2-1)+1*_NR)-1:(self.LDAT+(2-1)+2*_NR)-1]

    # Incident energies
    K = self.LDAT+(2-1)+2*_NR
    D['ne'] = _NE = int(self.xss[ K-1 ])
    K += 1
    D['Ein'] = self.xss[ (K+0*_NE)-1:(K+1*_NE)-1 ]   # Incident energy
    D['T']   = self.xss[ (K+1*_NE)-1:(K+2*_NE)-1 ]   # Tabulated T's
#   D['U']   = self.xss[ (K+2*_NE)-1 ]

    return K

def _common_law_in(LDAT, xss, D, **kwargs):
    """
    _common_la_inw will provide some common funcationality between the different
    laws for getting the information for the incident energies.
    """
    D['nr'] = _NR = self.xss[ (self.LDAT+(1-1))-1 ]  # Number of interpolation regions
    # ENDF interpolation parameters
    D['NBT'] = self.xss[ (self.LDAT+(2-1)+0*_NR)-1:(self.LDAT+(2-1)+1*_NR)-1]
    D['INT'] = self.xss[ (self.LDAT+(2-1)+1*_NR)-1:(self.LDAT+(2-1)+2*_NR)-1]

    # Incident energies
    K = self.LDAT+(2-1)+2*_NR
    D['ne'] = _NE = self.xss[ K-1 ]
    K += 1
    D['Ein'] = self.xss[ (K+0*_NE)-1:(K+1*_NE)-1 ]
    # Locations of distributions (relative to JXS[11] or JXS[19])
    D['loc'] = _L = self.xss[ (K+1*_NE)-1:(K+2*_NE)-1 ]

def _common_law_out(LDAT, xss, law, D, LDIS, **kwargs):
    """
    _common_la_inw will provide some common funcationality between the different
    laws for getting the information for the outgoing distributions.
    """
    # Data for the outgoing energies as a funtion of incident energies
    D['out'] = {}

    fields = 'interp E pdf cdf '
    if   law ==  4:  pass
    elif law == 44: fields += 'R A'
    elif law == 61: fields += 'angle_dist'
    out_data = collections.namedtuple('outgoing', fields)

    K = self.LDAT + (3+2*D['nr']+2*D['ne']) -1
    for e in D['Ein']:
        _INTT = self.xss[ K-1 ]                          # Interpolation scheme
#       if _INTT > 10: INTT = 10*_ND+_INTT
        _NP   = self.xss[ (K+1)-1 ]                      # Number of points
        _EOUT = self.xss[ (K+2+0*_NP)-1:(K+2+1*_NP)-1 ]  # Outgoing energies
        _PDF  = self.xss[ (K+2+1*_NP)-1:(K+2+2*_NP)-1 ]  # PDF
        _CDF  = self.xss[ (K+2+2*_NP)-1:(K+2+3*_NP)-1 ]  # CDF

        if law == 4:
            # Store outgoing data
            OD = out_data(_INTT, _EOUT, _PDF, _CDF)
            D['out'][e] = OD

            # Update K for next iteration
            K += 2+3*_NP

        elif law == 44:
            _R = self.xss[ (K+2+3*_NP)-1:(K+2+4*_NP)-1 ]  # Precompound fraction r
            _A = self.xss[ (K+2+4*_NP)-1:(K+2+5*_NP)-1 ]  # Angular dist. slope

            # Store outgoing data
            OD = out_data(_INTT,  _EOUT, _PDF, _CDF, _R, _A)
            D['out'][e] = OD

            # Update K for next iteration
            K += 2+5*_NP

        elif law == 61:
            _LC = self.xss[ (K+2+3*_NP)-1:(K+2+4*_NP)-1 ]
            K += 2+4*_NP

            # Angular distribution for outgoing energies
            angle_dist = []
            for loc in _LC:
                if loc > 0:
                    dist = {}       # Distribution for an outgoing energy

                    L = LDIS+abs(loc)-1
                    _JJ = self.xss[ (L+1)-1 ]
                    _aNP = self.xss[ (L+2)-1 ]
                    _CSOUT = self.xss[ (L+3+0*_aNP)-1:(L+3+1*_aNP)-1 ]

                # Isotropic outgoing
                elif loc == 0: dist = 'isotropic'

            # Store outgoing data
            OD = out_data(_INTT, _EOUT, _PDF, _CDF, angle_dist)
            D['out'][e] = OD

    return D

class XS(object):
    """
    XS is a class for storing and retrieving the data for a cross section.
    Perhaps the most important attribute of this class is the EXS array.  EXS is
    a numpy structured array with columns 'Energy' and 'XS' representing the
    energy and cross section values.  You can also access this data through
    with self.Energy and self.XS.
    """
    xs_dtype = numpy.dtype([('Energy', 'f8'), ("XS", 'f8')])

    def __init__(self, mt, Energy, xs, name=None):
        """
        mt: Reaction number for this cross section
        Energy: Energy bins for this cross section
        XS: Cross section values for this cross section
        name: Name for the cross section, i.e., 'total' (optional)
        """
        super(XS, self).__init__()

        self.mt = mt

        self.EXS = numpy.array([ (E, xs) for E, xs in zip(Energy, xs) ],
                dtype=self.xs_dtype)

#       self.Energy = Energy
#       self.XS = xs

        self.name = name

    def __getattr__(self, attr):
        if attr == "Energy":
            return self.EXS['Energy']
        elif attr == "XS":
            return self.EXS['XS']
        else:
            raise AttributeError(attr)

    def __repr__(self):
        s = "MT: {}".format(self.mt)
        if self.name:
            s += " name"

        return s

    def Sample(self, E):
        """
        Sample will return the cross section value for the given E.  If E is not
        in self.Energy, it will linearly interpolate between bin values. If E is
        outside the range if input energies, then it simply returns a 0.

        Note: THIS HAS NOT BEEN CAREFULLY TESTED

        E: Energy of incident neutron
        """
        if E in self.EXS['Energy']:
            return self.EXS['XS'][ self.EXS['Energy'] == E ]
        else:
            # Get bin edges
            index = self.EXS['Energy'].searchsorted(E)

            # I don't extrapolate
            if index == 0 or index >= len(self.EXS):
                return 0

            xs = self.EXS['XS'][(index-1):(index+1)]
            En = self.EXS['Energy'][(index-1):(index+1)]

            # Linearly interpolate
            xsR = (xs[1]-xs[0])/(En[1]-En[0])*(E - En[0]) + xs[0]

            return xsR

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
        # Process header arrays
        # self._processNXS()
        # self._processJXS()

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
#       words = self._file.readline().strip().split()

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

    def PlotXS(self, mt, ax, *args, **kwargs):
        """
        PlotXS will plot the cross section from the data file.  It returns the
        line object of the data that was plotted.

        mt: ENDF reaction number
        ax: Where to plot the data
        args: arguments passed to the plotting command
        kwargs: keyword arguments passed to the plotting command
        """

        if mt in self.xs:
            E, xs = self.xs[mt].Energy, self.xs[mt].XS
            line = ax.plot(E, xs, **kwargs)
            return line

        # Reaction doesn't exist
        else: raise(IndexError, "Reaction {} does not exist.".format(mt))

    def xsdirEntry(self, customize={}):
        """
        xsdirEntry will return an MCNP.xsdir.xsdir_entry object. Useful for
        creating and customizing the entry. The keys useful for customization
        are:
            'sza' or 'zaid': The SZA ID for the ZAID
            'suffix' or 'zaid_suffix': The suffix for the ZAID
            'atomic_weight_ratio': The atomic weight ratio
            'filename': The path to the file
            'access': Unused in Type 1 files
            'file_type': Type 1 is the only supported option
            'start_line': Line number to start with
            'table_length': How many entries in the XSS array
            'record_length': Unused in Type 1 files
            'num_entries': Unused in Type 1 files
            'temperature': Temperature
            'ptable': Wether there are probability tables or not
        Be careful when customizing the xsdirEntry. It may return something
        useless.
        """
        # Get original values
        custom_dict = copy.copy(self.__dict__)
        # Add some hard coded things which are not in the ACE file
        custom_dict['access'] = 0       # Unused in Type 1 files
        custom_dict['file_type'] = 1
#       custom_dict['address'] = self.start_line
        custom_dict['table_length'] = self.NXS[1]
        custom_dict['record_length'] = 0    # Unused in Type 1 files
        custom_dict['num_entries'] = 0    # Unused in Type 1 files
        if self.JXS[23] > 0:
            custom_dict['ptable'] = 'ptable'
        else:
            custom_dict['ptable'] = ''

        # Customize values
        custom_dict.update(customize)

        custom_dict['full_zaid'] = '{zaid}.{zaid_suffix}'.format(
            zaid=custom_dict.get('sza', custom_dict.get('zaid')),
            zaid_suffix=custom_dict.get('suffix',
            custom_dict.get('zaid_suffix', '')))

        entry = '{full_zaid} {atomic_weight_ratio} {filename} {access} ' \
                '{file_type} {start_line} {table_length} {record_length} ' \
                '{num_entries} {temperature:.6E} {ptable}'.format(**custom_dict)

        # Make fix if entry is greater than 80 characters
        entrylines = textwrap.wrap(entry, width=75, subsequent_indent=' '*10,
            break_on_hyphens=False)
        fixedEntry = ' + \n'.join(entrylines)

        return fixedEntry

class ce(ace):
    """
    ce is an object which parses and stores the data from a continuous energy
    (or discrete reaction) neutron ace file from njoy.

    When processed, the object instance will have many pieces of data.  An
    attempt has been made to maintain to name the data according to the variable
    names in the MCNP, Appendix F manual, but preceded by '_'.

    ***WARNING*** Photon data parsing has not yet been implemented.

    nubar data: If the nubar data exists (self.JXS[2] != 0) namedtuple's are
    created containing the data for prompt, delayed, and/or total nubar. The
    data is saved as self.prompt_nubar, self.del_nubar, self.total_nubar
    respectively. The fields of the namedtuple's are:
        nr: Number of interpolation regions
        NBT: ENDF interpolation parameters
        INT: ENDF interpolation parameters
        ne: Number of energies
        E: Energies
        nubar: Values of nubar
        secondary_dist: Distribution of secondary particles (mostly for delayed
            neutrons). Same type of object as is stored in self.secondary_dist.
    For more information on these values, see the MCNP documentation in table
    F.5 in Appendix F.

    Reaction cross sections: The cross sections for the reactions (other than
    elastic scattering) are stored in a dictionary, self.xs. The keys to the
    dictionary are the reaction numbers which are stored in self.mt. The values
    of the dictionary are also dictionaries with keys 'energy' and 'xs'
    referring to the energy values and cross sections, respectively.

    For reactions that have secondary neutrons---including elastic
    scattering---the angular distribution of the outgoing particles can be found
    in the dictionary self.angular_dist. The keys of the dictionary are the
    reaction numbers---and 'elastic' for elastic scattering. The values are
    again dictionaries with keys 'energy' and 'dist' for the energy values and
    the angle distribution for each energy value. 'dist' is, once again, a
    dictionary with keys:
        'type': What kind of distribution
        'interp': What kind of interpolation between bin edges
        'cod': Cosine bins
        'pdf': PDF
        'cdf': CDF
    For more information, see table F.12 in the MCNP manual.

    The distribution of secondary particles are given for each reaction (mt #)
    in the dictionary self.secondary_dist.  Each secondary distribution is
    given as a dictionary (no surprise there) with the following keys:
        'location': Location of next law
        'law': Name of the law
        'IDAT': Location of the data for the law
        'NBT': ENDF interpolation parameters
        'INT': ENDF interpolation parameters
        'energy': Array of energies
        'probability': Array of probabilities
        'LDAT': Class with attributes with all the data.  All law classes are
            daughter classes of the _Law_base class.
    For more information see table F.14 and the subsequent tables
    describing the all the laws.
    """
    def _processNXS(self):
        """
        _processNXS will process the array self._NXS and set other class
        attributes according to what is found
        """
        self.NXS['length of XSS']                           = self.NXS[1]
        self.NXS['zaid']                    = self._ZA      = self.NXS[2]
        self.NXS['num energies']            = self._NES     = self.NXS[3]
        self.NXS['num reactions']           = self._NTR     = self.NXS[4]
        self.NXS['num reactions with secondary neutrons'] \
                                            = self._NR      = self.NXS[5]
        self.NXS['num photon production reactions'] \
                                            = self._NTRP    = self.NXS[6]
#       self.NXS = self.NXS[7]
        self.NXS['num delayed neutron precursor families'] \
                                            = self._NPCR    = self.NXS[8]
        if self.isNewStyle:
            self.NXS['metastable state']    = self.S        = self.NXS[9]
            self.NXS['Z']                   = self.Z        = self.NXS[10]
            self.NXS['A']                   = self.A        = self.NXS[11]
#       self.NXS = self.NXS[12]
#       self.NXS = self.NXS[13]
#       self.NXS = self.NXS[14]
        self.NXS['num PIKMT reactions']     = self._NT      = self.NXS[15]
        self.NXS['photon production']       =            bool(self.NXS[16])

    def _processJXS(self):
        """
        _processJXS will process the array self._JXS and set other class
        attributes according to what is found
        """
        self.JXS['energy table location']               = self._ESZ     = self.JXS[1]
        self.JXS['fission nu data location']            = self._NU      = self.JXS[2]
        self.JXS['mt array location']                   = self._MTR     = self.JXS[3]
        self.JXS['Q-value array location']              = self._LQR     = self.JXS[4]
        self.JXS['reaction type array location']        = self._TYR     = self.JXS[5]
        self.JXS['xs  table locators']                  = self._LSIG    = self.JXS[6]
        self.JXS['xs location']                         = self._SIG     = self.JXS[7]
        self.JXS['angular distribution locators']       = self._LAND    = self.JXS[8]
        self.JXS['angular distributions location']      = self._AND     = self.JXS[9]
        self.JXS['energy distribution table locators']  = self._LDLW    = self.JXS[10]
        self.JXS['energy distribution location']        = self._DLW     = self.JXS[11]
        self.JXS['photon production data location']     = self._GPD     = self.JXS[12]
        self.JXS['photon production mt array location'] = self._MTRP    = self.JXS[13]
        self.JXS['photon production xs table locators'] = self._LSIGP   = self.JXS[14]
        self.JXS['photon production xs location']       = self._SIGP    = self.JXS[15]
        self.JXS['photon production angular distribution locators']  \
                                                        = self._LANDP   = self.JXS[16]
        self.JXS['photon production angular distributions location']      \
                                                        = self._ANDP    = self.JXS[17]
        self.JXS['photon production energy distribution table locators'] \
                                                        = self._LDLWP   = self.JXS[18]
        self.JXS['photon production energy distributions table location'] \
                                                        = self._DLWP    = self.JXS[19]
        self.JXS['yield multipliers table location']    = self._YP      = self.JXS[20]
        self.JXS['total fission xs location']           = self._FIS     = self.JXS[21]
        self.JXS['total fission xs table end']          = self._END     = self.JXS[22]
        self.JXS['probability table location']          = self._LUNR    = self.JXS[23]
        self.JXS['delayed nubar data location']         = self._DNU     = self.JXS[24]
        self.JXS['basic delayed data location']         = self._BDD     = self.JXS[25]
        self.JXS['energy distribution table locators']  = self._DNEDL   = self.JXS[26]
        self.JXS['energy distributions location']       = self._DNED    = self.JXS[27]
#       self.JXS  = self._JXS[28]
#       self.JXS  = self._JXS[29]
#       self.JXS  = self._JXS[30]
#       self.JXS  = self._JXS[31]
#       self.JXS  = self._JXS[32]

    def _processXSS(self):
        """
        _processXSS will process the array self._XSS and set other class
        attributes according to what is found. It will call other "processing"
        functions to process/identify different parts of the XSS array. These
        other functions serve primarily to simplify and isolate the coding.
        """
        assert len(self._XSS) == self.NXS[1], "XSS array has improper length."

        # Dictionary of secondary particle distributions
        self.secondary_dist = collections.OrderedDict()

        self._ESZBlock()    # Always exists
        self._NUBlock()     # Exists if self.JXS[2] != 0
        self._MTRBlock()    # Exists if self.NXS[4] != 0
        self._LQRBlock()    # Exists if self.NXS[4] != 0
        self._TYRBlock()    # Exists if self.NXS[4] != 0
        self._LSIGBlock()   # Exists if self.NXS[4] != 0
        self._SIGBlock()    # Exists if self.NXS[4] != 0
        self._LANDBlock()   # Always exists
        self._ANDBlock()    # Always exists

        # Define secondary particle distributions
        if self.NXS[5] != 0:
            NMT = self.NXS[5]
            # Get locaters for secondary particle distributions
            self.LOCC = self._LDLWBlock(LED=self.JXS[10], NMT=NMT).astype(int)
            # Get secondary particle distributions
            self._DLWBlock(
                JED=self.JXS[11], NMT=NMT, LDIS=self.JXS[11])

    def _ESZBlock(self):
        """
        ESZ Block—contains the main energy grid for the table and the total,
        absorption, and elastic cross sections as well as the average heating
        numbers. The ESZ Block always exists. See Table F.4.
        """
        start_index = self.JXS[1]-1
        N = self.NXS[3]       # Number of values
        self.energies = self._E = self._XSS[ start_index:N ]

        # Total cross section
        self.xs_total = self._sigma_t = self._XSS[start_index+N:start_index+2*N]
        self.xs['total'] = XS(1, self.energies, self.xs_total, name='total')
        self.xs[1] = self.xs['total']

        # Total absorption cross section
        self.xs_abs_total = self._sigma_a = \
                self._XSS[ start_index+2*N:start_index+3*N ]
        self.xs['abs'] = XS(102, self.energies, self.xs_abs_total, name='abs')
        self.xs[102] = self.xs['abs']

        # Elastic cross section
        self.xs_el = self._sigma_el = \
                self._XSS[ start_index+3*N:start_index+4*N ]
        self.xs['elastic'] = XS(2, self.energies, self.xs_el, name='elastic')
        self.xs[2] = self.xs['elastic']

        # Average heating numbers
        self.ave_heating = self._H_ave = \
                self._XSS[ start_index+4*N: start_index+5*N ]
        self.xs['heating'] = XS(301, self.energies, self.ave_heating,
            name='heating')
        self.xs[301] = self.xs['heating']

    def _NUBlock(self):
        """
        NU Block—contains prompt, delayed and/or total nubar as a function of
        incident neutron energy.  The NU Block exists only for fissionable
        isotopes (that is, if JXS(2) ≠ 0). See Table F.5.

        """
        def _poly_nubar():
            """
            _poly_nubar will process the nubar data if it is in the polynomial
            function form. See Table F.5, form a) of the MCNP manual.
            """
            pass

        def _tab_nubar():
            """
            _tab_nubar will process the nubar data if it is in the tabular data
            form. See Table F.5, forma b) of the MCNP manual. It returns a
            namedtuple containing all the information.
            """

            # Interpolation regions
            nr = int(self._XSS[int((KNU+1)-1)])
            # Interpolation parameters
            _INT = self.XSS[int((KNU+(2-1)+0*nr)-1):int((KNU+(2-1)+1*nr)-1)]
            _NBT = self.XSS[int((KNU+(2-1)+1*nr)-1):int((KNU+(2-1)+2*nr)-1)]

            # Tabulated data
            K = int(KNU+2+2*nr)
            ne = int(self._XSS[K-1])
            # nubar energies
            K += 1
            E = self._XSS[ (K+0*ne)-1:(K+1*ne)-1 ]
            # nubar values
            nu = self._XSS[ (K+1*ne)-1:(K+2*ne)-1 ]

            fields = 'nr NBT INT ne E nubar secondary_dist '
            self.nubar_tup = collections.namedtuple('nubar', field_names=fields)
            nubar = self.nubar_tup(nr, _NBT, _INT, ne, E, nu, None)

            return nubar

        def _extract_nubar():
            """
            _extract_nubar will extract the nubar data from the data file. It
            first determines what the format is and then calls the appropriate
            function to handle that form of the data.  It returns the nubar
            array which is extracted from XSS.
            """
            # What form for NU array
            self._LNU = self._XSS[int(KNU-1)]

            # Extract nubar data
            if self._LNU == 1: return _poly_nubar()
            if self._LNU == 2: return _tab_nubar()

        def _delayed_nubar():
            """
            _delayed_nubar will extract the delayed nubar data.
            """
            KNU = self.JXS[24]
            dnb_dist = _extract_nubar()

            # Locations of delayed nubar distributions
            LDIS = self.JXS[27]
            locations = self._LDLWBlock(LED=self.JXS[26], NMT=self.NXS[8])

            dnb_secondary = []

            # Iterate through all the delayed neutron groups
            K = self.JXS[25]
            for loc in locations:
                D = {}
                D['decay_constant'] = self.XSS[int((K+0)-1)]
                D['nr']  = _NR      = self.XSS[int((K+1)-1)]
                D['NBT'] = _NBT     = self.XSS[int((K+2+0*_NR)-1):int((K+2+1*_NR)-1)]
                D['INT'] = _INT     = self.XSS[int((K+2+1*_NR)-1):int((K+2+2*_NR)-1)]
                D['ne']  = _NE      = self.XSS[int((K+2+2*_NR)-1)]
                D['energy'] = \
                    self.XSS[int((K+3+2*_NR+0*_NE)-1):int((K+3+2*_NR+1*_NE)-1)]
                D['probability'] = \
                    self.XSS[int((K+3+2*_NR+1*_NE)-1):int((K+3+2*_NR+2*_NE)-1)]

                # Update K for next iteration
                K += 3+2*_NR+2*_NE

                # Get distribution of outgoing neutrons
                D['secondary_dist'] = self._secondary_dist(loc, LDIS, 0, 0)

                dnb_secondary.append(D)

            # Prepare data in correct format to return
            tmp_data = list(dnb_dist)
            tmp_data[-1] = dnb_secondary
            dnb = self.nubar_tup(*tmp_data)

            return dnb


        # NU Block only exists of self.JXS[2] != 0
        if self.JXS[2] != 0:

            # Either prompt or total nubar
            if self._XSS[ (self.JXS[2])-1 ] > 0:
                KNU = self.JXS[2]
                self.prompt_nubar = self.total_nubar = _extract_nubar()

            # Both prompt and total nubar
            elif self._XSS[ (self.JXS[2])-1 ] < 0:

                # Prompt nubar
                KNU = self.JXS[2]+1
                self.prompt_nubar = _extract_nubar()

                # Total nubar
                KNU = self.JXS[2]+abs(self._XSS[ self.JXS[2]-1 ])+1
                self.total_nubar = _extract_nubar()

            # Delayed nubar
            if self.JXS[24] > 0: self.delayed_nubar = _delayed_nubar()
            else: self.delayed_nubar = None
        else:
            self._KNU = None
            self._LNU = None

    def _MTRBlock(self):
        """
        MTR Block---contains a list of ENDF/B MT numbers for all neutron
        reactions other than elastic scattering. The MTR Block exists for all
        isotopes that have reactions other than elastic scattering (that is,
        all isotopes with NXS(4) ≠ 0). See Table F.6.
        """
        if self.NXS[4] != 0:
            LMT = int(self.JXS[3])    # Where my reactions are listed
            NMT = int(self.NXS[4])    # Number of mt reactions
            self.mt = (self._XSS[ (LMT-1):(LMT-1)+NMT ]).astype(int)
        else:
            self.mt = None

    def _LQRBlock(self):
        """
        LQR Block---contains a list of kinematic Q-values for all neutron
        reactions other than elastic scattering. The LQR Block exists if
        NXS(4) ≠ 0. See Table F.7.
        """
        if self.NXS[4] !=0:
            self.Q_values = \
                self._XSS[ (self.JXS[4]-1):(self.JXS[4]-1)+self.NXS[4] ]
        else: self.Q_values = None

    def _TYRBlock(self):
        """
        TYR Block---contains information about the type of reaction for all
        neutron reactions other than elastic scattering. Information for each
        reaction includes the number of secondary neutrons and whether
        secondary neutron angular distributions are in the laboratory or
        center-of-mass system. The TYR Block exists if NXS(4) ≠ 0. See Table
        F.8.
        """
        if self.NXS[4] !=0:
            self.rx_type = self.TYR = \
            self._XSS[int(self.JXS[5]-1):
                      int((self.JXS[5]-1)+self.NXS[4])].astype(int)
        else: self.rx_type = self.TYR = None


    def _LSIGBlock(self):
        """
        LSIG Block---contains a list of cross-section locators for all neutron
        reactions other than elastic scattering. The LSIG Block exists if
        NXS(4) ≠ 0. See Table F.9.
        """
        if self.NXS[4] !=0:
            LXS = int(self.JXS[6])
            NMT = int(self.NXS[4])
            self._LOCA = self._XSS[ (LXS-1):(LXS-1)+NMT ]
        else:
            self._LOCA = None

    def _SIGBlock(self):
        """
        SIG Block---contains cross sections for all reactions other than elastic
        scattering. The SIG Block exists if NXS(4) ≠ 0. See Table F.10.
        """
        if self.NXS[4] != 0:

            for mt, loc in zip(self.mt, self._LOCA):
                # Energy array
                E_index = self.XSS[int((self.JXS[7]+loc-1)-1)]
                len_E = self.XSS[int((self.JXS[7]+loc-1))]
                Energy = self.XSS[ int(E_index-1):int((E_index-1)+len_E) ]

                xs_start = self.JXS[7]+loc+1-1
                xs = self.XSS[ int(xs_start):int(xs_start+len_E) ]

                self.xs[int(mt)] = XS(mt, Energy, xs)

    def _LANDBlock(self):
        """
        LAND Block---contains a list of angular-distribution locators for all
        reactions producing secondary neutrons. The LAND Block always exists.
        See Table F.11.
        """
        self._LOCB = self.XSS[ (self.JXS[8]-1):
                               (self.JXS[8]-1)+(self.NXS[5]) ]

    def _ANDBlock(self):
        """
        AND Block---contains angular distributions for all reactions producing
        secondary neutrons.  The AND Block always exists. See Table F.12.
        """
        def _angular_data(locb):
            """
            _angular_data will extract the angular data, put it in a dictionary
            and return it.

            locb: Value of self._LOCB array associated with desired array.
                Relative (to self.JXS[9]) location of angular distribution
                array.
            """
            if locb > 0:
                angle_dist = {}

                start = jxs9+locb-1
                # Energies
                NE = self.XSS[ int(start-1) ]    # Number of energies
                E = self.XSS[ int(start):int((start)+NE) ]
                angle_dist['energy'] = E
                # Angular tables locations
                LC = self.XSS[ int(start+NE):int(start+2*NE) ]

                # Find angular distributions
                distributions = collections.OrderedDict()
                for e, lc in zip(E, LC):
                    dist = collections.OrderedDict()
                    distributions[e] = dist

                    if lc > 0:  # 32-equiprobable bin distribution
                        dist['type'] = 'equiprobable'
                        dist['interp'] = 1
                        dist['cos'] = self.XSS[ (jxs9+lc-1)-1 ]
                        dist['pdf'] = numpy.zeros(32, dtype='f8')+1/32
                        dist['cdf'] = dist['pdf'].cumsum()

                    elif lc == 0: # Isotropic
                        dist['type'] = 'isotropic'
                        dist['interp'] = 1
                        dist['cos'] = numpy.array([-1., 1.])
                        dist['pdf'] = numpy.array([1.0])
                        dist['cdf'] = numpy.array([1.0])

                    elif lc < 0:    # Tabular
                        dist['type'] = 'tabular'
                        LDAT = jxs9+abs(lc)-1
                        # Interpolation format
                        interp = self.XSS[int(LDAT-1)]
                        dist['interp'] = int(interp)

                        np = self.XSS[int((LDAT+1)-1)]
                        K = LDAT+2
                        dist['cos'] = self.XSS[int((K+0*np)-1):int((K+1*np)-1)]
                        dist['pdf'] = self.XSS[int((K+1*np)-1):int((K+2*np)-1)]
                        dist['cdf'] = self.XSS[int((K+2*np)-1):int((K+3*np)-1)]

                angle_dist['dist'] = distributions

                return angle_dist
            elif locb == 0: return 'isotropic'
            elif locb == -1: return -1
            else: return None

        self.angular_dist = collections.OrderedDict()
        jxs9 = self.JXS[9]  # For convenience

        # Elastic scattering
        if self.NXS[5] > 0:
            self.angular_dist['elastic'] = _angular_data(self._LOCB[0])
        else: self.angular_dist = None

        # Other reactions
        for i, loc in enumerate(self._LOCB[1:]):
            # Get reaction number
            mt = self.mt[i]
            self.angular_dist[mt] = _angular_data(loc)

    def _LDLWBlock(self, LED, NMT):
        """
        LDLW Block---contains a list of energy distribution locators for all
        reactions producing secondary neutrons except for elastic scattering.
        The LDLW Block exists if NXS(5) ≠ 0. See Table F.13.

        _LDLWBlock returns an array of locations.
        """
        return self.XSS[ (LED-1):(LED-1)+(NMT) ]

    def _DLWBlock(self, JED, NMT, LDIS):
        """
        DLW Block---contains energy distributions for all reactions producing
        secondary neutrons except for elastic scattering. The DLW Block exists
        if NXS(5) ≠ 0. See Table F.14.
        """

        # Iterate through
        for i, loc in enumerate(self.LOCC):
            mt = self.mt[i]
            TYR = self.TYR[i]

            if mt > 100:
                warnings.warn(
                    "I don't know how to handle reaction #: {}".format(mt))
            else:
                self.secondary_dist[mt] = self._secondary_dist(loc, LDIS, JED,
                                                               TYR)

    def _secondary_dist(self, locc, LDIS, JED, TYR):
        """
        _secondary_dist will extract the secondary distribution data for
        reactions with secondary neutrons (excluding elastic scattering).
        It will put the relevant information in a dictionary and return it.
        """
        secondary_dist = {}

        secondary_dist['LNW']  = _LNW = \
            self.XSS[int((LDIS+locc-1)-1)]       # Location of law

        # Raise a warning to a user that this isn't complete
        # When implementing this, an example is 8017.70c
        if _LNW != 0:
            warnings.warn("Haven't parsed data for more than one law.")

        # Law name based upon law location
        secondary_dist['law'] =_LAW = int(self.XSS[int((LDIS+locc)-1)])

        # Location of data
        secondary_dist['IDAT'] = _IDAT = self.XSS[int((LDIS+locc+1)-1)]
        _NR =self.XSS[int((LDIS+locc+2)-1)]       # Num interps.

        # ENDF interpolation parameters
        K = LDIS+locc+3
        secondary_dist['NBT'] = self.XSS[int((K+0*_NR)-1):int((K+1*_NR)-1)]
        secondary_dist['INT'] = self.XSS[int((K+1*_NR)-1):int((K+2*_NR)-1)]

        # Energy and probabilities
        K = int(LDIS+locc+3+2*_NR)
        ne = self.XSS[K-1]        # Number of energies
        K += 1
        secondary_dist['energy'] = self.XSS[int((K+0*ne)-1):int((K+1*ne)-1)]
        secondary_dist['probability'] = self.XSS[int((K+1*ne)-1):int((K+2*ne)-1)]

#           print('\t_LAW: {}'.format(_LAW))
        ldat = LDIS+_IDAT-1
        relative_loc = LDIS
        secondary_dist['LDAT'] = _laws.get(_LAW, _Law_base) \
            (ldat, LDIS, self.XSS, law=_LAW, A=self.A)
            # Arguments to law parsing function

        if TYR < 0:
            KY = int(JED + abs(TYR) - 101)
            secondary_dist['yield'] = self._energyDependentNeutronYield(KY)
        else:
            secondary_dist['yield'] = abs(TYR)

        return secondary_dist

    def _energyDependentNeutronYield(self, KY):
        """
        _energyDependentNeutronYield will extract the energy-dependent neutron
        yields for secondary neutrons. This is called when TYR < 0.
        """
        nYield = {}

        NR = int(self.XSS[KY-1])
        nYield['NBT'] = self.XSS[int((KY+1)-1):int((KY+1+NR)-1)]
        nYield['INT'] = self.XSS[int((KY+1+NR)-1):int((KY+1+2*NR)-1)]

        NE = int(self.XSS[(KY+1+2*NR)-1])
        nYield['energy'] = self.XSS[int((KY+2+2*NR)-1):
                                    int((KY+2+2*NR+NE)-1)]
        nYield['yields'] = self.XSS[int((KY+2+2*NR+NE)-1):
                                    int((KY+2+2*NR+2*NE)-1)]

        return nYield


class SaB(ace):
    """
    SaB is an object which parses and stores the data from an S(α,β) ace file
    from njoy.

    When processed, the object instance will have many pieces of data.  An
    attempt has been made to maintain to name the data according to the variable
    names in the MCNP, Appendix F manual, but preceded by '_'.

    When processing the ITXE block, a nested array is created called self.ITXE.
    Each element of the array contains the following information

    Ein: Incoming energy
    out: Outgoing information
        Eout 1->1 Outgoing energy
        mu 1->1: Array of equally-likely discrete cosines
        Eout 1->2
        mu 1->2
        ...
        Eout 1->N
        mu 1->N
        Eout 2->1
        mu 2-> 1
        ...
        Eout M->N
        mu M->N
    For specific information, see the documentation for S(α,β) data files.
    """

    def _processNXS(self):
        """
        _processNXS will process the array self._NXS and set other class
        attributes according to what is found
        """
        self.NXS['inelastic mode']                  = self._IDPNI = self.NXS[2]
        self.NXS['inelastic dimension parameter']   = self._NIL   = self.NXS[3]
        self.NXS['inel_num_exit_energy']            = self._NIEB  = self.NXS[4]
        self.NXS['elastic mode']                    = self._IDPNC = self.NXS[5]
        self.NXS['ellastic dimension parameter']    = self._NCL   = self.NXS[6]
        self.NXS['secondary energy mode']           = self._IFENG = self.NXS[7]
#       self.NXS      = self.NXS[8]
#       self.NXS      = self.NXS[9]
#       self.NXS      = self.NXS[10]
#       self.NXS      = self.NXS[11]
#       self.NXS      = self.NXS[12]
#       self.NXS      = self.NXS[13]
#       self.NXS      = self.NXS[14]
        self.NXS['num PIKMT reactions']             = self._NT    = self.NXS[15]
        self.NXS['photon production']               =         bool(self.NXS[16])

    def _processJXS(self):
        """
        _processJXS will process the array self._JXS and set other class
        attributes according to what is found
        """

        self.JXS['inelastic energy location']   = self._ITIE    = self.JXS[1]
        self.JXS['inelastic xs location']       = self._ITIX    = self.JXS[2]
        self.JXS['inelastic angle distribution location'] \
                                                = self._ITXE    = self.JXS[3]
        self.JXS['el energy location']          = self._ITCE    = self.JXS[4]
        self.JXS['elastic xs location']         = self._ITCX    = self.JXS[5]
        self.JXS['elastic angle distribution location'] \
                                                = self._ITCA    = self.JXS[6]
#       self.JXS     = self.JXS[7]
#       self.JXS     = self.JXS[8]
#       self.JXS     = self.JXS[9]
#       self.JXS     = self.JXS[10]
#       self.JXS     = self.JXS[11]
#       self.JXS     = self.JXS[12]
#       self.JXS     = self.JXS[13]
#       self.JXS     = self.JXS[14]
#       self.JXS     = self.JXS[15]
#       self.JXS     = self.JXS[16]
#       self.JXS     = self.JXS[17]
#       self.JXS     = self.JXS[18]
#       self.JXS     = self.JXS[19]
#       self.JXS     = self.JXS[20]
#       self.JXS     = self.JXS[21]
#       self.JXS     = self.JXS[22]
#       self.JXS     = self.JXS[23]
#       self.JXS     = self.JXS[24]
#       self.JXS     = self.JXS[25]
#       self.JXS     = self.JXS[26]
#       self.JXS     = self.JXS[27]
#       self.JXS     = self.JXS[28]
#       self.JXS     = self.JXS[29]
#       self.JXS     = self.JXS[30]
#       self.JXS     = self.JXS[31]
#       self.JXS     = self.JXS[32]

    def _processXSS(self):
        """
        """
        # ITIE Block
        self._ITIEBlock()

        # ITCEBlock
        self._ITCEBlock()

        # ITCABlock
        self._ITCABlock()

        # ITXE Block
        weight_type = self.NXS[7]
        if   weight_type == 1: self._ITXE_constant()
        elif weight_type == 2: self._ITXE_tabulated()
        else:
            warnings.warn(
                "A function to process for weighting type {} has not yet "\
                "been created.".format(weight_type))

    def _ITIEBlock(self):
        """
        # ITIE Block: energy-dependent inelastic scattering cross sections
        # --- Always exists
        """

        # Number of inelastic energies
        K = self.JXS[1]
        self.nEin = self._NEin = int(self.XSS[ (K+0)-1 ])
        # Inelastic energies
        self.Ein  = self._Ein  =     self.XSS[ (K+1)-1:(K+1+self._NEin)-1 ]
        # Inelastic cross sections
        xs = self._sigma_in = self.XSS[
                (K+1+1*self._NEin)-1:(K+1+2*self._NEin)-1 ]
        self.xs['inelastic'] = XS(4, self.Ein, xs, name='inelastic')
        self.xs[4] = self.xs['inelastic']

    def _ITCEBlock(self):
        """
        ITCE Block: energy-dependent elastic scattering cross sections
        --- exists if JXS(4) != 0
        """

        if self.JXS[4]:
            # Number of elastic energies
            self.nEel = self._NEel = int(self.XSS[ (self.JXS[4]) -1 ])
            # Elastic energies
            self.Eel = self._Eel = self.XSS[
                (self.JXS[4]+1+0*self._NEel)-1: (self.JXS[4]+1+1*self._NEel)-1 ]
            # Elastic cross sections
            self._sigma_el = self.XSS[
                (self.JXS[4]+1+1*self._NEel)-1: (self.JXS[4]+1+2*self._NEel)-1 ]

            # "un-normalize" if needed
            if self.NXS[5] == 4:
                Energies = [ self.Eel[0] ]
                xs = [ self._sigma_el[0]/self.Eel[0]*1E-2 ]

                for i in range(1, len(self.Eel)):
                    Energies.append(self.Eel[i-1])
                    Energies.append(self.Eel[i])

                    xs.append( self._sigma_el[i-1]/self.Eel[i-1] )
                    xs.append( self._sigma_el[i-1]/self.Eel[i] )
            else:
                Energies = self.Eel
                xs = self._sigma_el


            self.xs['elastic'] = XS(2, Energies, xs, name='elastic')
            self.xs[2] = self.xs['elastic']

        # Make sure things get set to None appropriately
        else:
            self.nEel = self._NEel = None
            self.Eel = self._Eel = None
            self._sigma_el = None

    def _ITCABlock(self):
        """
        ITCA Block: angular distributions for elastic scattering
        --- exists if JXS(4) != 0 and NXS(6) != -1
        """
        if self.JXS[4] and (self.NXS[6] != -1):
            num_mu = self.NXS[6]+1

            itca_dtype = numpy.dtype([('Eel', 'f8'), ('mu', 'f8', num_mu)])
            self.ITCA = self._ITCA = \
                numpy.zeros(self.nEel, dtype=itca_dtype)
            self.ITCA['Eel'] = self.Eel
            # Iterate through all incoming energies
            K = self.JXS[6]
            for i in range(self.nEel):
                self.ITCA['mu'][i] = self.XSS[
                    (K+i*num_mu)-1:(K+(i+1)*num_mu)-1 ]

        # Make sure things get set to None appropriately
        else: self.ITCA = self._ITCA = None

    def _ITXE_constant(self):
        """
        _ITXE_constant will process the array self._XSS if
        self.NXS[7] == 1.

        """
        # ITXE Block: coupled energy/angle distributions for inelastic scat.
        # --- Always exists
        num_mu = self.NXS[3] + 1  # Number of discrete cosines
        num_Eout = self.NXS[4]    # Number of energy out 'groups'

        itxe_dtype = numpy.dtype([('Ein','f8'),
            ('out', [ ('Eout', 'f8'), ('mu', 'f8', num_mu) ], num_Eout) ] )
        self.ITXE = self._ITXE = numpy.zeros(self._NEin, dtype=itxe_dtype)
        self.ITXE['Ein'] = self.Ein

        start = self.JXS[3]
        for i in range(1, self._NEin+1):    # Simulate indexing start at 1.
            # Get part of array we are currently interested in
            itxe_sub = self.ITXE[i-1]

            # Iterate over all the outgoing 'groups'
            for j in range(1, self.NXS[4]+1):
                Eout_index = j-1
                # Only worry about one part of the array
                itxe_Eout_sub = itxe_sub[1]

                # Outgoing energy
                index = start + (Eout_index)*(self.NXS[3]+2)
                Eout = itxe_Eout_sub['Eout'][j-1] = self.XSS[ index-1 ]

                # Equally-likely discrete cosines
                index += 1
                mu = self.XSS[(index)-1:(index+num_mu)-1]
                itxe_Eout_sub['mu'][Eout_index] = mu

            start = index+num_mu

    def _ITXE_tabulated(self):
        """
        _ITXE_tabulated will process the array self._XSS if
        self.NXS[7] == 2.
        """
        # ITXE Block: coupled energy/angle distributions for inelastic scat.
        # --- Always exists
        K = self.JXS[3] # Location of start of ITXE block?
        # Inelastic locations
        self.inel_out_loc = _locn = self.XSS[(K+0*self.nEin)-1:(K+1*self.nEin)-1]
        # Number of outgoing bins for every incident energy
        self.inel_out_bin = _nbin = self.XSS[(K+1*self.nEin)-1:(K+2*self.nEin)-1]

        num_mu = self.NXS[3]-1  # Number of outgoing cosine bins
        itxe_dtype = numpy.dtype([
            ('Eout','f8'), ('pdf','f8'), ('cdf','f8'), ('mu','f8',num_mu ) ])

        # Array length for each bin
        n = self.NXS[3]+2

        # Collect all outgoing information
        self.ITXE  = self._ITXE = collections.OrderedDict()
        for i, E in enumerate(self.Ein):
            l = self.inel_out_loc[i]

            # Get outgoing data for every incident energy
            N = int(self.inel_out_bin[i])
            self.ITXE[E] = numpy.zeros(N, dtype=itxe_dtype)
            for j in range(N):
                K = l+n*j
                self.ITXE[E][j]['Eout'] = self.XSS[ (K+0) ]
                self.ITXE[E][j]['pdf']  = self.XSS[ (K+1) ]
                self.ITXE[E][j]['cdf']  = self.XSS[ (K+2) ]
                self.ITXE[E]['mu'][j]   = self.XSS[ (K+3):(K+3)+(n-3) ]


class proton(ace):
    """
    proton is a object that parses and stores the data from a proton ACE file
    from NJOY. The extension is '.h'.
    """
    def _processNXS(self):
        """
        _processNXS will process the array self._NXS and set other class
        attributes according to what is found
        """
        self.NXS['length of XSS']                           = self.NXS[1]

    def _processJXS(self):
        """
        _processJXS will process the array self._JXS and set other class
        attributes according to what is found
        """
        pass


class photon(ace):
    """
    photon is a object that parses and stores the data from a photon ACE file
    from NJOY. The extension is '.p'.
    """
    def _processNXS(self):
        """
        _processNXS will process the array self._NXS and set other class
        attributes according to what is found
        """
        self.NXS['length of XSS']                           = self.NXS[1]

    def _processJXS(self):
        """
        _processJXS will process the array self._JXS and set other class
        attributes according to what is found
        """
        pass

if __name__ == "__main__":
    print("\n\nI'm an ACE!\n")

    path = '/Users/jlconlin/Documents/Data/type1/endf71x/Fm/100255.716nc'

    ACE = ce(filename=path, headerOnly=True)
