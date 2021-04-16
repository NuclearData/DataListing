#! /usr/bin/env python

import collections
import os
import pathlib

import xsdir

def printEntries(entries):
    """
    printEntries will print the ZAIDs and temperatures for every entry
    """
    for entry in entries:
        if entry.temperature:
            print("\t{:12} T={:.4E} (MeV)".format(
                    entry.zaid, entry.temperature))
        else:
            print("\t{:12}".format(entry.zaid))

def library(args, XSDIR):
    """
    List all the different 'libraries' in the XSDIR file.
    """
    libraries = {entry.filename.parts[0] for entry in XSDIR.entries}
    print("Libraries ({}) in XSDIR:\n\t{}".format(len(libraries), 
                                                  "\n\t".join(libraries)))

def ZA(args, XSDIR):
    """
    List all the ZAIDS for a given ZA
    """
    ZAs = [entry for entry in XSDIR.entries 
                  if entry.za == args.za]

    if args.temperature:
        ZAs = [entry for entry in ZAs 
                      if entry.temperature == args.temperature]
    print("ZAIDS for ZA={}:")
    printEntries(ZAs)

def SaB(args, XSDIR):
    """
    List all the thermal scattering materials. If a material name is given, 
    list all the ZAIDs for that material.
    """
    print("Thermal Scattering (S(α,β)) materials:")
    mats = [entry for entry in XSDIR.entries
                  if entry.lib_type == 't']
    if args.material:
        mats = [entry for entry in mats
                       if entry.za == args.material]

        printEntries(mats)
    else:
        mats = {entry.za for entry in mats}
        for mat in mats:
            print("\t{}".format(mat))

if __name__ == "__main__":

    import argparse

    description= "Listing the available ACE data"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--datapath', type=pathlib.Path, 
        default=pathlib.Path(os.environ['DATAPATH'], 'xsdir'),
        help="Path to xsdir file. Defaults to $DATAPATH/xsdir")

    listers = parser.add_subparsers(dest='lister',
        help="List available parameters in XSDIR")

    ZALister = listers.add_parser('ZAs',
        help="List all entries for a given ZA")
    ZALister.set_defaults(func=ZA)
    ZALister.add_argument('za', type=int,
        help="ZA identifier")
    ZALister.add_argument('-t', '--temperature', type=float,
        help="Temperature filter")

    sabLister = listers.add_parser('materials',
        help="List thermal scattering materials")
    sabLister.set_defaults(func=SaB)
    sabLister.add_argument('--material',
        help="Material identifier")
    sabLister.add_argument('-t', '--temperature', type=float,
        help="Temperature filter")

    libLister = listers.add_parser('libraries',
        help="List all libraries")
    libLister.set_defaults(func=library)

    args = parser.parse_args()

    XSDIR = xsdir.xsdir(args.datapath)

    print()
    if args.lister:
        args.func(args, XSDIR)
    else:
        parser.print_help()

