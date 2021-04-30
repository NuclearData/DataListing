#! /usr/bin/env python3
# vim: set fileencoding=utf-8 :

import collections
import os
import pathlib

import pandas as pd

import xsdir
import maps

def library_list(XSDIR):
    """
    Return a list of 'libraries' contained in the XSDIR
    """
    return {entry.parts[0] for entry in XSDIR['path']}

def library(args, XSDIR):
    """
    List all the different 'libraries' in the XSDIR file.
    """
    libraries = library_list(XSDIR)
    print("Libraries ({}) in XSDIR:".format(len(libraries)))

    return libraries

def ZAs_list(XSDIR, ZA=None, temperature=None, lib_type=[]):
    """
    Return a list of XSDIR entries for a given ZA and (optionally) temperature
    or lib_type.
    """
    if ZA:
        ZAs = XSDIR.query('ZA == @ZA')
    else:
        ZAs = XSDIR.ZA

    if temperature:
        ZAs = ZAs.query('temperature == @temperature')

    if lib_type:
        ZAs = ZAs[ZAs['lib_type'].isin(lib_type)]

    return ZAs

def printZA(args, XSDIR):
    """
    List all the ZAIDS for a given ZA
    """
    ZAs = ZAs_list(XSDIR, args.za, args.temperature, args.lib_type)
    print("ZAIDS for ZA={}:", args.za)
    return ZAs

def SaB_list(XSDIR, material, temperature=None):
    """
    Return a list of XSDIR entries for a given material and (optionally)
    temperature
    """
    mats = XSDIR.query('lib_type == "t"')
    if material:
        mats = mats.query('ZA == @material')
    else:
        mats = mats.ZA

    return mats

def printSaB(args, XSDIR):
    """
    List all the thermal scattering materials. If a material name is given, 
    list all the ZAIDs for that material.
    """
    print("Thermal Scattering (S(α,β)) materials:")
    mats = SaB_list(XSDIR, args.material, args.temperature)
    return mats

def parseXSDIR(datapath=None):
    """
    parseXSDIR will return a xsdir.xsdir object. If datapath isn't given, it
    uses the environment variable $DATAPATH
    """
    if not datapath:
        datapath = pathlib.Path(os.environ['DATAPATH'], 'xsdir')

    AWRs, entries = xsdir.readXSDIR(datapath)
    return entries

if __name__ == "__main__":

    import argparse

    description= "Listing the available ACE data"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--datapath', type=pathlib.Path, 
        default = pathlib.Path(os.environ['DATAPATH'], 'xsdir'),
        help="Path to xsdir file. Defaults to $DATAPATH/xsdir")

    listers = parser.add_subparsers(dest='lister',
        help="List available parameters in XSDIR")

    ZALister = listers.add_parser('ZAs',
        help="List all entries for a given ZA")
    ZALister.set_defaults(func=printZA)
    ZALister.add_argument('za', type=int,
        help="ZA identifier")
    ZALister.add_argument('-t', '--temperature', type=float,
        help="Temperature filter")
    ZALister.add_argument('-l', '--lib-type', type=str, default=None,
        choices = maps.DataTypes.keys(),
        help="Library type filter")

    sabLister = listers.add_parser('materials',
        help="List thermal scattering materials")
    sabLister.set_defaults(func=printSaB)
    sabLister.add_argument('--material',
        help="Material identifier")
    sabLister.add_argument('-t', '--temperature', type=float,
        help="Temperature filter")

    libLister = listers.add_parser('libraries',
        help="List all libraries")
    libLister.set_defaults(func=library)

    args = parser.parse_args()

    XSDIR = parseXSDIR(args.datapath)

    print()
    if args.lister:
        entries = list(args.func(args, XSDIR))
        print(entries)
    else:
        parser.print_help()
