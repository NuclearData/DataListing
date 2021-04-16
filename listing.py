import collections
import os
import pathlib

import xsdir

def libraryLister(args, XSDIR):
    """
    List all the different 'libraries' in the XSDIR file.
    """
    libraries = {entry.filename.parts[0] for entry in XSDIR.entries}
    print("Libraries in XSDIR:{}".format("\n\t".join(libraries)))


if __name__ == "__main__":

    import argparse

    description= "Listing the available ACE data"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--datapath', type=pathlib.Path, 
        default=pathlib.Path(os.environ['DATAPATH'], 'xsdir'),
        help="Path to xsdir file. Defaults to $DATAPATH/xsdir")

    subparsers = parser.add_subparsers()
    libLister = subparsers.add_parser('libraries',
        help="List all libraries")
    libLister.set_defaults(func=libraryLister)

    args = parser.parse_args()

    XSDIR = xsdir.xsdir(args.datapath)

    args.func(args, XSDIR)

