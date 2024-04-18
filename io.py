import argparse
import os
import os.path as path
import re

from PLI.__init__ import __version__


def parse_args():
    parser = argparse.ArgumentParser(description="A program to generate interaction profile between target "
                                                 "protein and ligand")

    parser.add_argument('-r', action='store', dest='receptor_path', default=False,
                        help='Single protein structure file [ex: prot.pdb]')

    # Input file settings
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-l', action='store', dest='ligand_path', default=False,
                       help='Single ligand structure in .pdbqt format [ex: somewhere/lig.pdbqt]')

    group.add_argument('-d', action='store', dest='lig_dir', default=False,
                       help='Direction of target ligands structure [ex: somewhere/ligs/]')

    # Output file settings
    parser.add_argument('-o', action='store', dest='output_path', default="./",
                        help='Output directory or file [ex : somewhere/ligs/]')

    # Special criteria of analysis
    parser.add_argument('--segment', nargs='*', action='store', dest='segid', default=False,
                        help='Assign which segment should be analysis [ex: SEG1 SEG2]')

    parser.add_argument('--active_site', nargs="*", dest='active_site', type=int, default=False,
                        help='Assign target active-site of protein (resid) [ex: 1 3 5 10]')

    parser.add_argument('--active_site_radius', dest='active_site_radius', type=int, default=20,
                        help='Determine the search space around to the target active-site [default: 10 A]')

    parser.add_argument('-v', '--version', action='version', version='Version ' + __version__)

    args = parser.parse_args()

    return args


def search_the_path_of_files(ligs_dir, lig_path):
    lig_files = []
    if ligs_dir:
        m = re.compile(r'.*.pdbqt')  # Selection criteria of compound name
        for root, dirs, files in os.walk(ligs_dir):
            for file in files:
                flag = m.match(file)
                if flag is not None:
                    full_path = path.join(root, file)
                    lig_files.append(full_path)  # Create ligand full-path list
    else:
        lig_files = [lig_path]

    return lig_files


def guess_output_filename(file_path, options):
    output_extension = ".csv"
    log_extension = ".log"
    
    name = path.basename(file_path).split('.')[0]
    logfile = path.join(options.output_path, name + log_extension)
    ofile = path.join(options.output_path, name + output_extension)
    

    return ofile, logfile
