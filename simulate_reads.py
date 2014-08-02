#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
#

import argparse


def main():

    # Parse args
    args = parse_arguments()




def parse_arguments():
    """ Load command line args.
    """

    parser = argparse.ArgumentParser(description='Simulate IGH MiSeq reads.')

    # Required args
    parser.add_argument('--BaseCallDir',
                        metavar='<dir>',
                        type=str,
                        required=True,
                        help='Directory containing base call intensitites')
    parser.add_argument('--RunParamXML',
                        metavar='<runParameters.xml>',
                        type=str,
                        required=True,
                        help='runParameters.xml file')
    parser.add_argument('--Lane',
                        metavar='<int>',
                        type=int,
                        required=True,
                        help='Lane number')
    parser.add_argument('--ReadStructure',
                        metavar='<str>',
                        type=str,
                        required=True,
                        help='Description of the logical structure of clusters in an Illumina Run, e.g. 151T8B8B')

    # Optional arguments
    parser.add_argument('--JarLoc',
                        metavar='<*.jar>',
                        type=str,
                        required=False,
                        default=os.path.join(root_dir, 'lib/picard-tools-1.115/IlluminaBasecallsToFastq.jar'),
                        help='Location of IlluminaBasecallsToFastq.jar (../lib/picard-tools-1.115/IlluminaBasecallsToFastq.jar)')
    parser.add_argument('--numCPU',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=1,
                        help='Number of CPUs to use. (1)')
    parser.add_argument('--readsPerTile',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=120000,
                        help='Max number of reads in RAM per tile, reduce if you have problems with memory. (120000)')
    parser.add_argument('--MaxInRam',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=500000,
                        help='Maximum number of records that are stored in the RAM. (500000)')
    parser.add_argument('--JavaRAM',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=2,
                        help='Amount of RAM allocated to the Java heap. (2)')


    return parser.parse_args()

if __name__ == '__main__':
    main()
