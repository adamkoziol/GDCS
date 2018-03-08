#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import printtime, make_path
import os
__author__ = 'adamkoziol'


class Format(object):

    def runner(self):
        self.parse()
        self.reformat()

    def parse(self):
        """
        Read in the .csv file
        """
        printtime('Parsing rMLST output file', self.start)
        # Open the rmlst file, and iterate through each line
        with open(self.rmlstfile, 'r') as rmlstfile:
            self.header = rmlstfile.readline()
            for line in rmlstfile:
                # Split the line on commas
                data = line.split(',')
                # Only append the data to the list of relevant sequence types if the line does not start with 'Strain'
                # (e.g. the header), and if the strain name is present. Certain strains without exact matches will
                # return multiple best hits, but only the first hit is required (this line will have the strain name)
                if data[0] != 'Strain' and data[0]:
                    # Replace the 'NA' value for Genus with the supplied organism name
                    data[1] = self.organism
                    # Add the data to the list that will be used to create the reformatted file
                    self.sequencetypes.append(data)

    def reformat(self):
        """
        Create the reformatted output file
        """
        printtime('Reformatting results', self.start)
        with open(os.path.join(self.outputpath, '{}_reformatted.csv'.format(self.organism)), 'w') as reformatted:
            reformatted.write(self.header)
            for line in self.sequencetypes:
                reformatted.write(','.join(line))

    def __init__(self, args):
        """
        :param args: command line arguments
        """
        # Initialise variables
        self.start = args.start
        # Define variables based on supplied arguments
        self.path = os.path.join(args.path)
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.rmlstfile = os.path.join(self.path, args.file)
        self.organism = args.organism
        self.outputpath = os.path.join(self.path, 'reformatted')
        make_path(self.outputpath)
        self.header = str()
        self.sequencetypes = list()
        self.runner()


if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    import time
    # Parser for arguments
    parser = ArgumentParser(description='Reformats a rMLST output file to be consistent with rmlst2gdcs.py. This '
                                        'involves adding a genus, and removing any lines that do not contain a '
                                        'strain name')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-f', '--file',
                        required=True,
                        help='Name of .csv file containing rMLST information. Must be within the supplied path')
    parser.add_argument('-o', '--organism',
                        help='Organism name to use when populating the new, formatted .csv file')

    arguments = parser.parse_args()
    arguments.pipeline = False
    # Define the start time
    arguments.start = time.time()

    # Run the script
    Format(arguments)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - arguments.start) + '\033[0m')

'''
/nas0/bio_requests/8318/rawoutput -f MLST_2017.03.17.18.10.01.csv -o Enterobacter
'''