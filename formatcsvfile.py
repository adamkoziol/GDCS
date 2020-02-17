#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import SetupLogging
from argparse import ArgumentParser
import logging
import os
__author__ = 'adamkoziol'


class Format(object):

    def runner(self):
        self.parse_confindr()
        self.parse_rmlst()
        self.reformat()

    def parse_confindr(self):
        """
        Read in the ConFindr report file. Ensure that the sequences used to create GDCS alleles are free from
        contamination
        """
        logging.info('Parsing ConFindr report')
        count = 0
        with open(self.confindrfile, 'r') as confindrfile:
            for line in confindrfile:
                count += 1
                data = line.rstrip().split(',')
                if data[1] in self.organisms:
                    if len(data) == 5:
                        if data[4] == 'Clean':
                            self.passing_strains.append(data[0])
                    elif len(data) == 6:
                        if data[3] == 'Clean':
                            self.passing_strains.append(data[0])

    def parse_rmlst(self):
        """
        Read in the rMLST .csv file
        """
        logging.info('Parsing rMLST output file')
        # Open the rmlst file, and iterate through each line
        with open(self.rmlstfile, 'r') as rmlstfile:
            self.header = rmlstfile.readline()
            for line in rmlstfile:
                # Split the line on commas
                data = line.split(',')
                # Only append the data to the list of relevant sequence types if the line does not start with 'Strain'
                # (e.g. the header), and if the strain name is present in the list of uncontaminated strains.
                # Certain strains without exact matches will return multiple best hits, but only the first hit is
                # required (this line will have the strain name)
                if data[0] != 'Strain' and data[0] not in self.checked_strains:
                    self.checked_strains.add(data[0])
                    if data[0] in self.passing_strains:
                        # Ensure that the genus is one of the target genera
                        if data[1] in self.organisms:
                            # Add the data to the list that will be used to create the reformatted file
                            self.strains.append(data)

    def reformat(self):
        """
        Create the reformatted output file
        """
        logging.info('Reformatting results')
        with open(os.path.join(self.path, 'combined_rmlst_reformatted.csv'), 'w') as reformatted:
            reformatted.write(self.header)
            for line in sorted(self.strains):
                reformatted.write(','.join(line))

    def __init__(self, args):
        """
        :param args: command line arguments
        """
        # Define variables based on supplied arguments
        self.path = os.path.abspath(os.path.join(args.path))
        assert os.path.isdir(self.path), 'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.reportpath = os.path.abspath(os.path.join(args.reportpath))
        assert os.path.isdir(self.reportpath), 'Could not locate the supplied report directory {0!r:s}'\
            .format(self.reportpath)
        self.rmlstfile = os.path.join(self.reportpath, 'rmlst.csv')
        assert os.path.isfile(self.rmlstfile), 'rMLST.csv not present in supplied report directory {0!r:s}'\
            .format(self.reportpath)
        self.confindrfile = os.path.join(self.reportpath, 'confindr_report.csv')
        assert os.path.isfile(self.confindrfile), 'confindr_report.csv not present in supplied report directory ' \
                                                  '{0!r:s}'.format(self.reportpath)
        self.organisms = [organism for organism in args.organisms.split(',')]
        self.header = str()
        self.strains = list()
        self.passing_strains = list()
        self.checked_strains = set()


if __name__ == '__main__':
    # Parser for arguments
    parser = ArgumentParser(description='Reformats a rMLST output file to be consistent with rmlst2gdcs.py. This '
                                        'involves adding a genus, and removing any lines that do not contain a '
                                        'strain name')
    parser.add_argument('-p', '--path',
                        required=True,
                        help='Specify working directory')
    parser.add_argument('-r', '--reportpath',
                        required=True,
                        help='Path of folder containing the rmlst.csv and confindr_report.csv files. These files'
                             'contain combined reports created by the reportaggregator package for all MiSeq runs')
    parser.add_argument('-o', '--organisms',
                        default='Bacillus,Campylobacter,Enterobacter,Escherichia,Listeria,Salmonella,Vibrio',
                        help='Comma-separated list of organisms for which the GDCS scheme is to be created. Default'
                             'is: Bacillus,Campylobacter,Enterobacter,Escherichia,Listeria,Salmonella,Vibrio')
    parser.add_argument('-v', '--verbose',
                        default=False,
                        action='store_true',
                        help='Enable debugging messages')
    arguments = parser.parse_args()
    SetupLogging(arguments.verbose)
    # Run the script
    format_csv = Format(arguments)
    format_csv.runner()
    # Print an exit statement
    logging.info('Analyses complete!')
