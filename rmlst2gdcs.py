#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject, SetupLogging
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import SeqIO
from argparse import ArgumentParser
from threading import Lock, Thread
from queue import Queue
import multiprocessing
import logging
import shutil
import numpy
import time
import csv
import os
__author__ = 'adamkoziol'


class GDCS(object):

    def runner(self):
        """
        Run the methods in the correct order
        """
        # Extract the alleles from the database
        self.alleleparser()
        # Get the all the alleles into a single file
        if self.complete:
            self.allelecomplement()
        # Create .fasta files of the the allele sequences
        self.alleleretriever()
        # Align the alleles
        self.allelealigner()
        # Find probes
        self.probefinder()
        # Choose the best probes for each gene
        self.probes()

    def alleleparser(self):
        """
        Parse a .csv file of rMLST alleles, and find all alleles for each gene for each organism of interest present in
        the file
        """
        logging.info('Parsing alleles')
        # Initialise each organism of interest as a sub-dictionary
        for organism in self.organisms:
            self.alleledict[organism] = dict()
            # Add an Enterobacteriaceae-specific entry
            if organism == 'Escherichia' or organism == 'Salmonella' or organism == 'Enterobacter':
                self.alleledict['Enterobacteriaceae'] = dict()
        # Get all the gene names into a list
        with open(self.rmlstfile, 'r') as rmlst:
            # Grab the header from the file
            header = rmlst.readline().rstrip()
            # Find all the gene names in the header
            self.genes = [gene for gene in header.split(',') if gene.startswith('BACT')]
        # Further prepare the dictionary to store a set of alleles for each gene
        for organism in self.organisms:
            for gene in self.genes:
                # Based on in silico testing, certain rMLST genes have alleles that pass the criteria for creating
                # GDCS probes, but fail when reference mapping reads to the probes
                include = True
                for excludegenus, excludegene in self.excludedict.items():
                    if organism == excludegenus and gene == excludegene:
                        include = False
                if include:
                    self.alleledict[organism][gene] = set()
                    if organism == 'Escherichia' or organism == 'Salmonella' or organism == 'Enterobacter':
                        self.alleledict['Enterobacteriaceae'][gene] = set()
        # Read the csv file into memory as a dictionary
        rmlstdict = csv.DictReader(open(self.rmlstfile))
        # Iterate through all the entries in the dictionary
        for row in rmlstdict:
            # Discard entries that do not match organisms of interest
            if row['Genus'] in self.organisms:
                # Find the allele for each gene
                for gene in self.genes:
                    # Alleles that are not present ('N') don't have a sequence, so they are ignored
                    if row[gene] != 'N':
                        # Remove the 'actual' allele from the sample profile e.g. for 10 692 (N), the (N) is removed.
                        # Additionally, the 10 692 are split into 10 and 692
                        try:
                            allele = row[gene].split(' (')[0]
                        except IndexError:
                            allele = row[gene]
                        # Split on the space between two alleles (10 692) and add both to the set
                        for alleles in allele.split(' '):
                            try:
                                # Add the integer of the allele to the set
                                self.alleledict[row['Genus']][gene].add(int(alleles))
                            # Skips adding the allele if it is an 'N' or if the organism/gene pair is one of the
                            # excluded sets
                            except (ValueError, KeyError):
                                pass
                            # Add all the Enterobacteriaceae to a combined entry
                            if row['Genus'] == 'Escherichia' or row['Genus'] == 'Salmonella' \
                                    or row['Genus'] == 'Enterobacter':
                                try:
                                    self.alleledict['Enterobacteriaceae'][gene].add(int(alleles))
                                except ValueError:
                                    pass

    def allelecomplement(self):
        """
        Retrieve the required alleles from a file of all alleles, and create complete allele files
        """
        logging.info('Retrieving complete list of alleles')
        # Read the csv file into memory as a dictionary
        completedict = csv.DictReader(open(self.rmlstfile))
        # Iterate through all the entries in the dictionary
        for row in completedict:
            # Grab all the alleles found for every organism in the list
            for gene in row:
                if gene.startswith('BACT'):
                    # Remove the 'actual' allele from the sample profile e.g. for 10 692 (N), the (N) is removed
                    # Additionally, the 10 692 are split into 10 and 692
                    allele = row[gene].split(' (')[0]
                    # Split on the space between two alleles (10 692) and add both to the set
                    for alleles in allele.split(' '):
                        try:
                            # Add the integer of the allele to the set
                            self.completedict[gene].add(int(alleles))
                        except KeyError:
                            self.completedict[gene] = set()
                            try:
                                self.completedict[gene].add(int(alleles))
                            # Skips adding the allele if it is an 'N'
                            except ValueError:
                                pass
                        except ValueError:
                            pass

        # Index all the records in the allele file
        logging.info('Loading complete rMLST records')
        recorddict = SeqIO.index(self.allelefile, 'fasta')
        logging.info('Creating complete allele output files')
        # Create the organism-specific files of alleles
        outpath = os.path.join(self.path, 'outputalleles', 'complete', '')
        # Delete and recreate the output path - as the files are appended to each time, they will be too large if
        # this script is run more than once
        try:
            shutil.rmtree(outpath)
        except OSError:
            pass
        make_path(outpath)
        combined = os.path.join(outpath, 'gdcs_alleles.fasta')
        allelefilelist = list()
        with open(combined, 'w') as combined:
            for gene, alleles in sorted(self.completedict.items()):
                # Open the file to append
                allelefiles = os.path.join(outpath, '{}.tfa'.format(gene))
                allelefilelist.append(allelefiles)
                with open(allelefiles, 'a') as allelefile:
                    # Write each allele record to the file
                    for allele in sorted(alleles):
                        # Skip adding alleles that are no longer in the database
                        try:
                            SeqIO.write(recorddict['{}_{}'.format(gene, allele)], allelefile, 'fasta')
                            SeqIO.write(recorddict['{}_{}'.format(gene, allele)], combined, 'fasta')
                        except KeyError:
                            pass

    def alleleretriever(self):
        """
        Retrieve the required alleles from a file of all alleles, and create organism-specific allele files
        """
        logging.info('Retrieving alleles')
        # Index all the records in the allele file
        logging.info('Loading rMLST records')
        recorddict = SeqIO.index(self.allelefile, 'fasta')
        logging.info('Creating allele output files')
        # Create the organism-specific files of alleles
        for organism in sorted(self.alleledict):
            # Make an object to store information for each strain
            metadata = MetadataObject()
            metadata.organism = organism
            metadata.path = self.path
            metadata.outpath = os.path.join(self.path, 'outputalleles', organism, '')
            # Delete and recreate the output path - as the files are appended to each time, they will be too large if
            # this script is run more than once
            try:
                shutil.rmtree(metadata.outpath)
            except OSError:
                pass
            make_path(metadata.outpath)
            metadata.combined = os.path.join(metadata.outpath, 'gdcs_alleles.fasta')
            metadata.allelefiles = list()
            with open(metadata.combined, 'w') as combined:
                for gene, alleles in sorted(self.alleledict[organism].items()):
                    # Open the file to append
                    allelefiles = os.path.join(metadata.outpath, '{}.tfa'.format(gene))
                    metadata.allelefiles.append(allelefiles)
                    with open(allelefiles, 'a') as allelefile:
                        # Write each allele record to the file
                        for allele in sorted(alleles):
                            # Skip adding alleles that are no longer in the database
                            try:
                                SeqIO.write(recorddict['{}_{}'.format(gene, allele)], allelefile, 'fasta')
                                SeqIO.write(recorddict['{}_{}'.format(gene, allele)], combined, 'fasta')
                            except KeyError:
                                pass
            # Add the populated metadata to the list
            self.samples.append(metadata)

    def allelealigner(self):
        """
        Perform a multiple sequence alignment of the allele sequences
        """

        logging.info('Aligning alleles')
        # Create the threads for the analysis
        for _ in range(self.cpus):
            threads = Thread(target=self.alignthreads, args=())
            threads.setDaemon(True)
            threads.start()
        for sample in self.samples:
            sample.alignpath = os.path.join(self.path, 'alignedalleles', sample.organism)
            make_path(sample.alignpath)
            # Create a list to store objects
            sample.alignedalleles = list()
            for outputfile in sample.allelefiles:
                aligned = os.path.join(sample.alignpath, os.path.basename(outputfile))
                sample.alignedalleles.append(aligned)
                # Create the command line call
                clustalomega = ClustalOmegaCommandline(infile=outputfile,
                                                       outfile=aligned,
                                                       threads=4,
                                                       auto=True)
                sample.clustalomega = str(clustalomega)
                self.queue.put((sample, clustalomega, outputfile, aligned))
        self.queue.join()

    def alignthreads(self):
        while True:
            sample, clustalomega, outputfile, aligned = self.queue.get()
            if not os.path.isfile(aligned):
                # Perform the alignments
                # noinspection PyBroadException
                try:
                    clustalomega()
                # Files with a single sequence cannot be aligned. Copy the original file over to the aligned folder
                except Exception:
                    shutil.copyfile(outputfile, aligned)
            self.queue.task_done()

    def probefinder(self):
        """
        Find the longest probe sequences
        """
        logging.info('Finding and filtering probe sequences')
        for sample in self.samples:
            # A list to store the metadata object for each alignment
            sample.gene = list()
            for align in sample.alignedalleles:
                # Create an object to store all the information for each alignment file
                metadata = GenObject()
                metadata.name = os.path.splitext(os.path.basename(align))[0]
                metadata.alignmentfile = align
                # Create an alignment object from the alignment file
                try:
                    metadata.alignment = AlignIO.read(align, 'fasta')
                except ValueError:
                    # If a ValueError: Sequences must all be the same length is raised, pad the shorter sequences
                    # to be the length of the longest sequence
                    # https://stackoverflow.com/questions/32833230/biopython-alignio-valueerror-says-strings-must-be-same-length
                    records = SeqIO.parse(align, 'fasta')
                    # Make a copy, otherwise our generator is exhausted after calculating maxlen
                    records = list(records)
                    # Calculate the length of the longest sequence
                    maxlen = max(len(record.seq) for record in records)
                    # Pad sequences so that they all have the same length
                    for record in records:
                        if len(record.seq) != maxlen:
                            sequence = str(record.seq).ljust(maxlen, '.')
                            record.seq = Seq(sequence)
                    assert all(len(record.seq) == maxlen for record in records)
                    # Write to file and do alignment
                    metadata.alignmentfile = '{}_padded.tfa'.format(os.path.splitext(align)[0])
                    with open(metadata.alignmentfile, 'w') as padded:
                        SeqIO.write(records, padded, 'fasta')
                    # Align the padded sequences
                    metadata.alignment = AlignIO.read(metadata.alignmentfile, 'fasta')

                metadata.summaryalign = AlignInfo.SummaryInfo(metadata.alignment)
                # The dumb consensus is a very simple consensus sequence calculated from the alignment. Default
                # parameters of threshold=.7, and ambiguous='X' are used
                consensus = metadata.summaryalign.dumb_consensus()
                metadata.consensus = str(consensus)
                # The position-specific scoring matrix (PSSM) stores the frequency of each based observed at each
                # location along the entire consensus sequence
                metadata.pssm = metadata.summaryalign.pos_specific_score_matrix(consensus)
                metadata.identity = list()
                # Find the prevalence of each base for every location along the sequence
                for line in metadata.pssm:
                    try:
                        bases = [line['A'], line['C'], line['G'], line['T'], line['-']]
                        # Calculate the frequency of the most common base - don't count gaps
                        metadata.identity.append(float('{:.2f}'.format(max(bases[:4]) / sum(bases) * 100)))
                    except KeyError:
                        bases = [line['A'], line['C'], line['G'], line['T']]
                        # Calculate the frequency of the most common base - don't count gaps
                        metadata.identity.append(float('{:.2f}'.format(max(bases) / sum(bases) * 100)))
                # List to store metadata objects
                metadata.windows = list()
                # Variable to store whether a suitable probe has been found for the current organism + gene pair.
                # As the probe sizes are evaluated in descending size, as soon as a probe has been discovered, the
                # search for more probes can stop, and subsequent probes will be smaller than the one(s) already found
                passing = False
                # Create sliding windows of size self.max - self.min from the list of identities for each column
                # of the alignment
                for i in reversed(range(self.min, self.max + 1)):
                    if not passing:
                        windowdata = MetadataObject()
                        windowdata.size = i
                        windowdata.max = 0
                        windowdata.sliding = list()
                        # Create a counter to store the starting location of the window in the sequence
                        n = 0
                        # Create sliding windows from the range of sizes for the list of identities
                        windows = self.window(metadata.identity, i)
                        # Go through each window from the collection of sliding windows to determine which window(s)
                        # has (have) the best results
                        for window in windows:
                            # Create another object to store all the data for the window
                            slidingdata = MetadataObject()
                            # Only consider the window if every position has a percent identity greater than the cutoff
                            if min(window) > self.cutoff:
                                # Populate the object with the necessary variables
                                slidingdata.location = '{}:{}'.format(n, n + i)
                                slidingdata.min = min(window)
                                slidingdata.mean = float('{:.2f}'.format(numpy.mean(window)))
                                slidingdata.sequence = str(consensus[n:n+i])
                                # Create attributes for evaluating windows. A greater/less windowdata.max/windowdata.min
                                #  means a better/less overall percent identity, respectively
                                windowdata.max = slidingdata.mean if slidingdata.mean >= windowdata.max \
                                    else windowdata.max
                                windowdata.min = slidingdata.mean if slidingdata.mean <= windowdata.max \
                                    else windowdata.min
                                # Add the object to the list of objects
                                windowdata.sliding.append(slidingdata)
                                passing = True
                            n += 1
                        # All the object to the list of objects
                        metadata.windows.append(windowdata)
                # All the object to the list of objects
                sample.gene.append(metadata)

    def probes(self):
        """
        Find the 'best' probes for each gene by evaluating the percent identity of the probe to the best recorded
        percent identity for that organism + gene pair
        """
        logging.info('Determining optimal probe sequences')
        for sample in self.samples:
            # Make a folder to store the probes
            sample.gdcsoutputpath = os.path.join(self.gdcsoutputpath, sample.organism)
            sample.gdcscombined = os.path.join(sample.gdcsoutputpath, '{}_gdcs_combined.fasta'.format(sample.organism))
            make_path(sample.gdcsoutputpath)
            with open(sample.gdcscombined, 'w') as combined:
                for gene in sample.gene:
                    # Open the file to append
                    gene.gdcsoutputfile = os.path.join(sample.gdcsoutputpath, '{}_gdcs.tfa'.format(gene.name))
                    with open(gene.gdcsoutputfile, 'w') as allelefile:
                        for window in gene.windows:
                            # Variable to record whether a probe has already been identified from this gene
                            passed = False
                            for sliding in window.sliding:
                                # Only consider the sequence if the sliding object has data, if the probe in question
                                # has a mean identity equal to the highest observed identity for that probe size, and
                                # if the mean identity is greater or equal than the lowest observed identity
                                if sliding.datastore and sliding.mean == window.max and sliding.mean >= window.min \
                                        and not passed:
                                    dnaseq = Seq(sliding.sequence, IUPAC.unambiguous_dna)
                                    # Create a sequence record using BioPython
                                    fasta = SeqRecord(dnaseq,
                                                      # Without this, the header will be improperly formatted
                                                      description='',
                                                      # Use the gene name as the header
                                                      id=gene.name)
                                    # Write each probe to the files
                                    SeqIO.write(fasta, allelefile, 'fasta')
                                    SeqIO.write(fasta, combined, 'fasta')
                                    passed = True

    @staticmethod
    def window(iterable, size):
        """
        https://coderwall.com/p/zvuvmg/sliding-window-in-python
        :param iterable: string from which sliding windows are to be created
        :param size: size of sliding window to create
        """
        i = iter(iterable)
        win = []
        for e in range(0, size):
            win.append(next(i))
        yield win
        for e in i:
            win = win[1:] + [e]
            yield win

    def __init__(self, args):
        """
        :param args: command line arguments
        """
        SetupLogging(debug=True)
        # Initialise variables
        self.start = args.start
        # Define variables based on supplied arguments
        self.path = os.path.join(args.path)
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.rmlstfile = os.path.join(self.path, args.file)
        self.organisms = args.organisms.split(',')
        self.allelefile = os.path.join(self.path, args.allelefile)
        self.gdcsoutputpath = os.path.join(self.path, 'gdcs')
        self.min = args.min
        self.max = args.max
        self.cutoff = args.cutoff
        self.complete = args.complete
        self.genes = list()
        self.alleledict = dict()
        self.completedict = dict()
        self.samples = list()
        self.cpus = multiprocessing.cpu_count()
        self.queue = Queue(maxsize=self.cpus)
        self.lock = Lock()
        self.excludedict = {'Listeria': 'BACT000014',
                            'Salmonella': 'BACT000062'}
        # Run the analyses
        self.runner()


if __name__ == '__main__':
    # Parser for arguments
    parser = ArgumentParser(description='For all organisms of interest, create .fasta files containing each allele'
                                        'found for every rMLST gene')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-f', '--file',
                        required=True,
                        help='Name of .csv file containing rMLST information. Must be within the supplied path')
    parser.add_argument('-o', '--organisms',
                        default='Bacillus,Campylobacter,Enterobacter,Escherichia,Listeria,Salmonella,Vibrio',
                        help='Comma-separated list of organisms of interest')
    parser.add_argument('-a', '--allelefile',
                        required=True,
                        help='File of combined rMLST alleles. This file must be within the supplied path')
    parser.add_argument('-m', '--min',
                        default=20,
                        help='Minimum size of probe to create')
    parser.add_argument('-M', '--max',
                        default=50,
                        help='Maximum size of probe to create')
    parser.add_argument('-c', '--cutoff',
                        default=70,
                        help='Cutoff percent identity of a nucleotide location to use')
    parser.add_argument('-C', '--complete',
                        action='store_true',
                        help='Optionally store all alleles found in the rMLST results file')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.pipeline = False
    # Define the start time
    arguments.start = time.time()

    # Run the script
    GDCS(arguments)

    # Print a bold, green exit statement
    logging.info('Analyses Complete!')

'''
/nas0/bio_requests/8318 -a rmlstcombinedalleles.fa -f rmlst.csv -C
'''
