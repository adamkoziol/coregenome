#!/usr/bin/env python
try:
    from accessoryFunctions import *
except ImportError:
    from SPAdesPipeline.OLCspades.accessoryFunctions import *
__author__ = 'adamkoziol'


class Core(object):

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """
        :param args: command line arguments
        :param pipelinecommit: pipeline commit or version
        :param startingtime: time the script was started
        :param scriptpath: home path of the script
        """
        import multiprocessing
        import coretyper
        import annotate
        # Initialise variables
        self.commit = str(pipelinecommit)
        self.start = startingtime
        self.homepath = scriptpath
        # Define variables based on supplied arguments
        self.args = args
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.sequencepath = os.path.join(args.sequencepath, '')
        assert os.path.isdir(self.sequencepath), u'Supplied sequence path is not a valid directory {0!r:s}' \
            .format(self.sequencepath)
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = int(args.numthreads if args.numthreads else multiprocessing.cpu_count())
        self.createdatabase = args.createdatabase
        self.dockerimage = args.dockerimage
        self.genus = args.genus
        self.species = args.species
        try:
            self.pipeline = args.pipeline
        except AttributeError:
            self.pipeline = False
        try:
            self.metadata = args.metadata
        except AttributeError:
            self.metadata = False

        if self.pipeline:
            self.profilelocation = args.profilelocation
            self.coregenelocation = args.coregenelocation
        if self.createdatabase:
            self.databasesequencepath = os.path.join(args.databasesequences, '')
            assert self.databasesequencepath, u'Please supply a valid path to the location of the files to be used ' \
                                              u'in creating the core genome database and profile'
            assert self.genus, u'Please provide the genus of the strains being used to construct the core database'
            assert self.species, u'Please provide the species of the strains being used to construct the core database'
            assert os.path.isdir(self.databasesequencepath), \
                u'Supplied path to the location of the files to be used in creating the core genome database and ' \
                u'profile is invalid {0!r:s}'.format(self.databasesequencepath)
        else:
            if args.annotateonly:
                self.databasesequencepath = self.sequencepath
                annotateonly = annotate.Annotate(self)
                annotateonly.annotatethreads()
                quit()
        # Run the analyses
        annotate.Annotate(self)
        # Run the core typing module
        coretyper.CoreTyper(self)


if __name__ == '__main__':
    import subprocess
    import time
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    # Parser for arguments
    parser = ArgumentParser(description='Calculate pan/core genome of a set of genomes')
    parser.add_argument('path',
                        help='Input directory. REQUIRED')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of assemblies to process. REQUIRED')
    parser.add_argument('-d', '--dockerimage',
                        required=True,
                        help='The name of the docker image to use. REQUIRED')
    parser.add_argument('-g', '--genus',
                        required=True,
                        help='Genus of the strains used to create the core database. REQUIRED')
    parser.add_argument('-n', '--numthreads',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-c', '--createdatabase',
                        action='store_true',
                        help='Calculate the core genes and profiles from a set of trusted sequences. Use these genes'
                             'and profiles to type your strain(s) of interest')
    parser.add_argument('-p', '--databasesequences',
                        help='If creating a database of core genes, and typing profiles, provide the path to the files'
                             'to use to create the database. REQUIRED if using -c flag')
    parser.add_argument('-S', '--species',
                        help='Species of the strains used to create the core database. REQUIRED if using -c flag')
    parser.add_argument('-a', '--annotateonly',
                        action='store_true',
                        help='Only perform annotations. Nothing else')

    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    start = time.time()

    # Run the script
    Core(arguments, commit, start, homepath)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m')


class PipelineInit(object):

    def __init__(self, inputobject, coremetadata):
        args = GenObject()
        pipelinecommit = inputobject.commit
        startingtime = inputobject.starttime
        scriptpath = inputobject.homepath
        args.path = inputobject.path
        args.sequencepath = inputobject.path
        args.dockerimage = '192.168.1.5:5000/coregenome'
        args.threads = inputobject.cpus
        args.genus = 'Escherichia'
        args.species = 'coli'
        args.numthreads = inputobject.cpus
        args.createdatabase = False
        args.pipeline = True
        args.metadata = coremetadata
        args.coregenelocation = os.path.join(inputobject.reffilepath, 'coregenome', 'Escherichia')
        args.profilelocation = args.coregenelocation
        Core(args, pipelinecommit, startingtime, scriptpath)

'''
-p
/nas0/bio_requests/8318/Genbank_Genomes/Listeria_Fasta/
-c

'''
