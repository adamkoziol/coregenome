#!/usr/bin/env python
# from SPAdesPipeline.OLCspades.accessoryFunctions import *
import SPAdesPipeline.OLCspades.metadataprinter as metadataprinter
from SPAdesPipeline.OLCspades.mMLST import *
from glob import glob
__author__ = 'adamkoziol'


class CoreTyper(object):

    def handler(self):
        """Run the required analyses"""
        printtime('Creating and populating objects', self.start)
        self.populate()
        printtime('Populating {} sequence profiles'.format(self.analysistype), self.start)
        self.profiler()
        # Annotate sequences with prokka
        self.annotatethreads()
        # Run the analyses
        self.cdsthreads()
        # Find core coding features
        self.cdssequencethreads()
        # Extract the sequence for each coding feature
        self.allelematchthreads()
        # Determine sequence types from the analyses
        printtime('Determining {} sequence types'.format(self.analysistype), self.start)
        self.sequencetyper()
        # Create reports
        printtime('Creating {} reports'.format(self.analysistype), self.start)
        self.reporter()

    def populate(self):
        import metagenomeFilter.createobject as createobject
        # Move the files to subfolders and create objects
        self.metadata = createobject.ObjectCreation(self)
        # Create and populate the .core attribute
        for sample in self.metadata.samples:
            setattr(sample, self.analysistype, GenObject())
            sample[self.analysistype].alleles = self.genes
            sample[self.analysistype].allelenames = [os.path.split(x)[1].split('.')[0] for x in self.genes]
            sample[self.analysistype].profile = self.profile
            sample[self.analysistype].alleledir = self.coregenelocation
            sample[self.analysistype].reportdir = os.path.join(sample.general.outputdirectory, self.analysistype)

    def profiler(self):
        """
        Creates a dictionary from the profile scheme(s)
        """
        # Initialise variables
        profiledata = defaultdict(make_dict)
        profileset = set()
        genedict = {}
        # Find all the unique profiles to use with a set
        for sample in self.metadata.samples:
            if sample[self.analysistype].profile != 'NA':
                profileset.add(sample[self.analysistype].profile[0])
        # Extract the profiles for each set
        for sequenceprofile in profileset:
            # Clear the list of genes
            genelist = []
            for sample in self.metadata.samples:
                if sequenceprofile == sample[self.analysistype].profile[0]:
                    genelist = [os.path.split(x)[1].split('.')[0] for x in sample[self.analysistype].alleles]
            try:
                # Open the sequence profile file as a dictionary
                profile = DictReader(open(sequenceprofile))
            # Revert to standard comma separated values
            except KeyError:
                # Open the sequence profile file as a dictionary
                profile = DictReader(open(sequenceprofile))
            # Iterate through the rows
            for row in profile:
                # Iterate through the genes
                for gene in genelist:
                    # Add the sequence profile, and type, the gene name and the allele number to the dictionary
                    try:
                        profiledata[sequenceprofile][row['ST']][gene] = row[gene]
                    except KeyError:
                        try:
                            profiledata[sequenceprofile][row['rST']][gene] = row[gene]
                        except KeyError:
                            raise

            # Add the gene list to a dictionary
            genedict[sequenceprofile] = sorted(genelist)
            # Add the profile data, and gene list to each sample
            for sample in self.metadata.samples:
                if sample.general.bestassembly != 'NA':
                    if sequenceprofile == sample[self.analysistype].profile[0]:
                        # Populate the metadata with the profile data
                        sample[self.analysistype].profiledata = profiledata[sample[self.analysistype].profile[0]]
                        # Add the allele directory to a list of directories used in this analysis
                        self.allelefolders.add(sample[self.analysistype].alleledir)
                        dotter()

    def annotatethreads(self):
        """
        Use prokka to annotate each strain
        """
        import metagenomeFilter.createobject as createobject
        # Move the files to subfolders and create objects
        self.runmetadata = createobject.ObjectCreation(self)
        # Fix headers
        self.headers()
        printtime('Performing prokka analyses', self.start)
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.annotate, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata.samples:
            # Create the prokka attribute in the metadata object
            setattr(sample, 'prokka', GenObject())
            # docker run -v /path/to/sequences:/path/to/sequences coregenome
            # prokka 2014-SEQ-0275.fasta --force --genus Escherichia --species coli --usegenus --addgenes
            # --prefix 2014-SEQ-0275 --locustag EC0275 --outputdir /path/to/sequences/2014-SEQ-0275/prokka
            sample.prokka.outputdir = os.path.join(sample.general.outputdirectory, 'prokka')
            # TODO Incorporate MASH/rMLST/user inputted genus, species results in the system call
            # Create the system call
            sample.prokka.command = 'docker run -v {}:{} {} ' \
                                    'prokka {} ' \
                                    '--force ' \
                                    '--genus {} ' \
                                    '--species {} ' \
                                    '--usegenus ' \
                                    '--addgenes ' \
                                    '--prefix {} ' \
                                    '--locustag {} ' \
                                    '--outdir {}' \
                .format(self.sequencepath, self.sequencepath, self.dockerimage, sample.general.fixedheaders,
                        self.genus, self.species, sample.name, sample.name, sample.prokka.outputdir)
            self.queue.put(sample)
        self.queue.join()

    def annotate(self):
        from subprocess import call
        while True:
            sample = self.queue.get()
            if not os.path.isfile('{}/{}.gff'.format(sample.prokka.outputdir, sample.name)):
                call(sample.prokka.command, shell=True, stdout=self.fnull, stderr=self.fnull)
            # List of the file extensions created with a prokka analysis
            files = ['err', 'faa', 'ffn', 'fna', 'fsa', 'gbk', 'gff', 'log', 'sqn', 'tbl', 'txt']
            # List of the files created for the sample by prokka
            prokkafiles = glob('{}/*'.format(sample.prokka.outputdir))
            # Find out which files have been created in the analysis
            for extension in files:
                # If the file was created, set the file path/name as the data for the attribute
                if extension in [prokka.split('.')[1] for prokka in prokkafiles]:
                    for output in prokkafiles:
                        setattr(sample.prokka, output.split('.')[1], output)
                # Otherwise, populate the attribute with 'NA'
                else:
                    setattr(sample.prokka, extension, 'NA')
            self.queue.task_done()

    def headers(self):
        """
        The contig ID must be twenty characters or fewer. The names of the headers created following SPAdes assembly
        are usually far too long. This renames them as the sample name
        """
        # from Bio import SeqIO
        for sample in self.metadata.samples:
            # Create an attribute to store the path/file name of the fasta file with fixed headers
            sample.general.fixedheaders = sample.general.bestassembly.replace('.fasta', '.ffn')
            # A list of contigs with modified record.id values
            fixedheaders = list()
            # Only do this if the file with fixed headers hasn't previously been created
            if not os.path.isfile(sample.general.fixedheaders):
                # Refseq genomes don't necessarily have underscores (or contig numbers) in the headers
                count = 0
                formatcount = '{:04d}'.format(count)
                for record in SeqIO.parse(open(sample.general.bestassembly, "rU"), "fasta"):
                    # Split off anything following the contig number
                    # >2013-SEQ-0129_1_length_303005_cov_13.1015_ID_17624 becomes
                    # >2013-SEQ-0129_1
                    record.id = record.id.split('_length')[0]
                    # print sample.name, record.id
                    # Prokka has a requirement that the header is unique and less than or equal to 20 characters
                    if len(record.id) > 20:
                        # Extract the contig number from the string - assumption is that this number is the final
                        # entry in the string, and that there are underscore separating the different components
                        contignumber = record.id.split('_')[-1] if '_' in record.id else formatcount
                        # Subtract the length of the contig number (and an additional one for the underscore) from
                        # 20 for the string slice, and add the contig number at the end
                        record.id = record.id[:(20 - len(contignumber) - 1)] + '_{}'.format(formatcount)
                    # Clear the name and description attributes of the record
                    record.name = ''
                    record.description = ''
                    # Add this record to our list
                    fixedheaders.append(record)
                # Open the filtered assembly file
                with open(sample.general.fixedheaders, 'wb') as formatted:
                    # Write the records in the list to the file
                    SeqIO.write(fixedheaders, formatted, 'fasta')

    def cdsthreads(self):
        """
        Determines which core genes from a pre-calculated database are present in each strain
        """
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.cds, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata.samples:
            #
            sample[self.analysistype].corepresence = dict()
            self.cdsqueue.put(sample)
        self.cdsqueue.join()

    def cds(self):
        while True:
            sample = self.cdsqueue.get()
            with open(sample.prokka.gff, 'rb') as gff:
                for feature in gff:
                    # Only interested in the sequence name if it is a CDS
                    if 'CDS' in feature:
                        # Extract the sequence name from the string. Example below
                        # 2013-SEQ-0123-2014_1	Prodigal:2.6	CDS	443	1741	.	+	0
                        # ID=0279_00002;Parent=0279_00002_gene;gene=kgtP_1;
                        # inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P0AEX3;
                        # locus_tag=0279_00002;product=Alpha-ketoglutarate permease
                        name = feature.split('ID=')[1].split(';')[0]
                        # Add number and names of genes to dictionaries
                        try:
                            gene = feature.split('gene=')[1].split(';')[0]
                            if gene in self.allelenames:
                                sample[self.analysistype].corepresence[name] = gene
                        except IndexError:
                                pass
            self.cdsqueue.task_done()

    def cdssequencethreads(self):
        """
        Extracts the sequence of each gene for each strain
        """
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.cdssequence, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata.samples:
            # Initialise a dictionary to store the sequence of each core gene
            sample[self.analysistype].coresequence = dict()
            self.sequencequeue.put(sample)
        self.sequencequeue.join()

    def cdssequence(self):
        while True:
            sample = self.sequencequeue.get()
            for record in SeqIO.parse(open(sample.prokka.ffn, 'rb'), 'fasta'):
                # If the gene name is present in the list of core genes, add the sequence to the dictionary
                if record.id in sample[self.analysistype].corepresence:
                    sample[self.analysistype].coresequence[record.id] = str(record.seq)
            self.sequencequeue.task_done()

    def allelematchthreads(self):
        """
        Determine allele of each gene
        """
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.allelematch, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata.samples:
            sample[self.analysistype].allelematches = dict()
            self.allelequeue.put(sample)
        self.allelequeue.join()

    def allelematch(self):
        while True:
            sample = self.allelequeue.get()
            # Iterate through all the core genes
            for name, gene in sample[self.analysistype].corepresence.items():
                # Iterate through all the alleles for the gene
                for record in SeqIO.parse(open(self.alleledict[gene], 'rb'), 'fasta'):
                    # If the current genes matches the database alleles
                    if str(record.seq) == sample[self.analysistype].coresequence[name]:
                        # Set the gene to the corresponding allele number
                        sample[self.analysistype].allelematches[gene] = record.id
            self.allelequeue.task_done()

    def sequencetyper(self):
        """
        Determines the sequence type of each strain based on comparisons to sequence type profiles
        """
        for sample in self.metadata.samples:
            if sample.general.bestassembly != 'NA':
                if type(sample[self.analysistype].allelenames) == list:
                    #
                    if sample[self.analysistype].profile != 'NA':
                        # Initialise dictionaries
                        sample[self.analysistype].profilematches = dict()
                        sample[self.analysistype].sequencetypematches = dict()
                        # Create the profiledata variable to avoid writing self.profiledata[self.analysistype]
                        profiledata = sample[self.analysistype].profiledata
                        # For each gene
                        for gene in sorted(sample[self.analysistype].allelenames):
                            try:
                                allelenumber = sample[self.analysistype].allelematches[gene].split('-')[1]
                                # Find the profile with the most alleles in common with the query genome
                                for sequencetype in profiledata:
                                    # refallele is the allele number of the sequence type
                                    refallele = profiledata[sequencetype][gene]
                                    if allelenumber == refallele:
                                        # Add matching alleles
                                        try:
                                            sample[self.analysistype].profilematches[sequencetype] += 1
                                            sample[self.analysistype].sequencetypematches[sequencetype].append(
                                                refallele)
                                        except KeyError:
                                            sample[self.analysistype].profilematches[sequencetype] = 1
                                            sample[self.analysistype].sequencetypematches[sequencetype] = list()
                                            sample[self.analysistype].sequencetypematches[sequencetype].append(
                                                refallele)
                            except KeyError:
                                pass

    def reporter(self):
        """
        Parse the results into a report
        """
        # Initialise variables
        header = ''
        row = ''
        for sample in self.metadata.samples:
            if sample[self.analysistype].reportdir != 'NA':
                if type(sample[self.analysistype].allelenames) == list:
                    # Populate the header with the appropriate data, including all the genes in the list of targets
                    header = 'Strain,SequenceType,Matches,Mismatches,NA,TotalGenes,{},\n' \
                        .format(','.join(sorted(sample[self.analysistype].allelenames)))
                sortedmatches = sorted(sample[self.analysistype].profilematches.items(),
                                       key=operator.itemgetter(1), reverse=True)[0]
                closestseqtype = sortedmatches[0]
                nummatches = int(sortedmatches[1])
                numna = 0
                queryallele = list()
                for gene, allele in sorted(sample[self.analysistype].profiledata[closestseqtype].items()):
                    #
                    try:
                        query = sample[self.analysistype].allelematches[gene].split('-')[1]
                        if allele == query:
                            queryallele.append(query)
                        else:
                            queryallele.append('{} ({})'.format(query, allele))
                    except KeyError:
                        queryallele.append('NA')
                        numna += 1
                mismatches = len(sample[self.analysistype].alleles) - nummatches - numna
                row += '{},{},{},{},{},{},{}'\
                    .format(sample.name, closestseqtype, nummatches, mismatches, numna,
                            len(sample[self.analysistype].alleles), ','.join(queryallele))
                row += '\n'
        # Create the report folder
        make_path(self.reportpath)
        # Create the report containing all the data from all samples
        with open('{}/{}.csv'.format(self.reportpath, self.analysistype), 'wb') \
                as combinedreport:
            # Write the results to this report
            combinedreport.write(header)
            combinedreport.write(row)

    def reporter1(self):
        """ Parse the results into a report"""
        # Initialise variables
        header = ''
        row = ''
        reportdirset = set()
        # Populate a set of all the report directories to use. A standard analysis will only have a single report
        # directory, while pipeline analyses will have as many report directories as there are assembled samples
        for sample in self.metadata.samples:
            # Ignore samples that lack a populated reportdir attribute
            if sample[self.analysistype].reportdir != 'NA':
                make_path(sample[self.analysistype].reportdir)
                # Add to the set - I probably could have used a counter here, but I decided against it
                reportdirset.add(sample[self.analysistype].reportdir)
        # Create a report for each sample from :self.resultprofile
        for sample in self.metadata.samples:
            if sample[self.analysistype].reportdir != 'NA':
                # row = ''
                if type(sample[self.analysistype].allelenames) == list:
                    # Populate the header with the appropriate data, including all the genes in the list of targets
                    header = 'Strain,SequenceType,Matches,TotalGenes,{},\n' \
                        .format(','.join(sorted(sample[self.analysistype].allelenames)))

                    # Set the sequence counter to 0. This will be used when a sample has multiple best sequence types.
                    # The name of the sample will not be written on subsequent rows in order to make the report clearer
                    seqcount = 0
                    # Iterate through the best sequence types for the sample (only occurs if update profile is disabled)
                    for seqtype in self.resultprofile[sample.name]:
                        sample[self.analysistype].sequencetype = seqtype
                        # The number of matches to the profile
                        matches = self.resultprofile[sample.name][seqtype].keys()[0]
                        # If this is the first of one or more sequence types, include the sample name
                        if seqcount == 0:
                            row += '{},{},{},{},'.format(sample.name, seqtype, matches,
                                                         len(sample[self.analysistype].alleles))
                        # Otherwise, skip the sample name
                        else:
                            row += ',{},{},{},'.format(seqtype, matches, len(sample[self.analysistype].alleles))
                        # Iterate through all the genes present in the analyses for the sample
                        for gene in sorted(sample[self.analysistype].allelenames):
                            # refallele = self.profiledata[self.analysistype][seqtype][gene]
                            refallele = sample[self.analysistype].profiledata[seqtype][gene]
                            # Set the allele and percent id from the dictionary's keys and values, respectively
                            allele = self.resultprofile[sample.name][seqtype][matches][gene].keys()[0]
                            percentid = self.resultprofile[sample.name][seqtype][matches][gene].values()[0]
                            if refallele and refallele != allele:
                                if 0 < float(percentid) < 100:
                                    row += '{} ({:.2f}%),'.format(allele, float(percentid))
                                else:
                                    row += '{} ({}),'.format(allele, refallele)
                            else:
                                # Add the allele and % id to the row (only add the percent identity if it is not 100%)
                                if 0 < float(percentid) < 100:
                                    row += '{} ({:.2f}%),'.format(allele, float(percentid))
                                else:
                                    row += '{},'.format(allele)
                        # Add a newline
                        row += '\n'
                        # Increment the number of sequence types observed for the sample
                        seqcount += 1
                    # If the length of the # of report directories is greater than 1 (script is being run as part of
                    # the assembly pipeline) make a report for each sample
                    if self.pipeline:
                        # Open the report
                        with open('{}{}_{}.csv'.format(sample[self.analysistype].reportdir, sample.name,
                                                       self.analysistype), 'wb') as report:
                            # Write the row to the report
                            report.write(row)
                dotter()
            # Create the report folder
            make_path(self.reportpath)
            # Create the report containing all the data from all samples
            with open('{}/{}.csv'.format(self.reportpath, self.analysistype), 'wb') \
                    as combinedreport:
                # Write the results to this report
                combinedreport.write(header)
                combinedreport.write(row)

    def __init__(self, inputobject):
        from Queue import Queue
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        self.start = inputobject.start
        self.cpus = inputobject.cpus
        self.genus = inputobject.genus
        self.species = inputobject.species
        self.dockerimage = inputobject.dockerimage
        self.pipeline = False
        self.metadata = MetadataObject()
        # Folders
        self.coregenelocation = os.path.join(self.path, 'coregenes', self.genus)
        self.profilelocation = os.path.join(self.path, 'profile', self.genus)
        self.reportpath = os.path.join(self.path, 'reports')
        # Class variables
        self.genes = sorted(glob('{}/*.fasta'.format(self.coregenelocation)))
        self.profile = glob('{}/*.txt'.format(self.profilelocation))
        self.analysistype = 'core'
        self.allelenames = sorted([os.path.basename(x).split('.')[0] for x in self.genes])
        self.alleledict = dict(zip(self.allelenames, self.genes))
        self.allelefolders = set()
        self.queue = Queue()
        self.dqueue = Queue()
        self.cdsqueue = Queue()
        self.sequencequeue = Queue()
        self.allelequeue = Queue()
        self.fnull = open(os.devnull, 'wb')
        self.resultprofile = defaultdict(make_dict)
        # Perform typing
        self.handler()
        # Remove the attributes from the object; they take up too much room on the .json report
        for sample in self.metadata.samples:
            try:
                delattr(sample[self.analysistype], "allelenames")
                delattr(sample[self.analysistype], "alleles")
                delattr(sample[self.analysistype], "profiledata")
            except KeyError:
                pass
        self.runmetadata = self.metadata
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)
