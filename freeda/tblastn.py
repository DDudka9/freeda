#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:21:53 2021

@author: damian

Runs blast using NCBI tblastn.

Fasta genome names need to be stored in "genomes.txt" file one per line using one underscore:
    xxxxx_xxxxx.fasta
    yyyyy_yyyyy.fasta
    
Gene names need to be stored in "genes.txt" file one per line:
    aaaaa
    bbbbb

Outputs tabulated text files per gene per genomes.

"""

from freeda import genomes_preprocessing, pyinstaller_compatibility
import subprocess
import os
import glob
import time
import logging

#logging.basicConfig(filename="tblast_temp.log", level=logging.INFO, format="%(message)s")


def run_blast(wdir, ref_species, all_genes):
    """Runs tblastn based on NCBI makedatabase routine."""

    database_path = wdir + "Genomes/"
    query_path = wdir + "Blast_input/"
    output_path = wdir + "Blast_output/"
    form = "6 qseqid means sseqid means qstart means qend means sstart means send means evalue means " \
           "bitscore length means pident means mismatch means gapopen means qlen means slen means"

    all_genomes = genomes_preprocessing.get_names(wdir, ref_species)
    genomes = [names[1] for names in all_genomes]

    # clear Blast_output folder
    all_old_blast_output_files = glob.glob(os.path.join(output_path, "*.txt"))
    for file in all_old_blast_output_files:
        os.remove(file)

    # make sure database for each genome already exists or is successfully built
    for genome in genomes:
        genome_file_database = check_genome_present(wdir, ref_species, database_path, genome, ref_genome=False)
        # failed to build database
        if genome_file_database is False:
            return None

    # PYINSTALLER: Add path to os path variable.
    tblastn_path = pyinstaller_compatibility.resource_path("tblastn")

    # perform blast
    for genome in genomes:
        for gene in all_genes:
            database = database_path + genome + ".fasta"
            query = query_path + gene + "_" + ref_species + "_protein.fasta"
            output = output_path + gene + "_" + genome + ".txt"
            to_blast = [tblastn_path, "-db", database, "-query", query, "-out", output, "-outfmt",
                        form, "-num_threads", "8"]
            message = "\nPerforming tblastn for gene: %s from genome: %s\n" % (gene, genome)
            #print("\nPerforming tblastn for gene: %s from genome: %s\n" % (gene, genome))
            logging.info(message)
            subprocess.call(to_blast)

    print("\ntblastn txt files have been generated.")
    
    return


def check_genome_present(wdir, ref_species, database_path, genome, ref_genome=False):
    """Checks if a given genome is present. Unpacks and unzips genomes downloaded from NCBI Assembly.
    Non-ncbi assemblies must be prepared as ".fasta" files conform with "genomes.txt" names.
    It also looks for reference genome if key-only argument reference_genome is invoked."""

    genome_file_database = True

    # if function called by input extractor module not tblastn module
    if ref_genome is True:
        ref_genome_path = database_path
    
    # define expected files
    zip_file = genome + ".zip"
    expected_genome_file = genome + ".fasta"
    # database extensions for tblastn
    expected_database_extensions = {#".nsq",    # cheetah genome has a bi different extensions
                                    #".nhr",
                                    #".nin",
                                    #".00.phr",
                                    #".00.pin",
                                    #".00.psq",
                                    #".01.phr",
                                    #".01.pin",
                                    #".01.psq",
                                    #".02.phr",
                                    #".02.pin",
                                    #".02.psq",
                                    ".pal"}
    
    # get info on all files available
    all_files = []
    for root, dirs, files in os.walk(database_path, topdown=False):
        for f in files:
            all_files.append(f)

    # there are no genomes in the folder
    if not all_files:
        all_files.append("empty")
    
    # check if genome database is present
    for file in all_files:
        # look for indices
        if expected_genome_file + ".pal" not in all_files and expected_genome_file + ".nin" not in all_files:
            genome_file_database = False

    # move on to next genome if database present
    if ref_genome is False and genome_file_database is True:
        
        message = "\nGenome : %s blast database already exists" % genome
        logging.info(message)
        #print("\nGenome : %s blast database already exists" % genome)
        
        return genome_file_database
    
    # build database on existing genome fasta file if present
    if ref_genome is False and expected_genome_file in all_files:

        message = "\nNOT detected database for : %s but %s file is present" \
                             " -> building database...\n" % (genome, expected_genome_file)
        logging.info(message)
        #print("\nNOT detected database for : %s but %s file is present" \
        #                     " -> building database...\n" % (genome, expected_genome_file))
        
        # make database
        make_blast_database(database_path, genome)
        # make sure to assign back this variable to TRUE
        genome_file_database = True
        
        return genome_file_database
    
    # download the genome as a tar file if both database and tar are missing
    if ref_genome is False and genome_file_database is False and zip_file not in all_files:

        message = "\nGenome : %s blast database does not exists" \
              " -> downloading and decompressing the genome (it might take couple of minutes)...\n" % genome
        logging.info(message)
        #print("\nGenome : %s blast database does not exists"
        #      " -> downloading and decompressing the genome (it might take couple of minutes)...\n" % genome)

        all_genomes = genomes_preprocessing.get_names(wdir, ref_species)
        accession_nr = [names[2] for names in all_genomes if genome in names][0]
        # download genome
        download_genome(genome, accession_nr, database_path)
        # decompress genome
        #unzip_genome(wdir, genome, accession_nr, database_path)
        # make database
        make_blast_database(database_path, genome)

        # check if genome was downloaded and unpacked successfully
        genome_found = False
        for root, dirs, files in os.walk(database_path, topdown=False):
            for file in files:
                if file == genome + ".fasta":
                    genome_found = True

        # exit pipeline if fasta genome absent
        if genome_found is False:
            genome_file_database = False
            message = "...FATAL_ERROR... : Genome : %s failed to download or decompress" \
                  " -> exciting the pipeline now...\n" % genome
            logging.info(message)
            #print("...FATAL_ERROR... : Genome : %s failed to download or decompress"
            #      " -> exciting the pipeline now...\n" % genome)
            return genome_file_database

        # validate that the database was generated
        all_files = []
        for root, dirs, files in os.walk(database_path, topdown=False):
            for f in files:
                all_files.append(f)

        for file in all_files:
            # look for indices
            if expected_genome_file + ".pal" not in all_files and expected_genome_file + ".nin" not in all_files:
                message = "\n...FATAL_ERROR... : Genome : %s database failed to build" \
                      " -> exciting the pipeline now...\n" % genome
                logging.info(message)
                #print("\n...FATAL_ERROR... : Genome : %s database failed to build"
                #      " -> exciting the pipeline now...\n" % genome)
                genome_file_database = False
                return genome_file_database

        genome_file_database = True

        # fasta file and databases are present for the genome
        if genome_found is True and genome_file_database is True:
            genome_file_database = True
            message = "\nGenome : %s was downloaded and decompressed successfully" % genome
            logging.info(message)
            #print("\nGenome : %s was downloaded and decompressed successfully" % genome)

            return genome_file_database


    # unpack and unzip if tar file present
    #if reference_genome is False and genome_file_database is False and tar_file in all_files:
        
    #    print("\nNOT detected database for : %s -> building database...\n " % genome)
        
    #    # unpack and decompress the genome into a fasta file
    #    genome_to_unzip = database_path + "/" + genome + ".tar"
    #    unzip_genome(database_path, genome, genome_to_unzip)
    #    # remove the tar file
    #    os.remove(database_path + tar_file)
    #    # make database
    #    make_blast_database(database_path, genome)
    #    # make sure to tick back this parameter to TRUE
    #    genome_file_database = True
    
    # unpack and decompress reference genome
    if ref_genome:
        
        if expected_genome_file in all_files:

            message = "\nReference genome : %s is present\n" % genome
            logging.info(message)
            #print("\nReference genome : %s is present\n" % genome)

            return True
            
        elif expected_genome_file not in all_files:

            message = "\nReference genome : %s does not exists" \
                              " -> downloading and decompressing it now...\n" % genome
            logging.info(message)
            #print("\nReference genome : %s does not exists" \
            #                  " -> downloading and decompressing it now...\n" % genome)

            # unpack and decompress the genome into a fasta file
            #genome_to_unzip = ref_genome_path + "/" + genome + "zip"
            all_genomes = genomes_preprocessing.get_names(wdir, ref_species, ref_genome=True)
            accession_nr = all_genomes[2]
            # download genome
            download_genome(genome, accession_nr, ref_genome_path)
            # decompress genome
            #unzip_genome(wdir, genome, accession_nr, ref_genome_path)
            # remove the tar file
            #os.remove(ref_genome_path + tar_file)
            return True
        
        #elif expected_genome_file not in all_files and tar_file not in all_files:y
        #    print("\n...FATAL_ERROR... : Reference genome : %s does not exists and archive file is missing" \
        #                      " -> exiting the pipeline now...\n" % genome)
        #    return False

    else:
        message = "Else statement when making blast database"
        logging.info(message)
        print("Else statement when making blast database")
        
    return genome_file_database


def download_genome(genome, accession_nr, database_path):
    """Downloads and unzips genomes using NCBI Datasets CLI, makes one fasta file"""

    # download genome
    #print("         Downloading and unzipping genome : %s - it might take a couple of min ..." % genome)
    start_time = time.time()
    filepath_1 = database_path + genome + ".zip"
    # need to exclude genomic cds because its also ".fna" which confuses the concatenation
    cmd1 = ["datasets", "download", "genome", "accession", accession_nr, "--exclude-genomic-cds", "--filename",
            filepath_1]
    subprocess.call(cmd1, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb')) # mute output and cautions

    # unzip all chromosomes into single file
    filepath_2 = database_path + genome + ".fasta"
    cmd2 = ["unzip", "-pq", filepath_1, "cat", "*/*.fna", filepath_2]
    with open(filepath_2, "w") as outfile:
        subprocess.call(cmd2, stdout=outfile, stderr=open(os.devnull, 'wb')) # redirect output to file and mute cautions

    stop_time = time.time()
    message = "         -> Done : in %s min" % ((stop_time - start_time) / 60)
    logging.info(message)
    #print("         -> Done : in %s min" % ((stop_time - start_time) / 60))

    os.remove(filepath_1)


def make_blast_database(database_path, genome):
    """ Makes a blast database based on genome fasta file."""

    genome_file = genome + ".fasta"
    make_database_nucl = ["makeblastdb", "-in", database_path + genome_file, "-dbtype", "nucl"]
    make_database_prot = ["makeblastdb", "-in", database_path + genome_file, "-dbtype", "prot"]

    message = "\n                  Building blast database for genome : %s ..." % genome
    logging.info(message)
    #print("\n                  Building blast database for genome : %s ..." % genome)
    subprocess.call(make_database_nucl, stdout=open(os.devnull, 'wb')) # mute log info
    subprocess.call(make_database_prot, stdout=open(os.devnull, 'wb')) # mute log info


"""

def download_genome(wdir, genome, accession_nr, database_path, ref_genome=False):
    #Downloads preselected genomes from NCBO Assemblies using NCBI API : ncbi.datasets.

    # Almost empty zip downloaded for : MUSCULUS_genome, SPRETUS_genome,

    # INSTALLING pip install ncbi-datasets-pylib==12.13.1 requires -> pip install wheel
    # and probably also pip install google-cloud-bigquery and pip install pandas-gbq

    import ncbi.datasets.openapi
    from pprint import pprint
    from ncbi.datasets.openapi.api import gene_api
    from ncbi.datasets.openapi.model.v1_fasta import V1Fasta

    # Defining the host is optional and defaults to https://api.ncbi.nlm.nih.gov/datasets/v1
    # See configuration.py for a list of all supported configuration parameters.
    configuration = ncbi.datasets.openapi.Configuration(host="https://api.ncbi.nlm.nih.gov/datasets/v1")

    # The client must configure the authentication and authorization parameters
    # in accordance with the API server security policy.
    # Examples for each auth method are provided below, use the example that
    # satisfies your auth use case.

    # Configure API key authorization: ApiKeyAuthHeader
    configuration.api_key['ApiKeyAuthHeader'] = 'YOUR_API_KEY'

    # Uncomment below to setup prefix (e.g. Bearer) for API key, if needed
    # configuration.api_key_prefix['ApiKeyAuthHeader'] = 'Bearer'

    # Enter a context with an instance of the API client
    with ncbi.datasets.openapi.ApiClient(configuration) as api_client:
        # Create an instance of the API class
        api_instance = gene_api.GeneApi(api_client)
        gene_ids = [59067,]
        # [int] | NCBI gene ids

    include_annotation_type = [V1Fasta("FASTA_UNSPECIFIED"),]
    # [V1Fasta] | Select additional types of annotation to include in the data package.  If unset, no annotation is provided. (optional)

    fasta_filter = ["fasta_filter_example",]
    # [str] | Limit the FASTA sequences in the datasets package to these transcript and protein accessions (optional)

    filename = "ncbi_dataset.zip"
    # str | Output file name. (optional) (default to "ncbi_dataset.zip")

    try:
        # Get a gene dataset by gene ID
        api_response = api_instance.download_gene_package(gene_ids, include_annotation_type=include_annotation_type,
                                                          fasta_filter=fasta_filter, filename=filename)
        pprint(api_response)
    except ncbi.datasets.openapi.ApiException as e:
        print("Exception when calling GeneApi->download_gene_package: %s\n" % e)

    api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())

    # get genomes
    assembly_accessions = [accession_nr]
    print(accession_nr)
    exclude_sequence = False
    #include_annotation_type = ['PROT_FASTA']
    api_response = api_instance.download_assembly_package(assembly_accessions,
                                                        exclude_sequence=exclude_sequence,
                                                        #include_annotation_type=include_annotation_type,
                                                        # Because we are streaming back the results to disk,
                                                        # we should defer reading/decoding the response
                                                        _preload_content=False)

    # write genome to a zip file
    with open(wdir + genome + ".zip", 'wb') as f:
        f.write(api_response.data)

    # move the zip file into either Genomes folder or Reference_genomes folder
    #if ref_genome is False:
    #    shutil.move(wdir + genome + ".zip", database_path + genome + ".zip")
    #else:
    #    shutil.move(wdir + genome + ".zip", wdir + "/Reference_genomes/" + genome + ".zip")

def unzip_genome(wdir, genome, accession_nr, database_path):
    #Unzips NCBI API downloaded genomes in zip format using ... gzip and shutil modules

    import os
    import shutil
    from zipfile import ZipFile

    # unzip the genome and move into Genomes
    to_unzip = wdir + genome + ".zip"
    with ZipFile(to_unzip, 'r') as z:
        # makes a list of files inside a zipped file
        all_names = z.namelist()
        for file_name in all_names:
            # find the genome fasta (.fna) file and extract it into temporary folder
            if file_name.endswith('.fna'):
                temp = database_path + genome
                z.extract(file_name, temp)

    # get the genome file, rename and move directly into Genomes folder
    for root, dirs, files in os.walk(temp, topdown=False):
        for file in files:
            if accession_nr in file:
                #print(file)
                #print(accession_nr)
                #print(root)
                shutil.move(root + "/" + file, database_path + genome + ".fasta")

    # removes directory and its content
    #shutil.rmtree(temp)
    #os.remove(to_unzip)


("\n"
 "            # to avoid duplications\n"
 "            if [string for string in all_files if pattern in string] == []:\n"
 "                \n"
 "                gene = gene.lstrip(\"\n\").rstrip(\"\n\")\n"
 "                database = database_path + genome\n"
 "                query = query_path + gene + \"_\" + ref_species + \"_protein.fasta\"\n"
 "                output = output_path + gene + \"_\" + genome + \".txt\"\n"
 "                to_blast = [\"tblastn\", \"-db\", database, \"-query\", query, \"-out\", output, \"-outfmt\", form, \"-num_threads\", \"5\"]\n"
 "                subprocess.call(to_blast)\n"
 "            \n"
 "            else:\n"
 "                print(\"\nGene: %s from genome: %s -> was already blasted.\" % (gene, genome))    \n"
 "\n"
 "\n"
 "\n"
 " import tarfile
 " import shutil
 " import gzip
 "def unzip_genome(database_path, genome, genome_to_unzip):\n"
 "    #Unzips tar archive files using gzip and shutil modules\n"
 "    \n"
 "    tar = tarfile.open(genome_to_unzip)\n"
 "    for tarinfo in tar:\n"
 "        if os.path.splitext(tarinfo.name)[1] == \".gz\":\n"
 "            source_directory = database_path + tarinfo.name\n"
 "            filename = tarinfo.name.split(\"/\")[1]\n"
 "    \n"
 "    tar.extractall(path=database_path, members=find_packed_genome(tar, genome))\n"
 "    tar.close()\n"
 "    \n"
 "    shutil.move(source_directory, database_path)\n"
 "    # remove the now empty directory\n"
 "    moved_gz_filename = source_directory.split(\"/\")[-1]\n"
 "    empty_directory = source_directory.rstrip(moved_gz_filename)\n"
 "    os.rmdir(empty_directory)\n"
 "    # open the genome content in text mode (\"rt\")\n"
 "    with gzip.open(database_path + filename, 'rt') as f_in:\n"
 "        print(\"\n        Detected compressed genome file for : %s -> decompressing...\n \" % genome)\n"
 "\n"
 "        # write it into a fasta file\n"
 "        with open(database_path + genome + '.fasta', 'wt') as f_out:\n"
 "            shutil.copyfileobj(f_in, f_out)\n"
 "    \n"
 "    # remove compressed file\n"
 "    os.remove(database_path + filename)\n"
 "            \n"
 "\n"
 "\n"
 "def find_packed_genome(members, genome):\n"
 "    #Finds a packed genome in the archive file using tarfile module.\n"
 "    \n"
 "    for tarinfo in members:\n"
 "        if os.path.splitext(tarinfo.name)[1] == \".gz\":\n"
 "            print(\"\n    Detected archive genome file for : %s -> unpacking...\n \" % genome)\n"
 "            \n"
 "            yield tarinfo\n"
 "\n"
 "DEPRECATED WAY TO OPEN MANUALLY DOWNLOADED ASSEMBLIES\n"
 "\n"
 "import subprocess\n"
 "import os\n"
 "import glob\n"
 "import tarfile\n"
 "import shutil\n"
 "import gzip\n"
 "\n"
 "def check_genome_present(database_path, genome, reference_genome=False):\n"
 "    \n"
 "    genome_file_database = True\n"
 "\n"
 "    if reference_genome == True:\n"
 "        reference_genome_path = database_path\n"
 "\n"
 "    # define expected files\n"
 "    tar_file = genome + \".tar\"\n"
 "    expected_genome_file = genome + \".fasta\"\n"
 "    expected_database_extensions = {\".00.phr\",\n"
 "                           \".00.pin\",\n"
 "                           \".00.psq\",\n"
 "                           \".01.phr\",\n"
 "                           \".01.pin\",\n"
 "                           \".01.psq\",\n"
 "                           \".02.phr\",\n"
 "                           \".02.pin\",\n"
 "                           \".02.psq\",\n"
 "                           \".nhr\",\n"
 "                           \".nin\",\n"
 "                           \".nsq\",\n"
 "                           \".pal\"}\n"
 "\n"
 "    # get info on all tar files available\n"
 "    all_files = []\n"
 "    for root, dirs, files in os.walk(database_path, topdown=False):\n"
 "        for f in files:\n"
 "            all_files.append(f)\n"
 "\n"
 "    # check if genome database is present\n"
 "    for file in all_files:\n"
 "        for extension in expected_database_extensions:\n"
 "            if expected_genome_file + extension not in all_files:\n"
 "                genome_file_database = False\n"
 "\n"
 "\n"
 "    # move on to next genome if database present\n"
 "    if reference_genome == False and genome_file_database == True:\n"
 "\n"
 "        print(\"\nGenome : %s blast database already exists\" % genome)\n"
 "\n"
 "        return genome_file_database\n"
 "\n"
 "    # build datbase on existing genome fasta file if present\n"
 "    elif reference_genome == False and expected_genome_file in all_files:\n"
 "\n"
 "        print(\"\nNOT detected database for : %s but %s file is present\" \\n"
 "                             \" -> building database...\n\" % (genome, expected_genome_file))\n"
 "\n"
 "        # make database\n"
 "        make_blast_database(database_path, genome)\n"
 "        # make sure to tick back this parameter to TRUE\n"
 "        genome_file_database = True\n"
 "\n"
 "        return genome_file_database\n"
 "\n"
 "    # exit pipeline if both database and tar are missing\n"
 "    elif reference_genome == False and genome_file_database == False and tar_file not in all_files:\n"
 "\n"
 "        print(\"\nGenome : %s blast database does not exists and archive file is missing\" \\n"
 "                            \" -> exiting the pipeline now...\n\" % genome)\n"
 "\n"
 "        return genome_file_database\n"
 "\n"
 "    # unpack and unzip if tar file present\n"
 "    elif reference_genome == False and genome_file_database == False and tar_file in all_files:\n"
 "\n"
 "        print(\"\nNOT detected database for : %s -> building database...\n \" % genome)\n"
 "\n"
 "        # unpack and decompress the genome into a fasta file\n"
 "        genome_to_unzip = database_path + \"/\" + genome + \".tar\"\n"
 "        unzip_genome(database_path, genome, genome_to_unzip)\n"
 "        # remove the tar file\n"
 "        os.remove(database_path + tar_file)\n"
 "        # make database\n"
 "        make_blast_database(database_path, genome)\n"
 "        # make sure to tick back this parameter to TRUE\n"
 "        genome_file_database = True\n"
 "\n"
 "    # unpack and decompress reference genome\n"
 "    elif reference_genome == True:\n"
 "\n"
 "        if expected_genome_file in all_files:\n"
 "            print(\"\nReference genome : %s is present\n\" % genome)\n"
 "            return True\n"
 "\n"
 "        elif expected_genome_file not in all_files and tar_file in all_files:\n"
 "\n"
 "            print(\"\nReference genome : %s does not exists but archive file is present\" \\n"
 "                              \" -> unpacking and decompressing...\n\" % genome)          \n"
 "\n"
 "            # unpack and decompress the genome into a fasta file\n"
 "            genome_to_unzip = reference_genome_path + \"/\" + genome + \".tar\"\n"
 "            unzip_genome(reference_genome_path, genome, genome_to_unzip)\n"
 "            # remove the tar file\n"
 "            os.remove(reference_genome_path + tar_file)\n"
 "            return True\n"
 "\n"
 "        elif expected_genome_file not in all_files and tar_file not in all_files:\n"
 "            print(\"\nReference genome : %s does not exists and archive file is missing\" \\n"
 "                              \" -> exiting the pipeline now...\n\" % genome)\n"
 "            return False\n"
 "\n"
 "    else:\n"
 "        print(\"Else statement when making blast database\")\n"
 "\n"
 "    return genome_file_database\n"
 "\n"
 "\n"
 "def unzip_genome(database_path, genome, genome_to_unzip):\n"
 "    #Unzips tar archive files using gzip and shutil modules\n"
 "\n"
 "    tar = tarfile.open(genome_to_unzip)\n"
 "    for tarinfo in tar:\n"
 "        if os.path.splitext(tarinfo.name)[1] == \".gz\":\n"
 "            source_directory = database_path + tarinfo.name\n"
 "            filename = tarinfo.name.split(\"/\")[1]\n"
 "\n"
 "    tar.extractall(path=database_path, members=find_packed_genome(tar, genome))\n"
 "    tar.close()\n"
 "\n"
 "    shutil.move(source_directory, database_path)\n"
 "    # remove the now empty directory\n"
 "    moved_gz_filename = source_directory.split(\"/\")[-1]\n"
 "    empty_directory = source_directory.rstrip(moved_gz_filename)\n"
 "    os.rmdir(empty_directory)\n"
 "    # open the genome content in text mode (\"rt\")\n"
 "    with gzip.open(database_path + filename, 'rt') as f_in:\n"
 "        print(\"\n        Detected compressed genome file for : %s -> decompressing...\n \" % genome)\n"
 "\n"
 "        # write it into a fasta file\n"
 "        with open(database_path + genome + '.fasta', 'wt') as f_out:\n"
 "            shutil.copyfileobj(f_in, f_out)\n"
 "\n"
 "    # remove compressed file\n"
 "    os.remove(database_path + filename)\n"
 "\n"
 "\n"
 "def find_packed_genome(members, genome):\n"
 "    #Finds a packed genome in the archive file using tarfile module\n"
 "\n"
 "    for tarinfo in members:\n"
 "        if os.path.splitext(tarinfo.name)[1] == \".gz\":\n"
 "            print(\"\n    Detected archive genome file for : %s -> unpacking...\n \" % genome)\n"
 "\n"
 "            yield tarinfo\n"
 "\n"
 "\n"
 "def make_blast_database(database_path, genome):\n"
 "    #Makes a blast database based on genome fasta file\n"
 "\n"
 "    genome_file = genome + \".fasta\"\n"
 "    make_database_nucl = [\"makeblastdb\", \"-in\", database_path + genome_file, \"-dbtype\", \"nucl\"]\n"
 "    subprocess.call(make_database_nucl)\n"
 "    make_database_prot = [\"makeblastdb\", \"-in\", database_path + genome_file, \"-dbtype\", \"prot\"]\n"
 "    subprocess.call(make_database_prot)\n"
 "    \n")"""