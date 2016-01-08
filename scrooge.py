#!/usr/bin/env python
# -*- coding: utf-8 -*-


# ***************************************************************** #
# **         SCROOGE                                             ** #
# **         estimate genome size from single copy gene coverage ** #
# ***************************************************************** #



from Bio.Blast import NCBIXML
from Bio import SeqIO
import argparse
import sys
import os.path
import pysam

from scroogeclasses import *

def print_error(errormsg):
    print >> sys.stderr,"ERROR: %s"%errormsg
    parser.print_usage(sys.stderr)
    exit(1)

def tmpfile(filename):
    return os.path.join(args.tmpdir,filename)


def main():
    
    # ***************************************************************** #
    # **         command line parsing                                ** #
    # ***************************************************************** #
    global parser
    parser = argparse.ArgumentParser(description = "SCROOGE (Single Copy gene ROOted Genome size Estimation)\n  try to infer genome size from the coverage of single copy genes\n  by extrapolating the relation:\n\n  Genome size = (number reads) * (read length) / (coverage depth)",formatter_class = argparse.RawTextHelpFormatter)

    # global options: parameters for external files, temporary directory, etc ...
    parser.add_argument("-O","--optionfile",default="config.externalprograms.xml",
			help="XML file containing default options for external programs\n(default: 'config.externalprograms.xml')")
    parser.add_argument("-t","--tmpdir",default="tmp",
			help="Directory for temporary files\n(default: './tmp')")
    parser.add_argument("-T","--trashtmp",default=False,action="store_true",
			help="Trash temporary files\n(default: keep them)")
    parser.add_argument("-V","--nonverbose",action="store_true",default=False,
			help="Do not write information about current step to screen\n(default: write info)")

    # database for mining SCGs
    parser_db = parser.add_mutually_exclusive_group()
    parser_db.add_argument("-d","--dbfile",default=None,
			  help="Basename of database files for searching single copy genes\n(default: 'BUSCO')")
    parser_db.add_argument("-D","--dbseqfile",default=None,
			  help="Sequence file (FASTA format) to create database for searching single copy genes")

    # parameters of the algorithm itself
    parser.add_argument("-l","--cutofflength",type=int,default=100,
			help = "Cutoff length for SCGs\n (default:100)")
    parser.add_argument("-L","--readminlength",type=int,default=30,
			help = "Min length for reads\n(default: 30)")


    # datafiles for various stages of the algorithm
    parser.add_argument("-q","--query",default=None,
			help="FASTA file with sequences of single copy genes")
    parser.add_argument("-r","--reads",default=None,
			help="Reads from sequencing run")
    parser.add_argument("-c","--coveragefile",default="coverage.out",
			help="output file to write coverage depth")

    global args
    args = parser.parse_args()


    # ***************************************************************** #
    # **         create temporary data structure                     ** #
    # ***************************************************************** #
    if os.path.isdir(args.tmpdir):
	print >> sys.stderr,"WARNING: temporary directory '%s' exists. Files will be overwritten!"%args.tmpdir
    else:
	os.mkdir("tmp")



    # ***************************************************************** #
    # **         initialize blast search object                      ** #
    # ***************************************************************** #
    try:
	blastsearch = externalprogram(args.optionfile,"scgminingsearch",not args.nonverbose)
    except IOError:
	print_error("Could not find option file '%s'"%args.optionfile)
    except ValueError:
	print_error("Could not find options for step 'scgminingsearch' in file '%s'"%args.optionfile)
    except:
	print_error("weird error!")


    # ***************************************************************** #
    # **         create blast db                                     ** #
    # ***************************************************************** #
    if args.dbseqfile != None:
	try:
	    blastdb = externalprogram(args.optionfile,"scgminingcreatedb",not args.nonverbose)
	except IOError:
	    print_error("Could not find option file '%s'"%args.optionfile)
	except ValueError:
	    print_error("Could not find options for step 'scgminingcreatedb' in file '%s'"%args.optionfile)
	except:
	    print_error("weird error!")
	
	if os.path.isfile(args.dbseqfile):
	    blastdb.set_option("in",args.dbseqfile)
	    blastdb.set_option("out",tmpfile("SCGdb"),outfile=True)
	    blastsearch.set_option("db",tmpfile("SCGdb"))
	    blastdb.set_stderr(tmpfile("stderr.blastdb"))
	    blastdb.set_stdout(tmpfile("stdout.blastdb"))
	    
	    blastdb.execute()
	else:
	    print_error("could not find sequence file '%s' to create DB"%args.dbseqfile)
    # ***************************************************************** #
    # **         blast search                                        ** #
    # ***************************************************************** #
    elif args.dbfile != None:
	blastsearch.set_option("db",args.dbfile)

    # check for db files
    if not file_exists(blastsearch.get_files("db")):
	raise IOError
	    #print_error("could not find files")


    if args.query != None:
	if os.path.isfile(args.query):
	    blastsearch.set_option("query",args.query)
	else:
	    print_error("could not find query file")
    else:
	print_error("need query file")

    blastsearch.set_stderr(tmpfile("stderr.blastsearch"))
    blastsearch.set_stdout(tmpfile("stdout.blastsearch"))
    blastsearch.set_option("out",tmpfile(blastsearch.get_option("out")))
    blastsearch.execute()



    # ***************************************************************** #
    # **         initialize building hashfild for bowtie mapping     ** #
    # ***************************************************************** #
    try:
	bowtiebuild = externalprogram(args.optionfile,"generatehashfile",not args.nonverbose)
    except IOError:
	print_error("Could not find option file '%s'"%args.optionfile)
    except ValueError:
	print_error("Could not find options for step 'generatehashfile' in file '%s'"%args.optionfile)
    except:
	print_error("weird error!")

    bowtiebuild.add_parameter(tmpfile("GenomeSCG.fasta"))
    bowtiebuild.add_parameter(tmpfile("GenomeSCG.idx"))
    bowtiebuild.set_stderr(tmpfile("stderr.bowtiebuild"))
    bowtiebuild.set_stdout(tmpfile("stdout.bowtiebuild"))


    # ***************************************************************** #
    # **         building list of single copy genes                  ** #
    # **         by assigning correct identifiers for SCGs and       ** #
    # **         introducing cutoff in length                        ** #
    # ***************************************************************** #
    fa = open(args.query) # needed for correct names
    fblastxml = open(blastsearch.get_option("out")) # XML output from blast

    assembly_dict = SeqIO.to_dict(SeqIO.parse(fa,"fasta"))
    scg_blast_records = NCBIXML.parse(fblastxml)

    scg = SingleCopyGeneList(scglength = args.cutofflength, readlength = args.readminlength)

    for record in scg_blast_records:
	qid = record.query.split()[0]
	n=1
	for alignment in record.alignments:
	    for hsp in alignment.hsps:
		s = min(hsp.query_start,hsp.query_end)
		e = max(hsp.query_start,hsp.query_end)
		scg.add_sequence(qid+str(n),assembly_dict[qid].seq[s:e])
		n+=1

    fa.close()
    fblastxml.close()
    scg.write_sequence_file(bowtiebuild.get_parameters()[0])



    # ***************************************************************** #
    # **         generate hashfile for bowtie mapping                ** #
    # ***************************************************************** #

    bowtiebuild.execute()

    # ***************************************************************** #
    # **         bowtie mapping                                      ** #
    # ***************************************************************** #
    try:
	bowtie = externalprogram(args.optionfile,"mapping",not args.nonverbose)
    except IOError:
	print_error("Could not find option file '%s'"%args.optionfile)
    except ValueError:
	print_error("Could not find options for step 'generatehashfile' in file '%s'"%args.optionfile)
    except:
	print_error("weird error!")

    bowtie.set_option("x",bowtiebuild.get_parameters()[1])
    bowtie.set_option("U",args.reads)
    bowtie.set_option("S",tmpfile("mapping.sam"))

    bowtie.set_stderr(tmpfile("stderr.bowtie"))
    bowtie.set_stdout(tmpfile("stdout.bowtie"))

    bowtie.execute()


    # ***************************************************************** #
    # **         get coverage                                        ** #
    # ***************************************************************** #
    samfile = pysam.Samfile(bowtie.get_option("S"),"r")

    for alignment in samfile.fetch():
	if alignment.tid >= 0: # (id == -1) -> read not mapped
	    scg[alignment.reference_name].add_coverage(alignment.reference_start,alignment.reference_end)

    scg.write_coverage_file(args.coveragefile)



    # ***************************************************************** #
    # **         delete temporary folder with all files in it        ** #
    # ***************************************************************** #
    if args.trashtmp:
	from shutil import rmtree
	rmtree(args.tmpdir,ignore_errors=True)




if __name__ == "__main__":
  main()


