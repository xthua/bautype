#!/usr/bin/python
#-*- coding: UTF-8 -*-

import os
import sys
import time
import subprocess

def makeblastdb(infile):
    print ("start makeblastdb for input fasta at "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    cmd = "makeblastdb -dbtype nucl -in "+infile
    subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL)

def blastn(query,db,outdir):
    outfile=outdir+"/blast.reference.out"
    print("start blast to reference fasta at "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    os.system("blastn -task blastn -db " + db + ' -query ' + query + ' -out ' + outfile +'  -evalue 1e-5 -num_threads 2 -outfmt \"6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq sseq\"')
    return outfile

def tblastn(query,db,outdir):
    outfile=outdir+"/blast.gene.reference.out"
    print("start blast to gene reference fasta at "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    os.system("tblastn -db_gencode 11 -seg no -db " + db + ' -query ' + query + ' -out ' + outfile + ' -evalue 1e-5  -outfmt \"6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq sseq\"')
    return outfile
