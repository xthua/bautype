#!/usr/bin/python
#-*- coding: UTF-8 -*-
"""
Author:liangqian at 20190806 ;lq4977@tkgeneclub.com
Name:ABtypingTools
Usage: This program is the pipeline of LPS or CPS typing for Acinetobacter baumannii

-i,--input  the input file of a acinetobacter baumannii genome fasta
-r,--reference  Prefix reference loci file in reference_data directory; 
              Abaumannii_KL_reference(KPS) 
              Abaumannii_OCL_reference(CPS)
-o,--outdir	output directory
-m,----tempdir	temp directory
"""

import os
import sys
import argparse
from argparse import ArgumentParser
import datetime
import time
from Bio import SeqIO
from plot_gene import plot_gene_coverage
from Blast_format import blastn,tblastn,makeblastdb
import os.path
import shutil


def main(args):
    bindir = sys.path[0]
    outdir = args.outdir;make_dir(outdir);outdir = os.path.abspath(outdir)
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    temp_dir = args.tempdir;make_dir(temp_dir);temp_dir = os.path.abspath(temp_dir)
    infasta = args.input
    indir = args.indir
    ref_fasta = args.reference+".fasta"
    ref_fasta_len = args.reference+".fasta.len"
    ref_gene = args.reference+".pep.fasta"
    ref_gene_pos = args.reference+".gene.pos"
    if infasta and indir:
        print("Please choose one parameter\(-d or -i\) to run")
        sys.exit(0)
    if infasta:
        in_genedir = outdir+"/ingene";make_dir(in_genedir);
        out_genedir = outdir+"/outgene";make_dir(out_genedir);
        runFlag = checkFileFormat(infasta)
        if runFlag:
            try:
                print ("Analysing file",infasta,"at "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
                run(bindir,outdir,infasta,in_genedir,out_genedir,temp_dir,indir,ref_fasta,ref_fasta_len,ref_gene,ref_gene_pos)
            except:
                print("Please check file: ",infasta)
        else:
            print("Please check file: ",infasta)
    if indir:
        for file in os.listdir(indir):
            infile1 = os.path.join(indir,file)
            name = file.split('.')[0]
            outdir1 = outdir + '/' + time.strftime('%Y%m%d%H%M%S',time.localtime())+'_'+name
            make_dir(outdir1)
            infasta = outdir+'/'+file
            shutil.copyfile(infile1,infasta)
            in_genedir = outdir1+"/ingene";make_dir(in_genedir);
            out_genedir = outdir1+"/outgene";make_dir(out_genedir);
            runFlag = checkFileFormat(infasta)
            if runFlag:
                try:
                    print ("Analysing file",infasta," at "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
                    run(bindir,outdir1,infasta,in_genedir,out_genedir,temp_dir,indir,ref_fasta,ref_fasta_len,ref_gene,ref_gene_pos)
                except:
                    print("Please check file: ",infasta)
            else:
                print("Please check file: ",infasta)

def is_chinese(string):
    for ch in string:
        if u'\u4e00' <= ch <= u'\u9fff':
            return True
    return False

def checkFileFormat(infile):
    try:
        for rec in SeqIO.parse(infile,'fasta'):
            id = rec.description
            isChinese = is_chinese(id)
            if isChinese:
                print("invalid character found")
                return False
            seq = rec.seq
            if seq:
                return True
    except:
        return False
    return False

def run(bindir,outdir,infasta,in_genedir,out_genedir,temp_dir,indir,ref_fasta,ref_fasta_len,ref_gene,ref_gene_pos):
    ## makeblastdb input fasta
    makeblastdb(infasta)
    keyname=os.path.basename(infasta)
    ## blastn to locus reference sequence
    refblast=blastn(ref_fasta,infasta,temp_dir)
    geneblast=tblastn(ref_gene,infasta,temp_dir)
    hitout,allhits,q_len=format_blast(refblast)
    hitout_gene,allhits_gene,q_gene_len=format_blast(geneblast)
    best_k_ref,sign1,sign2,output_match,best_blast=get_best_k_match(allhits,q_len)
    (match_seq,match_blastout)=output_locus_match_input(best_k_ref,best_blast,infasta)
    #print(best_k_ref)
    print(output_match)
    savefile(outdir+"/output_match_blast.txt",match_blastout)
    savefile(outdir+"/output_match_seq.fa",match_seq)
    ingene,ingene_cov,other_gene,all_ingene=get_gene_match(allhits_gene,q_gene_len,best_k_ref)
    out_match_format=out_format(keyname,best_k_ref,ref_gene_pos,output_match,ingene_cov)
    savefile(outdir+"/output_match_result.txt",out_match_format)
    if all_ingene !="None":
        gene_content_all="Gene Name\tCoverage(%)\tIdentity(%)\n"+"\n".join("\t".join(str(s) for s in ss) for ss in all_ingene)
        savefile(outdir+"/output_allInLocus_genematch.txt",gene_content_all)
    if ingene !="None":
        gene_content="Gene Name\tLocus\tIdentity(%)\tMatch seq in Gene\tGene Length(bp)\tGenome ID\tPosition in Genome\n"+"\n".join("\t".join(str(s) for s in ss) for ss in ingene)
        savefile(outdir+"/output_InLocus_genematch.txt",gene_content)
        plot_gene_coverage(best_k_ref,ref_gene_pos,all_ingene,ingene_cov,outdir)
        get_geneSeq(ref_gene,ingene,in_genedir)
    if other_gene !="None":
        gene_content="Gene Name\tLocus\tIdentity(%)\tMatch seq in Gene\tGene Length(bp)\tGenome ID \tPosition in Genome\n"+"\n".join("\t".join(str(s) for s in ss) for ss in other_gene)
        savefile(outdir+"/output_otherLocus_genematch.txt",gene_content)
        get_geneSeq(ref_gene,other_gene,out_genedir)
    ###plot gene Collinea when the Sign is perfect
    tag=out_match_format.split('\n')[1].split('\t')[2]
    if tag=="perfect":
        Collineardir=outdir+"/Collinear"
        make_dir(Collineardir)
        os.system("python "+bindir+"/plot_gene_Collinear.py "+outdir+"/output_match_blast.txt "+temp_dir+"/blast.gene.reference.out "+ref_gene_pos+" "+ref_fasta_len+" "+Collineardir+" Match Sequence")

def get_geneSeq(ref_gene,genearray,outdir):
    gene_key={}
    for ss in genearray:
        key,cov,id,pos,lenght,subjectid,subject_ranges=ss
        gene_key[key]=ss 
    with open(ref_gene,'r') as f1:
        for record in SeqIO.parse(f1,"fasta"):
             if str(record.id) in gene_key.keys():
                 filename=outdir+"/"+str(record.id)+".fasta"
                 SeqIO.write(record, filename, "fasta")

def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
def tid_maker():
    return '{0:%H%M%S}'.format(datetime.datetime.now())

def savefile(filename,contents):
    fh = open(filename, 'w')
    fh.write(contents)
    fh.close()

def out_format(inputname,locus,ref_gene_pos,out_locus,ingene):
    format_out=''
    sumnum=0
    ingene_num=0
    gene_cov=[]
    gene_id=[]
    with open(ref_gene_pos,'r') as handle:
        for line in handle.readlines():
            line=line.strip('\n')
            cut=line.split('\t')
            if cut[1] == locus:sumnum+=1
    if ingene !="None":
        ingene_num=len(ingene)
        for ss in ingene:
            cov_t=float(ss[1])
            gene_cov.append(cov_t)
            gene_id.append(float(ss[2]))
    min_gene_cov=''
    min_gene_id=''
    if ingene_num>0:
        min_gene_cov=str(min(gene_cov))
        min_gene_id=str(min(gene_id))
        #print(min_gene_cov,min_gene_id)
    else:
        min_gene_cov='NA'
        min_gene_id='NA'
    if sumnum>0:
        fre=int(ingene_num)*100/int(sumnum)
        fre='%.1f' % fre+"%"
    else:
        fre=''
    Gene_tag=str(ingene_num)+"/"+str(sumnum)+","+fre
    cut=out_locus.split('\n')
    header="Input file\t"+cut[0]+"\tMatch Gene\tMin Gene Cov(%)\tMin Gene Identity(%)\n";
    context=inputname+"\t"+cut[1]+"\t"+Gene_tag+"\t"+min_gene_cov+"\t"+min_gene_id
    format_out=  header+context+"\n"
    return format_out
    
def DNA_complement(sequence):
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    return sequence.upper()
def DNA_reverse(sequence):
    sequence = sequence.upper()
    return sequence[::-1]

def contains(a, other):
    contained = False
    other_start, other_end = other
    for this_range in a:
        start,end = this_range
        if (int(other_start) >= int(start) and int(other_end) <= int(end)):
            contained = True
            break
    return contained 



def combined(ranges):
    starts_ends = [(x[0], 1) for x in ranges]
    starts_ends += [(x[1], -1) for x in ranges]
    starts_ends.sort(key=lambda z: int(z[0]))
    current_sum = 0
    cumulative_sum = []
    for start_end in starts_ends:
        current_sum += start_end[1]
        cumulative_sum.append((start_end[0], current_sum))
    prev_depth = 0
    start = 0
    combined = []
    for pos, depth in cumulative_sum:
        if prev_depth == 0:
            start = pos
        elif depth == 0:
            combined.append((start, pos))
        prev_depth = depth
    return combined

def coverage(ranges,length):
    ranges=combined(ranges)
    total_len=sum([int(x[1])-int(x[0])+1 for x in ranges])
    cov=100.0 *total_len/float(length)
    return '%.2f' % cov

def mean_identity(ranges,identity):
    identity_sum = 0.0
    length_sum = 0
    for i in range(len(ranges)):
         start,end=ranges[i]
         length_sum+=int(end)-int(start)+1
         identity_sum+=(int(end)-int(start)+1)*float(identity[i])
    if identity_sum == 0.0:
         return 0.0
    else:
        id=identity_sum / length_sum
        return '%.2f' % id

class Hits(object):
    def __init__(self,line):
        cut=line.split("\t")
        self.query=cut[0]
        self.subject=cut[1]
        self.ranges=(cut[2],cut[3])
        self.subject_ranges=cut[4]+".."+cut[5]
        self.length=cut[10]
        self.q_length=cut[11]
        self.identity=cut[9]
        self.blasthit=line
class CleanHit(object):
    def __init__(self,line):
        cut=line.split("\t")
        self.query=cut[0]
        self.subject=cut[1]
        self.q_ranges=(cut[2],cut[3])
        self.s_ranges=(cut[4],cut[5])
        self.s_strand=cut[6]
        self.coverage=cut[7]
        self.identity=cut[9]

    def printout(self):
        return self.query+"\t"+self.subject+"("+self.s_ranges[0]+":"+self.s_ranges[1]+")\t"+self.s_strand+"Coverage:"+self.coverage+"(%)\t"+self.identity+"(%)\n"

def format_blast(blastout):
    out={}
    ranges={}
    identity={}
    q_len={}
    allhits={}
    with open (blastout,'r') as handle:
        for line in handle:
            line=line.strip()
            cut=line.split('\t')
            strand='+'
            match_len=int(cut[3])-int(cut[2])+1
            s_start=cut[4]
            s_end=cut[5]
            if int(cut[4])>int(cut[5]):
                strand='-'
                s_start=cut[5]
                s_end=cut[4]
            subject_len=int(s_end)-int(s_start)+1
            cov_q=100.0* match_len/int(cut[10])
            cov_q="%.2f" %cov_q
            cov_s=100.0* match_len/subject_len
            if match_len>subject_len:
                cov_s=100.0
            cov_s="%.2f" %cov_s
            iden="%.2f" % float(cut[9])
            line=cut[0]+"\t"+cut[1]+"\t"+cut[2]+"\t"+cut[3]+"\t"+s_start+"\t"+s_end+"\t"+strand+"\t"+str(cov_q)+"\t"+str(cov_s)+"\t"+iden+"\t"+cut[8]+"\t"+cut[10]
            out.setdefault(cut[0],[]).append(line)
            temp_Hit=Hits(line)
            allhits.setdefault(cut[0],[]).append(temp_Hit)
            q_len[cut[0]]=cut[10]
    
    return (out,allhits,q_len)


def get_best_k_sign(best_k_ref,cov,identity):
    if cov==100 and identity==100:
        return ("","perfect")
    if cov>=99 and identity>=95:
        return ("*","High")
    else:
        return ("?","Low")

def get_best_k_match(allhits,q_len):
    best_k_ref=''
    best_cov = 0.0
    best_id = 0.0
    best_blast={}
    best_info={}
    for key in allhits:
        fixed_ranges=[]
        fixed_iden=[]
        for hit in allhits[key]:
            range_t=hit.ranges
            if not contains(fixed_ranges, range_t):
                fixed_ranges.append(range_t)
                fixed_iden.append(hit.identity)
        length=q_len[key]             
        cov=float(coverage(fixed_ranges,length))
        mean_id=float(mean_identity(fixed_ranges,fixed_iden))
        best_info.setdefault(key,[]).append(fixed_ranges)
        best_info.setdefault(key,[]).append(fixed_iden)
        if cov > best_cov:
            best_cov = cov
            best_k_ref = key
            best_id =mean_id
        elif cov == best_cov and best_k_ref and mean_id>best_id:
            best_k_ref = key
            best_id=mean_id
    sign1,sign2=get_best_k_sign(best_k_ref,best_cov,best_id)
    for hit in allhits[best_k_ref]:
        for i in  range(len(best_info[best_k_ref][0])):
            if hit.__dict__['ranges']==best_info[best_k_ref][0][i] and hit.__dict__['identity']==best_info[best_k_ref][1][i]:
                temp_hit=CleanHit(hit.__dict__['blasthit'])
                best_blast.setdefault(best_k_ref,[]).append(temp_hit)
                break        
    
    out="Locus_name\tSign\tCoverage(%)\tIdentity(%)\n"+best_k_ref+sign1+"\t"+sign2+"\t"+str(best_cov)+"\t"+str(best_id)+"\n"
    return (best_k_ref,sign1,sign2,out,best_blast)

def get_gene_match(allhits_gene,q_gene_len,locus):
    in_gene=[]
    ingene_cov=[]
    other_gene=[]
    other_gene_temp=[]
    all_ingene=[]
    tag_ingene={}
    for gene in allhits_gene:
        fixed_ranges=[]
        fixed_iden=[]
        for hit in allhits_gene[gene]:
            range_t=hit.ranges
            if not contains(fixed_ranges, range_t):
                fixed_ranges.append(range_t)
                fixed_iden.append(hit.identity)
        length=q_gene_len[gene]
        cov=float(coverage(fixed_ranges,length))
        mean_id=float(mean_identity(fixed_ranges,fixed_iden))
        key=gene.split("_")[0]
        if key == locus:
            all_ingene.append((gene,cov,mean_id))
            genename=gene.split('_')[-1]
            tag_ingene[genename]=1
            if cov >= args.min_gene_cov and float(mean_id) >=args.min_gene_id:
                ingene_cov.append((gene,cov,mean_id))
                for hit in allhits_gene[gene]:
                    tag={}
                    for i in range(len(fixed_ranges)):
                        if hit.__dict__['ranges']==fixed_ranges[i]:
                            rang_t=hit.ranges[0]+".."+hit.ranges[1]
                            name_t=gene+rang_t
                            if name_t not in tag.keys():
                                tag[name_t]=1
                                in_gene.append((gene,key,hit.identity,rang_t,hit.q_length,hit.subject,hit.subject_ranges))
        else:
            if cov >= args.min_gene_cov and float(mean_id) >=args.min_gene_id:
                for hit in allhits_gene[gene]:
                    tag={}
                    for i in range(len(fixed_ranges)):
                        if hit.__dict__['ranges']==fixed_ranges[i]:
                            rang_t=hit.ranges[0]+".."+hit.ranges[1]
                            name_t=gene+rang_t
                            if name_t not in tag.keys():
                                tag[name_t]=1
                                other_gene_temp.append((gene,key,hit.identity,rang_t,hit.q_length,hit.subject,hit.subject_ranges))
    for ss in other_gene_temp:
        genename=ss[0].split('_')[-1]
        if genename not in tag_ingene.keys():
            other_gene.append((ss))
    if not all_ingene:all_ingene="None"
    if not other_gene:other_gene="None"
    if not in_gene:in_gene ="None"
    return(in_gene,ingene_cov,other_gene,all_ingene)

def get_fasta(seq,start,end,strand):
   start=int(start)
   start=start-1
   end=int(end)
   if strand == "+":
       return seq[start:end]
   if strand == "-":
       return DNA_reverse(DNA_complement(seq[start:end]))

def output_locus_match_input(best_k_ref,best_blast,assembly_fasta):
    fasta={}
    seq=""
    out_Fragment="Sequence ID\tLocus Name\tStart\tEnd\tstrand\tCoverage(%)\tIdentity(%)\n"
    for seq_record in SeqIO.parse(assembly_fasta, "fasta"):
        fasta[seq_record.id]=str(seq_record.seq)
    for hit in best_blast[best_k_ref]:
        start,end=hit.s_ranges
        seq+=">"+hit.subject+"|"+start+":"+end+"|"+hit.s_strand+"\n"+get_fasta(fasta[hit.subject],start,end,hit.s_strand)+"\n"
        out_Fragment+=hit.subject+"\t"+best_k_ref+"\t"+start+"\t"+end+"\t"+hit.s_strand+"\t"+hit.coverage+"\t"+hit.identity+"\n"
    return (seq,out_Fragment)

def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    

if __name__ == '__main__':
    parser=ArgumentParser(description="Author:liangqian at 20190806 ;This program is the  pipeline of LPS or CPS typing for Acinetobacter baumannii",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str,help='FASTA file for genomes ')
    parser.add_argument('-d', '--indir', type=str,help='input dirname for multi FASTA files')
    parser.add_argument('-r', '--reference', type=str,required=True,help='Prefix reference loci file in reference_data directory')
    parser.add_argument('-o', '--outdir', type=str, required=False, default='./result',help='output result directory')
    parser.add_argument('-m', '--tempdir', type=str, required=False, default='./tempdir',help='output tempfiles directory')
    parser.add_argument('--min_locus_cov', type=float, required=False, default=95.0,help='minimum coverage for locus sequence match')
    parser.add_argument('--min_locus_id', type=float, required=False, default=90.0,help='minimum identity for locus sequence match')
    parser.add_argument('--min_gene_cov', type=float, required=False, default=99.0,help='minimum coverage for locus genes match')
    parser.add_argument('--min_gene_id', type=float, required=False, default=99.0,help='minimum  identity for locus genes match')
    args = parser.parse_args()
    main(args)
