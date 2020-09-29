#!/usr/bin/python
#-*- conding: UTF-8 -*-
#author: liangqian at 20190905
'''
Author: liangqian at 20190905
This program is ploting genome with locus  typeing result,only for perfect sign*
'''
import sys
import os
import svgwrite
#import cairosvg
import numpy as np
import re
from collections import OrderedDict
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM

def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def sorted_twolist(a,b):
    a=[ int(x) for x in a ]
    b=[ int(x) for x in b ]
    c=[]
    d=[]
    c=np.sort(a)
    x=np.array(a)
    y=x.argsort()
    d=[b[i] for i in y]
    return (c,d)
    
def deal_blast(match_blast,gene_blast):
    genes={}
    genes=OrderedDict()
    loucs=''
    contig=''
    ranges=[]
    strand=[]
    contig_info={}
    start1=[];end1=[];match_info={};
    num=0
    for line in open(match_blast,'r').readlines():
        num=num+1
        if num >1:
            line=line.strip('\n')
            cut=line.split('\t')
            locus=cut[1]
            locus=locus.replace('*','')
            locus=locus.replace('?','')
            contig=cut[0]
            strand.append(cut[4])
            start1.append(cut[2])
            end1.append(cut[3])
    start,end= sorted_twolist(start1,end1)
    #print(start1)
    #print(end1)
    #print(start)
    #print(end)
    for i in range(len(start)):
        ranges.append(str(start[i])+".."+str(end[i]))
    #print(ranges)
    flag={}
    for line in open(gene_blast,'r').readlines():
        line=line.strip('\n')
        cut=line.split('\t')
        key=cut[0].split('_')[0]
        if key == locus and cut[1] == contig:
            st='+'
            s=cut[4]
            e=cut[5]
            if int(cut[4])>int(cut[5]):
                st='-';s=cut[5];e=cut[4]
            contig_start=(ranges[0].split('..'))[0]
            contig_end=(ranges[-1].split('..'))[1]
            if int(s)>=int(contig_start) and int(e)<=int(contig_end) and cut[0] not in flag.keys():
                genes[cut[0]]=cut[2]+"|"+cut[3]+"|"+st+"|"+s+"|"+e
                flag[cut[0]]=1
    #print(genes)
    return(locus,contig,ranges,genes)

def plot_Collinear(locus,locus_lenfile,contig,ranges,geneinfo,matchgenes,output,key_name):
    #print("***********test**********"+locus)
    geneinfos={}
    locus_len=0
    with open(locus_lenfile,'r') as handle:
        for line in handle.readlines():
            line=line.strip('\n')
            cut=line.split('\t')
            if cut[0] == locus:
                locus_len=cut[1]
    have_geneTag='False'
    with open(geneinfo,'r') as handle:
        for line in handle.readlines():
            line=line.strip('\n')
            cut=line.split('\t')
            if cut[1] == locus:
                geneinfos[cut[0]]=[cut[1],cut[2],cut[3]]
                locus_len=int(cut[4])
                have_geneTag='True'
    if have_geneTag == 'True':
        match_len=0
        contig_start=(ranges[0].split('..'))[0]
        contig_end=(ranges[-1].split('..'))[1]
        for i in range(len(ranges)):
            (start1,end1)=ranges[i].split('..')
            match_len+=int(end1)-int(start1)+1
        contig_len=int(contig_end)-int(contig_start)+1
        edgeTop,edgeBottom,edgeLeft,edgeRight = (80,40,160,80);
        W,H=(1200,100)
        width=W+edgeLeft+edgeRight
        height=H+edgeTop+edgeBottom
        scaleX=int(W)/int(locus_len)
        len_tag=abs(contig_len-locus_len)
        if len_tag>5000:
            exit("warnings: sequence length difference >5000")
        if locus_len <=contig_len:
            scaleX=int(W)/int(contig_len)
        dwg = svgwrite.Drawing(output+"/out.Collinear.svg",size=(width,height))
        dele_number=len(geneinfos.keys())-len(matchgenes.keys())
        ret_list = [item for item in geneinfos.keys() if item not in matchgenes.keys()]
        #print(key_name+"\tMissing Genes:\t"+str(dele_number))
        #print(",".join(str(i) for i in ret_list)) 
        out=open(output+"/Missing_genes.txt",'w')
        out.write(key_name+"\t"+contig+"\t"+",".join(str(i) for i in ret_list)+"\n")
        out.close()
        ###**********plot locus and gene info**************************************
        x1=edgeLeft
        x2=scaleX*int(locus_len)+edgeLeft;y1=y2=H+edgeTop
        dwg.add(dwg.line((x1, y1), (x2, y2), stroke="black",stroke_width=2))
        y2=H+edgeTop+8
        dwg.add(dwg.line((x1, y1), (x1, y2), stroke="black",stroke_width=2))
        dwg.add(dwg.line((x2, y1), (x2, y2), stroke="black",stroke_width=2))
        x1=edgeLeft-10;y1=H+edgeTop+3
        dwg.add(dwg.text(locus, insert=(x1,y1),fill="red",style="font-family:Times New Roman,Times,serif;font-size:18px;font-weight:bold;text-anchor:end"))
        x1=edgeLeft-4;y1=y1=H+edgeTop+30;
        dwg.add(dwg.text("1", insert=(x1,y1),fill="rgb(0,0,0)",style="font-family:Times New Roman,Times,serif;font-size:16px;font-weight:bold"))
        x1=scaleX*int(locus_len)+edgeLeft-20
        y1=H+edgeTop+30;text1=str(locus_len)+"bp"
        dwg.add(dwg.text(text1, insert=(x1,y1),fill="rgb(0,0,0)",style="font-family:Times New Roman,Times,serif;font-size:16px;font-weight:bold"))
        center_x=scaleX*int(locus_len/2)+edgeLeft;y1=H+edgeTop;y2=H+edgeTop+8
        dwg.add(dwg.line((center_x, y1), (center_x, y2), stroke="black",stroke_width=2))
        x1=center_x-10;y1=H+edgeTop+30;text2=str(int(locus_len/2))+"bp"
        dwg.add(dwg.text(text2, insert=(x1,y1),fill="rgb(0,0,0)",style="font-family:Times New Roman,Times,serif;font-size:16px;font-weight:bold"))
        x1=edgeLeft;x2=scaleX*int(locus_len)+edgeLeft;y1=y2=H+edgeTop-2-24;
        dwg.add(dwg.line((x1, y1), (x2, y2), stroke="black",stroke_width=1))
        gene_lens={}
        for key,value in geneinfos.items():
            ref,pos,strand=value
            start,end=pos.split(':')
            start=int(start);end=int(end)
            length=end-start+1
            gene_lens[key]=length
            x2=edgeLeft+start*scaleX;x1=edgeLeft+start*scaleX+length*2*scaleX/3;x3=edgeLeft+end*scaleX
            y1=H+edgeTop-2;y2=y1-12;y3=y1-40;y4=y3-12;y5=y3+16;
            if strand == "-":
                x2=edgeLeft+end*scaleX;x3=edgeLeft+start*scaleX;x1=edgeLeft+start*scaleX+length*1*scaleX/3;
            color1="#ffffff"
            if key not in matchgenes.keys():
                color1="#c2c2c2"
            stroke_color="#b1b1b1"
            dwg.add(dwg.path(d='M{0},{1} L{2},{3} L{4},{5} L{6},{7} L{8},{9} L{10},{11} L{12},{13} L{14},{15} Z'.format(x1,y1,x1,y2,x2,y2,x2,y3,x1,y3,x1,y4,x3,y5,x1,y1),fill=color1, stroke=stroke_color))
        ###********************************plot contig infomations and match ranges*******************************
        x1=edgeLeft
        x2=scaleX*int(contig_len)+edgeLeft;y1=edgeTop-30;y2=y1+10;
        dwg.add(dwg.path(d='M{0},{1} L{2},{3} M{4},{5} L{6},{7} '.format(x1,y1,x2,y1,x2,y2,x1,y2),fill="#ffffff", stroke="#000",stroke_width=2))
        w=scaleX*int(contig_len)
        num=int(w/10)
        reverse=0;flag_s0=0;flag_s1=0;
        tempgenes=list(matchgenes.keys())
        flag_s0=matchgenes[tempgenes[0]].split('|')[3];flag_s1=matchgenes[tempgenes[-1]].split('|')[3]
        #print(flag_s0,flag_s1)
        if int(flag_s0)>int(flag_s1):
            reverse=1
        if reverse==1 :
            temp=contig_start
            contig_start=contig_end
            contig_end=temp
            
        for i in range(num):
            x11=edgeLeft+i*10;x21=edgeLeft+i*10+10;
            dwg.add(dwg.line((x21, y1), (x11, y2), stroke="#666",stroke_width=1))
        x11=x1-10;y11=y1+10;
        text=str(key_name)
        dwg.add(dwg.text(text, insert=(x11,y11),fill="red",style="font-family:Times New Roman,Times,serif;font-size:18px;font-weight:bold;text-anchor:end"))
        x11=x1-5;y11=y1-10;
        text1=str(contig_start)
        dwg.add(dwg.text(text1, insert=(x11,y11),fill="rgb(0,0,0)",style="font-family:Times New Roman,Times,serif;font-size:16px;font-weight:bold"))
        x12=x2-10;y12=y11
        text2=str(contig_end)
        dwg.add(dwg.text(text2, insert=(x12,y12),fill="rgb(0,0,0)",style="font-family:Times New Roman,Times,serif;font-size:16px;font-weight:bold"))
        out=open(output+"/contig.info.txt",'w')
        out.write(key_name+"\t"+contig+"\t"+contig_start+"\t"+contig_end+"\n")
        out.close()
        out=open(output+"/match.gene.txt",'w')
        out.write("Gene\tstart\tend\tstrand\tcontig\tcontig_start\tcontig_end\tgene_length\n")
        for gene in matchgenes.keys():
            cut=matchgenes[gene].split('|')
            color2="#f44336"
            ref,pos,strand_t=geneinfos[gene]
            gene_start,gene_end=pos.split(':')
            start=int(cut[0])*3+int(gene_start)-1;end=int(cut[1])*3+int(gene_start)-1;
            if cut[2]=="-":
                color2="#00aa54"
                start=int(gene_end)-int(cut[1])*3+1;end=int(gene_end)-int(cut[0])*3+1;
            #print(gene+str(cut)+str(gene_lens[gene]))
            out.write(gene+"\t"+cut[0]+"\t"+cut[1]+"\t"+cut[2]+"\t"+contig+"\t"+cut[3]+"\t"+cut[4]+"\t"+str(gene_lens[gene])+"\n")
            start2=int(cut[3])-int(contig_start)+1;end2=int(cut[4])-int(contig_start)+1;
            if reverse==1:
                start2=int(contig_start)-int(cut[4])+1;end2=int(contig_start)-int(cut[3])+1;
            x1=edgeLeft+start*scaleX;x2=edgeLeft+end*scaleX
            x3=edgeLeft+start2*scaleX;x4=edgeLeft+end2*scaleX
            y1=H+edgeTop;y2=edgeTop-20
            dwg.add(dwg.path(d='M{0},{1} L{2},{3} L{4},{5} L{6},{7} Z'.format(x1,y1,x2,y1,x4,y2,x3,y2),fill=color2, stroke="none",opacity="0.6"))
        dwg.save()
        out.close()
        #cairosvg.svg2png(url=output+"/out.Collinear.svg", write_to=output+'/out.Collinear.png')
        drawing = svg2rlg(output+"/out.Collinear.svg")
        renderPM.drawToFile(drawing, output+'/out.Collinear.png', fmt="PNG")
    else:
        return("The locus "+locus+" have no reference genes!")
 

if __name__ == '__main__':
    match_blast=sys.argv[1]
    blast_gene=sys.argv[2]
    gene_info=sys.argv[3]
    locus_lenfile=sys.argv[4]
    outdir=sys.argv[5]
    key=sys.argv[6]
    make_dir(outdir)
    locus,contig,ranges,match_genes=deal_blast(match_blast,blast_gene)
    plot_Collinear(locus,locus_lenfile,contig,ranges,gene_info,match_genes,outdir,key)
