#!/usr/bin/python
#-*- conding: UTF-8 -*-
#author: liangqian at 20190805

import sys
import os
import svgwrite
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM

def plot_gene_coverage(locus,geneinfo,allingene,matchGenelist,output):
    geneinfos={}
    locus_len=0
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
        match_gene={} 
        for ss in matchGenelist:
            name,cov,id=ss
            match_gene[name]=[cov,id]
   
        LOWgenelist={}
        for ss in allingene:
            name,cov,id=ss
            if name not in match_gene.keys():
                LOWgenelist[name]=[cov,id]
    
        edgeTop,edgeBottom,edgeLeft,edgeRight = (150,50,20,140);
        W,H=(1200,200)
        width=W+edgeLeft+edgeRight
        height=H+edgeTop+edgeBottom
        dwg = svgwrite.Drawing(output+"/locus_gene.coverage.svg",size=(width,height))
        x1=edgeLeft
        x2=W+edgeLeft;y1=y2=H+edgeTop
        dwg.add(dwg.line((x1, y1), (x2, y2), stroke="black",stroke_width=2))
        y2=H+edgeTop+8
        dwg.add(dwg.line((x1, y1), (x1, y2), stroke="black",stroke_width=2))
        dwg.add(dwg.line((x2, y1), (x2, y2), stroke="black",stroke_width=2))
        x1=W+edgeLeft+10;y1=H+edgeTop+3;
        dwg.add(dwg.text(locus, insert=(x1,y1),fill="#a203aa",style="font-family:Times New Roman,Times,serif;font-size:18px;font-weight:bold"))
        x1=edgeLeft-4;y1=y1=H+edgeTop+30;
        dwg.add(dwg.text("1", insert=(x1,y1),fill="rgb(0,0,0)",style="font-family:Times New Roman,Times,serif;font-size:16px;font-weight:bold"))
        x1=W+edgeLeft-20
        y1=H+edgeTop+30;text1=str(locus_len)+"bp"
        dwg.add(dwg.text(text1, insert=(x1,y1),fill="rgb(0,0,0)",style="font-family:Times New Roman,Times,serif;font-size:16px;font-weight:bold"))
        subx=float(W/locus_len)
        center_x=float(W/locus_len)*int(locus_len/2)+edgeLeft;y1=H+edgeTop;y2=H+edgeTop+8
        dwg.add(dwg.line((center_x, y1), (center_x, y2), stroke="black",stroke_width=2))
        x1=center_x-10;y1=H+edgeTop+30;text2=str(int(locus_len/2))+"bp"
        dwg.add(dwg.text(text2, insert=(x1,y1),fill="rgb(0,0,0)",style="font-family:Times New Roman,Times,serif;font-size:16px;font-weight:bold"))
        x1=edgeLeft;x2=W+edgeLeft;y1=y2=H+edgeTop-2-25;
        dwg.add(dwg.line((x1, y1), (x2, y2), stroke="black",stroke_width=1))
        for key,value in geneinfos.items():
            ref,pos,strand=value
            start,end=pos.split(':')
            start=int(start);end=int(end)
            length=end-start+1
            x2=edgeLeft+start*subx;x1=edgeLeft+start*subx+length*2*subx/3;x3=edgeLeft+end*subx
            y1=H+edgeTop-2;y2=y1-12;y3=y1-40;y4=y3-12;y5=y3+16;
            if strand == "-":
                x2=edgeLeft+end*subx;x3=edgeLeft+start*subx;x1=edgeLeft+start*subx+length*1*subx/3
            color=""
            stroke_color="#777777"
            cov=0;id=0
            if key in match_gene.keys():
                cov,id=match_gene[key]
                if cov== 100:
                    color="#d00b1e"
                elif cov >99:
                    color="#d00b0e"
                else:
                    color="#d0390b"
            elif key in LOWgenelist.keys():
                color="#ffbb63"
                cov,id=LOWgenelist[key]
            else: 
                color="#ffffff"
            dwg.add(dwg.path(d='M{0},{1} L{2},{3} L{4},{5} L{6},{7} L{8},{9} L{10},{11} L{12},{13} L{14},{15} Z'.format(x1,y1,x1,y2,x2,y2,x2,y3,x1,y3,x1,y4,x3,y5,x1,y1),fill=color, stroke=stroke_color,stroke_width=1))
            if cov !=0 and id !=0:
                x1=edgeLeft+start*subx+20;y1=H+edgeTop-60;text=key+"(Cov:"+str(cov)+"%,Id:"+str(id)+"%)"
                rotate1="rotate"+"(-50,"+str(x1)+","+str(y1)+")"
                dwg.add(dwg.text(text, insert=(x1,y1),fill="#0017ec", transform=rotate1,style="font-size:14px;font-weight:bold"))
        dwg.save()
        #cairosvg.svg2png(url=output+"/locus_gene.coverage.svg", write_to=output+'/locus_gene.coverage.png')
        drawing = svg2rlg(output+"/locus_gene.coverage.svg")
        renderPM.drawToFile(drawing, output+'/locus_gene.coverage.png', fmt="PNG")
    else:
        return("The locus "+locus+" have no reference genes!")
 
    
