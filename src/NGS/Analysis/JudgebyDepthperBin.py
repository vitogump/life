# -*- coding: UTF-8 -*-

'''
Created on 2014-6-16

@author: liurui
'''

from optparse import OptionParser
import pickle
import re,copy
from NGS.BasicUtil import Util, VCFutil, Caculators
import src.NGS.BasicUtil.DBManager as dbm


primaryID = "chrID"
chrlength="chrlength"
parser = OptionParser()
parser.add_option("-c", "--genomeCoveragefile", dest="genomedepth", help="genomeCoveragefile")
parser.add_option("-p", "--percentageofCovered", dest="percentageofCovered", help="percentageofCovered")
parser.add_option("-a", "--averagedepthThreshold", dest="averagedepthThreshold", help="averagedepthThreshold")


parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
parser.add_option("-s","--slidesize",dest="slidesize",help="default infile2_infile1")#
parser.add_option("-o","--outfile",dest="outfile",help="outfile")#
parser.add_option("-m","--mindepth",dest="mindepth",help="outfile")#
parser.add_option("-n","--speciesnames",dest="speciesnames",action="append",default=[],help="folder name ")#
(options, args) = parser.parse_args()

outfile=open(options.outfile,'w')
outfilewithvalue=open(options.outfile+"_withvalue",'w')
percentage = float(options.percentageofCovered)
averagedepth=int(options.averagedepthThreshold)

chromtable = Util.pekingduckchromtable
windowWidth=int(options.winwidth)
slideSize=int(options.slidesize)
mindepth=int(options.mindepth)
print(percentage)
if __name__ == '__main__':
    
    depthbinmap={}
    mywin = Util.Window()
    dbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    depthfile = Util.GATK_depthfile(options.genomedepth, options.genomedepth + ".index")
    if "" in depthfile.title:
        depthfile.title.remove("")
    print(depthfile.title,len(depthfile.title)-3)
    if options.speciesnames==[]:
        print("chrom","start_pos","end_pos",*depthfile.title[3:],sep="\t",file=outfile)
        print("chrom","start_pos","end_pos",*depthfile.title[3:],sep="\t",file=outfilewithvalue)
        caculator=Caculators.Caculate_depth_judge(len(depthfile.title)-3,windowWidth,mindepth)
    else:
        outtitle=[]
        sampleidxlisttocount={}
        for speciesname in options.speciesnames:
            sampleidxlisttocount[speciesname]=[]
            outtitle.append("Depth_for_"+speciesname.strip())
            for title in depthfile.title:
                if re.search(r""+speciesname,title)!=None:
                    sampleidxlisttocount[speciesname].append(depthfile.title.index(title.strip())-2)
        print("chrom","start_pos","end_pos",*outtitle[:],sep="\t",file=outfile)
        print("chrom","start_pos","end_pos",*outtitle[:],sep="\t",file=outfilewithvalue)
        caculator=Caculators.Caculate_depth_judge(len(depthfile.title)-3,windowWidth,mindepth,options.speciesnames,sampleidxlisttocount)
    
    
    totalChroms = dbtools.operateDB("select", "select count(*) from " + chromtable)[0][0]
    m = depthfile.covfileidx
    m.pop("title")
    m.pop("chromOrder")
    chrmapListOrderByDepthfilepos = sorted(m.items(), key=lambda m:m[1])
    print(chrmapListOrderByDepthfilepos[-1][0])
    print(*chrmapListOrderByDepthfilepos,sep="\n")
    
    for i in range(0, totalChroms, 20):
        currentsql = "select * from " + chromtable + " order by " + primaryID + " limit " + str(i) + ",20"
        result = dbtools.operateDB("select", currentsql)
        for row in result:
            currentchrID = row[0]
            currentchrLen = int(row[1])
            if currentchrID in depthfile.covfileidx:
                depthfile.depthfilefp.seek(depthfile.covfileidx[currentchrID])
                if chrmapListOrderByDepthfilepos[-1][0]==currentchrID:
                    alllines_a_chrom=depthfile.depthfilefp.read().strip()
                else:
                    alllines_a_chrom=depthfile.depthfilefp.read(chrmapListOrderByDepthfilepos[chrmapListOrderByDepthfilepos.index((currentchrID, depthfile.covfileidx[currentchrID])) + 1][1] - 1 - depthfile.covfileidx[currentchrID]).strip()
                alllist_a_chrom=re.split(r'\n',alllines_a_chrom)
                for lineNo in range(len(alllist_a_chrom)):
                    linelist = re.split(r"\s+", alllist_a_chrom[lineNo].strip())
                    pos = int(re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2))
                    alllist_a_chrom[lineNo]=[pos]+linelist[3:]
                    
                mywin.slidWindowOverlap(alllist_a_chrom, currentchrLen, windowWidth, slideSize, caculator)
                depthbinmap[currentchrID] = copy.deepcopy(mywin.winValueL)
            else:
                fillNA = [(0, 0, 0,([0]*(len(depthfile.title)-3),[0]*(len(depthfile.title)-3)))]
                for i in range(int((currentchrLen - windowWidth) / slideSize)):
                    fillNA.append((0, 0, 0,([0]*(len(depthfile.title)-3),[0]*(len(depthfile.title)-3))))
                depthbinmap[currentchrID] = fillNA
            #print all wins of a chrom into file
            """
            win=[startpos,endpos,([sample1_percentage_cover,sample2_percentage_cover,,,],[sample1_average_depth,sample2_average_depth,,,])]
            or 
            win=[startpos,endpos,([species1_percentage_cover,species2_percentage_cover,,,,,],[species1_average_depth,species2_average_depth,,,,])]
            """
            for winNo in range(len(depthbinmap[currentchrID])):
                judgelist=[]
                valuelist=[]
                for i in range(len(depthbinmap[currentchrID][winNo][3][0])):
                    valuelist.append(str(depthbinmap[currentchrID][winNo][3][0][i])+","+str(depthbinmap[currentchrID][winNo][3][1][i]))
                    if depthbinmap[currentchrID][winNo][3][0][i]>=percentage and int(depthbinmap[currentchrID][winNo][3][1][i]) >=averagedepth:
                        judgelist.append("passed")
                    else:
                        judgelist.append("filtered")
                print(currentchrID,winNo,depthbinmap[currentchrID][winNo][0],depthbinmap[currentchrID][winNo][1],*judgelist,sep="\t",file=outfile)
                print(currentchrID,winNo,depthbinmap[currentchrID][winNo][0],depthbinmap[currentchrID][winNo][1],*valuelist,sep="\t",file=outfilewithvalue)
            depthbinmap.clear()
    depthfile.closedepthfile()
    outfile.close()
    outfilewithvalue.close()
#     species_idx = depthfile.title.index("Depth_for_" + options.species)
