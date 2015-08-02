'''
Created on 2015-8-1

@author: liurui
'''
import fractions
from optparse import OptionParser
import re

from NGS.BasicUtil import *
import NGS.BasicUtil.Util
import src.NGS.BasicUtil.DBManager as dbm


primaryID = "chrID"
parser = OptionParser()
# parser.add_option("-c", "--chromtable", dest="chromtable",# action="callback",type="string",callback=useoptionvalue_previous2,
#                   help="write report to FILE")
parser.add_option("-T","--targetpopvcffile_withdepth",dest="targetpopvcffile_withdepth",action="append",nargs=2,help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-R","--refpopvcffile_withdepth",dest="refpopvcffile_withdepth",action="append",nargs=2,help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-t","--topleveltaple",dest="topleveltaple",help="R(r)/G(g)")
parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
parser.add_option("-s","--slideSize",dest="slideSize",help="default infile2_infile1")#
parser.add_option("-m","--minlength",dest="minlength")
parser.add_option("-n","--numberofindvdoftargetpop_todividintobin",dest="numberofindvdoftargetpop_todividintobin",default="o")
parser.add_option("-o","--outfileprewithpath",dest="outfileprewithpath")

parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

windowWidth=int(options.winwidth)
slideSize=int(options.slideSize)
if __name__ == '__main__':
    d_increase=fractions.Fraction(1, (2*int(options.numberofindvdoftargetpop_todividintobin)))
    d_increase=round(d_increase,11)
    intervalFileName=options.outfileprewithpath+".interval"
    intervalfile=open(intervalFileName,"w")
    minvalue=0.000000000000
    freq_fixednumberRelation={}
    
    while minvalue+d_increase<=1:
        freq_fixednumberRelation[(minvalue,minvalue+d_increase+0.00000000004)]=[]
        print('%.12f'%minvalue,'%.12f'%(minvalue+d_increase+0.00000000004),sep="\t",file=intervalfile)
        minvalue+=d_increase

    intervalfile.close()
    outputname=options.outfileprewithpath
    poplist=[];methodlist=[];listofpopvcfmapOfAChr=[];depthobjmap={}
    depth_idxname_mapbyfilenames={}#{ vcftablename1:[depthfilename1,name1,name2] , vcftablename2:[depthfilename2,name1,name2] } or {vcftablename1:None, vcftablename2:None}
    for vcf,depthconfig in options.targetpopvcffile_withdepth[:]+options.refpopvcffile_withdepth[:]:
        listofpopvcfmapOfAChr.append({})
        outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcf).group(0))[0])
        poplist.append(VCFutil.VCF_Data(vcf))  # new a class
        if depthconfig.lower()!="none":
            depth_idxname_mapbyfilenames[vcf]=[]
            fp=open(depthconfig,'r')
            for line in fp:
                depthfile_obj=re.search(r"depthfilename=(.*)",line.strip())
                if depthfile_obj!=None:
                    depth_idxname_mapbyfilenames[vcf].append(depthfile_obj.group(1).strip())#append depthfilename1
                    depthobjmap[vcf]=Util.GATK_depthfile(depth_idxname_mapbyfilenames[vcf][0],depth_idxname_mapbyfilenames[vcf][0]+".index")
                elif line.split():
                    depth_idxname_mapbyfilenames[vcf].append(line.strip())#append idxname
            fp.close()
        if re.search(r"indvd[^/]+",vcf)!=None:
            methodlist.append("indvd")
        elif re.search(r"pool[^/]+",vcf)!=None:
            methodlist.append("pool")
        
    outfile = open(outputname + ".het"+str(windowWidth)+"_"+str(slideSize), 'w')   
    genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    
    
    print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile)
    print("select","select count(*) from "+Util.pekingduckchromtable + " where chrlength>="+options.minlength)
    N_of_targetpop=len(options.targetpopvcffile_withdepth)
    N_of_refpop=len(options.refpopvcffile_withdepth)
 
#make correlation file     
    totalChroms = genomedbtools.operateDB("select","select count(*) from "+Util.pekingduckchromtable + " where chrlength>="+options.minlength)[0][0]
    for i in range(0,totalChroms,20):
        currentsql="select * from " + Util.pekingduckchromtable+" where chrlength>="+options.minlength+" order by "+primaryID+" limit "+str(i)+",20"
        result=genomedbtools.operateDB("select",currentsql)
        
        for row in result:
            currentchrID=row[0].strip()
            currentchrLen=int(row[1])
            for vcfobj in poplist:
                if currentchrID in vcfobj.VcfIndexMap:
                    break
            else:
                continue
            #this chr exist in one of the vcffile
            for vcfobj_idx in range(len(poplist)):
                listofpopvcfmapOfAChr[vcfobj_idx]={}
                listofpopvcfmapOfAChr[vcfobj_idx][currentchrID]=poplist[vcfobj_idx].getVcfListByChrom(currentchrID)
            target_ref_SNPs=Util.alinmultPopSnpPos(listofpopvcfmapOfAChr, "o")
            for snp in target_ref_SNPs[currentchrID]:
                for i in range(3,N_of_targetpop):
                snp=dbvariantstools.operateDB("select","select * from "+options.topleveltablejudgeancestral+" where chrID='"+currentchrID+"' and snp_pos='"+str(curpos)+"'")
    #slide window to caculate S
    for i in range(0,totalChroms,20):
        currentsql="select * from " + Util.pekingduckchromtable+" where chrlength>="+options.minlength+" order by "+primaryID+" limit "+str(i)+",20"
        result=genomedbtools.operateDB("select",currentsql)
        for row in result:
            currentchrID=row[0].strip()
            currentchrLen=int(row[1])
            for vcfobj in poplist:
                if currentchrID in vcfobj.VcfIndexMap:
                    break
            else:
                continue
            #this chr exist in one of the vcffile
            for vcfobj_idx in range(len(poplist)):
                listofpopvcfmapOfAChr[vcfobj_idx]={}
                listofpopvcfmapOfAChr[vcfobj_idx][currentchrID]=poplist[vcfobj_idx].getVcfListByChrom(currentchrID)
                if depthobjmap!={}:
                    depthobjmap[(options.targetpopvcffile_withdepth[:]+options.refpopvcffile_withdepth[:])[vcfobj_idx][0]].depthfilefp.seek(depthobjmap[(options.targetpopvcffile_withdepth[:]+options.refpopvcffile_withdepth[:])[vcfobj_idx][0]].covfileidx[currentchrID])
            target_ref_SNPs=Util      

    
    
    