# -*- coding: UTF-8 -*-
from src.NGS.BasicUtil import *
from itertools import combinations

import numpy,pickle,re,sys,copy
from optparse import OptionParser
import src.NGS.BasicUtil.DBManager as dbm
import time
SLEEP_FOR_NEXT_TRY=10

'''
Created on 2013-6-30

@author: rui
'''
parser = OptionParser()
parser.add_option("-c", "--chromtable", dest="chromtable",default="1",help="1 pekingduck,2 D2Bduck")
parser.add_option("-v","--vcffile",dest="vcffile_withdepth",action="append",nargs=2,help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-t","--fsttype",dest="fsttype",help="R(r)/G(g)")
parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
parser.add_option("-s","--slideSize",dest="slideSize",help="default infile2_infile1")#
parser.add_option("-m","--minlength",dest="minlength")
parser.add_option("-j","--outjoin_innerjoin",dest="outjoin_innerjoin",default="o")
parser.add_option("-o","--outputpath",dest="outputpath")
# parser.add_option("-d","--mindepth",dest="mindepth",help="mindepth to judge fixed")
parser.add_option("-F","--considerFixdifferentinFSTcaculate",dest="considerFixdifferentinFSTcaculate",action="store_true",default=False)
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()
print(options.considerFixdifferentinFSTcaculate)

# if len(sys.argv) < 7:
#     print("python CaculateFst.py [vcf1] [vcf2] [vcf3]....[globe_Fst(G)/reletivepaire_Fsts(R)] [winwidth] [slidesize] [chromtable]")
#     exit(-1)

outputpath=options.outputpath.strip()
minlength=options.minlength
windowWidth=int(options.winwidth)
slideSize=int(options.slideSize)
tttttt=False  
if options.chromtable=="1":
    chromtable = Util.pekingduckchromtable
elif options.chromtable=="2":
    chromtable = Util.D2Bduckchromtable
fsttype=options.fsttype
primaryID = "chrID"

vcffileslist=[];depthfilenames={}#{ vcftablename1:[depthfilename1,name1,name2] , vcftablename2:[depthfilename2,name1,name2] } or {vcftablename1:None, vcftablename2:None}
for vcffilename,namefile in options.vcffile_withdepth:
    vcffileslist.append(vcffilename)
    
    if namefile.lower()!="none":
        depthfilenames[vcffilename]=[]
        fp=open(namefile,'r')
        for line in fp:
            depthfile_obj=re.search(r"depthfilename=(.*)",line.strip())
            if depthfile_obj!=None:
                depthfilenames[vcffilename].append(depthfile_obj.group(1).strip())
            elif line.split():
                depthfilenames[vcffilename].append(line.strip())
        fp.close()
    else:
        pass
#         depthfilenames[vcffilename]=None
sql = "select * from " + chromtable

class Fst():
    def __init__(self,vcffileslist,dbtoolsforgenome,depthfilenames):
        super().__init__()
#        self.doubleVcfMap = {}
#         self.FstMapByChrom = {}  # {chr:[(first_snp_pos,last_snp_pos,fst),(),()],chr:[],chr:[]}
        self.distMap = {}
        self.mode="g"
        self.vcffileslist=copy.deepcopy(vcffileslist)
        self.depthfilenames=copy.deepcopy(depthfilenames)
        self.dbtools=dbtoolsforgenome
    def setMode(self,mode,i_o,outputpath,windowWidth,slideSize,chromtable,minlength,tempdpip=None,tempdbusername=None,tempdbpw=None,tempdpname=None):
        self.mode=mode
        self.chromtable=chromtable
        self.minlength=minlength
        self.windowWidth=windowWidth
        self.slideSize=slideSize
        self.outputpath=outputpath
        self.i_o=i_o
        if self.mode=='R' or self.mode=='r':
            self.allspeices=[]
            self.treearrayprename=""
            for pathtoname in self.vcffileslist[:]:
                self.allspeices.append(re.search(r"[^/]*$",pathtoname).group(0).replace('.','_'))
                self.treearrayprename+=re.search(r"[^/]*$",pathtoname).group(0)[0]
            self.phyliparrayinfile=open(outputpath+self.treearrayprename+"phylip.arrayin"+str(self.windowWidth)+"_"+str(self.slideSize),'w')
            print("mysqltablename: "+self.treearrayprename+"treearray")
            arraytitle=""
            for name in self.allspeices:
                arraytitle+=(name+"\t")
            print("\t"+arraytitle+"\n")
            for namerow in self.allspeices:
                print(namerow[0:8]+"\n")        
                
            self.allkindofpaire = list(combinations(self.vcffileslist[:], 2))
            
            self.tempdbtools = dbm.DBTools(tempdpip, tempdbusername, tempdbpw, tempdpname)
            TABLES = {}
            TABLES[self.treearrayprename+"treearray"] = (
                "CREATE TABLE "+self.treearrayprename+"treearray ("
                " `chrID` varchar(128) NOT NULL ,"
                " `winNo` int(18) NOT NULL,"
                " PRIMARY KEY (`chrID`,`winNo`)"
                ")engine=innodb default charset=utf8"
                )
            self.tempdbtools.drop_table(self.treearrayprename+"treearray")
            time.sleep(SLEEP_FOR_NEXT_TRY)
            self.tempdbtools.create_table(TABLES)
        elif self.mode == 'G' or self.mode == 'g':
            self.globalFstMapByChrom={}        
            self.specisnum=str(len(self.vcffileslist[:]))
            
    def executeCaculate(self):
        if self.mode=='R' or self.mode=='r':
            alldistMap={}
            tableindextoarrayindex=[]
            for fstpaire in self.allkindofpaire:
    
                fstpaire1name = re.search(r"[^/]*$",fstpaire[0]).group(0).replace('.','_')
                fstpaire2name = re.search(r"[^/]*$", fstpaire[1]).group(0).replace('.','_')  # for linux
                tableindextoarrayindex.append((self.allspeices.index(fstpaire1name),self.allspeices.index(fstpaire2name)))
                outfilefixdif = open(self.outputpath+fstpaire1name + fstpaire2name + ".df"+str(self.windowWidth)+"_"+str(self.slideSize), 'w')
                outfile = open(self.outputpath+fstpaire1name + fstpaire2name + ".fst"+str(self.windowWidth)+"_"+str(self.slideSize), 'w')
                print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile)
                print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfilefixdif)
                if re.search(r"indvd[^/]*$",fstpaire[0])!=None:
                    MethodToSeqpop1="indvd"
                elif re.search(r"pool[^/]+",fstpaire[0])!=None:
                    MethodToSeqpop1="pool"
                if re.search(r"indvd[^/]+",fstpaire[1])!=None:
                    MethodToSeqpop2="indvd"
                elif re.search(r"pool[^/]+",fstpaire[1])!=None:
                    MethodToSeqpop2="pool"
                fst_caculator = Caculators.Caculate_Fst(MethodToSeqpop1=MethodToSeqpop1, MethodToSeqpop2=MethodToSeqpop2)

                if options.considerFixdifferentinFSTcaculate:
                    fst_caculator.considerfixdiffinfst=True
#                 fst = Fst() 
                FstMapByChrom={}
                self.tempdbtools.operateDB("callproc", "mysql_sp_add_column", data=("ninglabvariantdata_tmp", self.treearrayprename+"treearray", (fstpaire1name[0:5]+fstpaire2name[0:5]), "text", "default null"))
                    
                print("startcaculatefst:\n", fstpaire1name,fstpaire[0],'\n', fstpaire2name,fstpaire[1])
                FstMapByChrom=self.caculateFstAccordingdb(FstMapByChrom,dbtools, self.chromtable, fstpaire[0], fstpaire[1], fst_caculator, self.windowWidth,self.slideSize,self.minlength)
    
                winCrossGenome = []
                for chrom in FstMapByChrom.keys():
                    for i in range(len(FstMapByChrom[chrom])):
                        if FstMapByChrom[chrom][i][3] != "NA":
                            winCrossGenome.append(FstMapByChrom[chrom][i][3])
                exception = numpy.mean(winCrossGenome)
                std0 = numpy.std(winCrossGenome, ddof=0)
                std1 = numpy.std(winCrossGenome, ddof=1)
                del winCrossGenome
                
                totalChroms = self.dbtools.operateDB("select","select count(*) from "+self.chromtable+" where chrlength>="+self.minlength)[0][0]
                for i in range(0,totalChroms,20):
                    currentsql=sql+" where chrlength>="+self.minlength+" order by "+primaryID+" limit "+str(i)+",20"
                    result=self.dbtools.operateDB("select",currentsql)
                    for row in result:
                        currentchrID=row[0]
                        currentchrLen=int(row[1])
                        if currentchrID in FstMapByChrom:
                            for i in range(len(FstMapByChrom[currentchrID])):
                                if FstMapByChrom[currentchrID][i][2][1] != "NA" and FstMapByChrom[currentchrID][i][3] != "NA":
#                                     if i+1==len(FstMapByChrom[currentchrID]):
#                                         lastlength=currentchrLen-(self.windowWidth+self.slideSize*(i-1))
#                                         print(currentchrID + "\t" + str(i) + "\t" + str(FstMapByChrom[currentchrID][i][0]) + "\t" + str(FstMapByChrom[currentchrID][i][1]) +"\t" + str(FstMapByChrom[currentchrID][i][2][1]) + "\t" + 'NA' + "\t" + '%.15f'%(FstMapByChrom[currentchrID][i][2][1]/lastlength), file=outfilefixdif)
#                                     else:
                                    print(currentchrID + "\t" + str(i) + "\t" + str(FstMapByChrom[currentchrID][i][0]) + "\t" + str(FstMapByChrom[currentchrID][i][1]) +"\t" + str(FstMapByChrom[currentchrID][i][2][1]) + "\t" + 'NA' + "\t" + '%.15f'%(FstMapByChrom[currentchrID][i][2][1]/self.windowWidth), file=outfilefixdif)
                                else:
                                    print(currentchrID + "\t" + str(i) + "\t" + str(FstMapByChrom[currentchrID][i][0]) + "\t" + str(FstMapByChrom[currentchrID][i][1]) +"\t" + str(FstMapByChrom[currentchrID][i][2][1]) + "\t" + 'NA' + "\t" + 'NA', file=outfilefixdif)
                                if FstMapByChrom[currentchrID][i][3] != "NA":
                                    zFst = (FstMapByChrom[currentchrID][i][3] - exception) / std1
                                    print(currentchrID + "\t" + str(i) + "\t" + str(FstMapByChrom[currentchrID][i][0]) + "\t" + str(FstMapByChrom[currentchrID][i][1])+ "\t" + str(FstMapByChrom[currentchrID][i][2][0]) + "\t" + '%.15f'%(FstMapByChrom[currentchrID][i][3]) + "\t" + '%.12f'%(zFst), file=outfile) 
                                else:
                                    zFst = "NA"
                                    print(currentchrID + "\t" + str(i) + "\t" + str(FstMapByChrom[currentchrID][i][0]) + "\t" + str(FstMapByChrom[currentchrID][i][1])+ "\t" + str(FstMapByChrom[currentchrID][i][2][0]) + "\t" + str(FstMapByChrom[currentchrID][i][3]) + "\t" + str(zFst), file=outfile)                    
    #            for chrom in sorted(fst.FstMapByChrom.keys()):        
                sum = 0
                Number = 0
                for chrom in sorted(FstMapByChrom.keys()):
                    for i in range(len(FstMapByChrom[chrom])):
                        if FstMapByChrom[chrom][i][3] != 'NA':
                            Number += 1
                            sum += FstMapByChrom[chrom][i][3]
    #                     if chrom=="KB820867.1" or tttttt:
    #                         tttttt=True
                        print("insert into "+self.treearrayprename+"treearray(chrID,winNo,"+fstpaire1name[0:5]+fstpaire2name[0:5]+") values(%s,%s,%s) on duplicate key update "+fstpaire1name[0:5]+fstpaire2name[0:5]+" = '"+str(FstMapByChrom[chrom][i][3])+"'",chrom,str(i),str(FstMapByChrom[chrom][i][3]))
                        self.tempdbtools.operateDB("insert","insert into "+self.treearrayprename+"treearray(chrID,winNo,"+fstpaire1name[0:5]+fstpaire2name[0:5]+") values(%s,%s,%s) on duplicate key update "+fstpaire1name[0:5]+fstpaire2name[0:5]+" = '"+str(FstMapByChrom[chrom][i][3])+"'",data=(chrom,str(i),str(FstMapByChrom[chrom][i][3])))
                alldistMap[fstpaire1name+fstpaire2name] = (sum / Number,self.allspeices.index(fstpaire1name))
                outfile.close()
                outfilefixdif.close()
            
            for n in alldistMap.keys():
                print(n + "\t" + str(alldistMap[n]), file=open(self.outputpath+"testdist.txt", 'a'))
            print("    "+str(len(self.allspeices)),file=open(self.outputpath+"testdist.txt", 'a'))
            for row_sname in self.allspeices:
                print(row_sname[0:8],end="  ",file=open(self.outputpath+"testdist.txt", 'a'))
                for col_sname in self.allspeices:
                    
                    if row_sname+col_sname in alldistMap:
                        print(alldistMap[row_sname+col_sname][0],end="\t",file=open(self.outputpath+"testdist.txt", 'a'))
                    elif col_sname+row_sname in alldistMap:
                        print(alldistMap[col_sname+row_sname][0],end="\t",file=open(self.outputpath+"testdist.txt", 'a'))
                    elif row_sname==col_sname:
                        print("0\t",end="\t",file=open(self.outputpath+"testdist.txt", 'a'))
                print("",file=open(self.outputpath+"testdist.txt", 'a'))
            tatalwins = self.tempdbtools.operateDB("select", "select count(*) from "+self.treearrayprename+"treearray")[0][0]
            for j in range(0, tatalwins, 100):
                wins = self.tempdbtools.operateDB("select","select * from "+self.treearrayprename+"treearray order by chrID asc,winNo asc limit "+str(j) +",100")
                for win in wins:
                    abandonthisWin=False
                    tmparray=[[0 for x in range(len(self.allspeices))] for y in range(len(self.allspeices))]
                    
    #                print("\t"+arraytitle,file=phyliparrayinfile)
                    for i in range(len(win[2:])):
                        tmparray[tableindextoarrayindex[i][0]][tableindextoarrayindex[i][1]]=str(win[i+2])
                        tmparray[tableindextoarrayindex[i][1]][tableindextoarrayindex[i][0]]=str(win[i+2])
                        if win[i+2]==None or win[i+2]=="NA" or win[i+2]=='NULL':
                            abandonthisWin=True
                    if abandonthisWin:
                        continue
                    print("    "+str(len(self.allspeices)),file=self.phyliparrayinfile)
                    for i in range(len(self.allspeices)):
                        tmparray[i][i]='0'
                        try:
                            print(self.allspeices[i][0:8]+"  "+"\t".join(tmparray[i]),file=self.phyliparrayinfile)
                        except TypeError:
                            print(i,self.allspeices,tmparray)
                            print('Error:when making phyliparrayinfile')
                            exit(-1)
            self.tempdbtools.disconnect()
            self.phyliparrayinfile.close()
        elif self.mode=="G" or self.mode=="g":
            for majorpop in self.vcffileslist[:]:
                MethodToSeqpop1=None
                fstlist=[]
                vcfname=re.search(r"[^/]*$",majorpop).group(0)
                vcfname+="VS"
                for othrpop in self.vcffileslist[:]:
                    MethodToSeqpop2=None
                    if majorpop == othrpop:
                        continue
    
                    print("startcaculatefst", majorpop, othrpop)
#                     fstlist.append(Fst())
                    if re.search(r"indvd[^/]+",majorpop)!=None:
                        MethodToSeqpop1="indvd"
                    elif re.search(r"pool[^/]+",majorpop)!=None:
                        MethodToSeqpop1="pool"
                    if re.search(r"indvd[^/]+",othrpop)!=None:
                        MethodToSeqpop2="indvd"
                    elif re.search(r"pool[^/]+",othrpop)!=None:
                        MethodToSeqpop2="pool"
                    FstMapByChrom={}
                    fst_caculator = Caculators.Caculate_Fst(MethodToSeqpop1=MethodToSeqpop1, MethodToSeqpop2=MethodToSeqpop2)

                    if options.considerFixdifferentinFSTcaculate:
                        fst_caculator.considerfixdiffinfst=True
                    fstlist.append(self.caculateFstAccordingdb(FstMapByChrom,self.dbtools, self.chromtable, majorpop, othrpop, fst_caculator, self.windowWidth,self.slideSize,self.minlength))    
                    print("fstlist[-1]",fstlist[-1])
                    vcfname+=("_"+re.search(r"[^/]*$",othrpop).group(0)[0])
                outfile=open(self.outputpath+vcfname+'.gfst'+str(self.windowWidth)+"_"+str(self.slideSize)+"_"+self.specisnum,'w')
                outfilefixdif=open(self.outputpath+vcfname+'.df'+str(self.windowWidth)+"_"+str(self.slideSize)+"_"+self.specisnum,'w')
                print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile)
                print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfilefixdif)
                if len(fstlist) != 0:
                    for chrom in fstlist[0].keys():
                        self.globalFstMapByChrom[chrom]=[]
                        for winNo in range(0,len(fstlist[0][chrom])):
                            sumFstInAWin=0
                            Number=0
                            minNumberOfsnp=10000000000000000000000000000000000000000000
                            for i in range(0,len(fstlist)):
                                try:
    
                                    if fstlist[i][chrom][winNo][3]!= 'NA':
                                        Number+=1
                                        minNumberOfsnp=min(minNumberOfsnp,fstlist[i][chrom][winNo][2][0])
                                        fstlist[i][chrom][winNo][2][0]=minNumberOfsnp
                                        sumFstInAWin+=fstlist[i][chrom][winNo][3]
                                except IndexError:
                                    for j in range(0,len(fstlist)):
                                        print(str(i),str(j),self.vcffileslist[j],chrom,str(winNo),str(len(fstlist[j][chrom])),fstlist[i][chrom][winNo-1])
                                    continue# always in the last position,and the value is caculate any way,so can't mispostion.
                            try:
                                gfst=sumFstInAWin/Number
                            except ZeroDivisionError:
                                gfst="NA"
                            self.globalFstMapByChrom[chrom].append((fstlist[0][chrom][winNo][0],fstlist[0][chrom][winNo][1],copy.deepcopy(fstlist[i][chrom][winNo][2]),gfst))
    
                    winCrossGenome = []
                    for chrom in self.globalFstMapByChrom.keys():
                        for i in range(len(self.globalFstMapByChrom[chrom])):
                            if self.globalFstMapByChrom[chrom][i][3] != "NA":
                                winCrossGenome.append(self.globalFstMapByChrom[chrom][i][3])
                    exception = numpy.mean(winCrossGenome)
                    std0 = numpy.std(winCrossGenome, ddof=0)
                    std1 = numpy.std(winCrossGenome, ddof=1)
                    del winCrossGenome
    
                    totalChroms = dbtools.operateDB("select","select count(*) from "+self.chromtable+" where chrlength>="+self.minlength)[0][0]
                    for i in range(0,totalChroms,20):
                        currentsql=sql+" where chrlength>="+self.minlength+" order by "+primaryID+" limit "+str(i)+",20"
                        result=dbtools.operateDB("select",currentsql)
                        for row in result:
                            currentchrID=row[0]
                            currentchrLen=int(row[1])
                            if currentchrID in self.globalFstMapByChrom:                                
    #                for chrom in sorted(globalFstMapByChrom.keys()):
                                for i in range(len(self.globalFstMapByChrom[currentchrID])):
                                    if self.globalFstMapByChrom[currentchrID][i][2][1] != "NA" and self.globalFstMapByChrom[currentchrID][i][3] != "NA":
#                                         if i+1==len(self.globalFstMapByChrom[currentchrID]):
#                                             lastlength=currentchrLen-(self.windowWidth+self.slideSize*(i-1))
#                                             print(currentchrID + "\t" + str(i) + "\t" + str(self.globalFstMapByChrom[currentchrID][i][0]) + "\t" + str(self.globalFstMapByChrom[currentchrID][i][1]) +"\t" + str(self.globalFstMapByChrom[currentchrID][i][2][1]) + "\t" + 'NA' + "\t" + '%.15f'%(self.globalFstMapByChrom[currentchrID][i][2][1]/lastlength), file=outfilefixdif)
#                                         else:
                                        print(currentchrID + "\t" + str(i) + "\t" + str(self.globalFstMapByChrom[currentchrID][i][0]) + "\t" + str(self.globalFstMapByChrom[currentchrID][i][1]) +"\t" + str(self.globalFstMapByChrom[currentchrID][i][2][1]) + "\t" + 'NA' + "\t" + '%.15f'%(self.globalFstMapByChrom[currentchrID][i][2][1]/self.windowWidth), file=outfilefixdif)
                                    else:
                                        print(currentchrID + "\t" + str(i) + "\t" + str(self.globalFstMapByChrom[currentchrID][i][0]) + "\t" + str(self.globalFstMapByChrom[currentchrID][i][1]) +"\t" + str(self.globalFstMapByChrom[currentchrID][i][2][1]) + "\t" + 'NA' + "\t" + 'NA', file=outfilefixdif)
                                    if self.globalFstMapByChrom[currentchrID][i][3] != "NA":
                                        zgFst = (self.globalFstMapByChrom[currentchrID][i][3] - exception) / std1
    #                                     print(globalFstMapByChrom[currentchrID][i][3])
                                        print(currentchrID + "\t" + str(i) + "\t" + str(self.globalFstMapByChrom[currentchrID][i][0]) + "\t" + str(self.globalFstMapByChrom[currentchrID][i][1]) +"\t" + str(self.globalFstMapByChrom[currentchrID][i][2][0]) + "\t" + '%.15f'%(self.globalFstMapByChrom[currentchrID][i][3]) + "\t" + '%.12f'%(zgFst), file=outfile)
                                    else:
                                        zgFst = "NA"
                                        print(currentchrID + "\t" + str(i) + "\t" + str(self.globalFstMapByChrom[currentchrID][i][0]) + "\t" + str(self.globalFstMapByChrom[currentchrID][i][1]) +"\t" + str(self.globalFstMapByChrom[currentchrID][i][2][0]) + "\t" + self.globalFstMapByChrom[currentchrID][i][3] + "\t" + zgFst, file=outfile)
                outfile.close()
                outfilefixdif.close()
                if len(self.vcffileslist)==2:
                    break
    def caculateFstAccordingdb(self,FstMapByChrom,dbtools,chromstable,vcfNAME_POP1,vcfNAME_POP2,fst_caculator,winwidth,slideSize,minlengthOfchrom):
        pop1 = VCFutil.VCF_Data(vcfNAME_POP1)  # new a class
        pop2 = VCFutil.VCF_Data(vcfNAME_POP2)  # new a class
        if fst_caculator.MethodToSeqpop1=="indvd":
            fst_caculator.pop1_indvdsormediandepth=2*(len(pop1.VcfIndexMap["title"])-9)
        elif fst_caculator.MethodToSeqpop1=="pool":
#             tempf=open(self.depthfilenames[vcfNAME_POP1][0]+".sample_summary","r")
#             for templine in tempf:
#                 if 
#             tempf.
            fst_caculator.pop1_indvdsormediandepth=20
        if fst_caculator.MethodToSeqpop2=="indvd":
            fst_caculator.pop2_indvdsormediandepth=2*(len(pop2.VcfIndexMap["title"])-9)
        elif fst_caculator.MethodToSeqpop2=="pool":
            fst_caculator.pop2_indvdsormediandepth=20
#         print("select count(*) from "+chromstable+" where chrlength>="+minlengthOfchrom)
        totalChroms = dbtools.operateDB("select","select count(*) from "+chromstable+" where chrlength>="+minlengthOfchrom)[0][0]
        ########################### caculate Fst across all vcf file and fill in self.FstMapByChrom 
        for i in range(0,totalChroms,20):
            currentsql=sql+" where chrlength>="+minlengthOfchrom+" order by "+primaryID+" limit "+str(i)+",20"
            result=dbtools.operateDB("select",currentsql)

            for row in result:
                currentchrID=row[0]
                currentchrLen=int(row[1])
                if currentchrID in pop1.VcfIndexMap:
                    pop1SeqOfAChr={}
                    pop2SeqOfAChr={}
                    pop1SeqOfAChr[currentchrID]=pop1.getVcfListByChrom(currentchrID)
                    pop2SeqOfAChr[currentchrID]=pop2.getVcfListByChrom(currentchrID)
                    if self.i_o=="o" and self.depthfilenames!={}:
                        depthobjmap={};fst_caculator.species_idx_map={}
                        depthobjmap[vcfNAME_POP1]=Util.GATK_depthfile(self.depthfilenames[vcfNAME_POP1][0],self.depthfilenames[vcfNAME_POP1][0]+".index")
#                         depthobjmap[vcfNAME_POP1].depthfilefp.seek(depthobjmap[vcfNAME_POP1].covfileidx[currentchrID.strip()])
                        
                        fst_caculator.species_idx_map["vcfpop1_ref"]=[]
                        for name in self.depthfilenames[vcfNAME_POP1][1:]:
                            fst_caculator.species_idx_map["vcfpop1_ref"].append(depthobjmap[vcfNAME_POP1].title.index("Depth_for_"+name))
                        #process vcfpop2

                        depthobjmap[vcfNAME_POP2]=Util.GATK_depthfile(self.depthfilenames[vcfNAME_POP2][0],self.depthfilenames[vcfNAME_POP2][0]+".index")
#                         depthobjmap[vcfNAME_POP2].depthfilefp.seek(depthobjmap[vcfNAME_POP2].covfileidx[currentchrID.strip()])

                        fst_caculator.species_idx_map["vcfpop2"]=[]
                        for name in self.depthfilenames[vcfNAME_POP2][1:]:
                            fst_caculator.species_idx_map["vcfpop2"].append(depthobjmap[vcfNAME_POP2].title.index("Depth_for_"+name))
                        fst_caculator.depthobjmap={}    
                        fst_caculator.depthobjmap["vcfpop1_ref"]=depthobjmap[vcfNAME_POP1]
                        fst_caculator.depthobjmap["vcfpop2"]=depthobjmap[vcfNAME_POP2]
                        fst_caculator.currentchrID=currentchrID
#                         fst_caculator.depthforcurrentchrom=depthfileContentForcurrentchrID
                        
                    newFstMapByChrom=self.caculateFst(FstMapByChrom,pop1SeqOfAChr,pop2SeqOfAChr, fst_caculator,currentchrID,currentchrLen,winwidth,slideSize)
                else:#pop1 don't contation the current chromosome
                    fillNA=[(0,0,[0,'NA'],'NA')]
                    for i in range(int((currentchrLen-winwidth)/slideSize)):
                        fillNA.append((0,0,[0,'NA'],'NA'))
                    newFstMapByChrom[currentchrID]=fillNA
        return copy.deepcopy(newFstMapByChrom)
                                 
    def caculateFst(self, FstMapByChrom,vcfMap1_ref, vcfMap2, caculator,currentchrID,currentchrLen, winwidth, slideSize):
        win = Util.Window()
        try:
            if self.i_o=="i":
                doubleVcfMap=Util.alinmultPopSnpPos([vcfMap1_ref, vcfMap2],"i")
            elif self.i_o=="o":
                doubleVcfMap = Util.alinmultPopSnpPos([vcfMap1_ref, vcfMap2],"o")#produce self.doubleVcfMap{}
#            for currentChrom in self.doubleVcfMap.keys():
    #             self.FstMapByChrom[currentChrom]=[]
            win.winValueL = []
            print("after alinmultPopSnpPos ,caculateFst value in "+currentchrID,len(doubleVcfMap[currentchrID]))
                
            win.slidWindowOverlap(doubleVcfMap[currentchrID], currentchrLen,winwidth, slideSize, caculator)
            FstMapByChrom[currentchrID] = copy.deepcopy(win.winValueL)
        except TypeError:#vcfMap2(pop2) don't contation the current chromosome
            print("caculateFst TypeError")
            fillNA=[(0,0,[max(len(vcfMap1_ref[currentchrID]),len(vcfMap2[currentchrID])),'NA'],'NA')]
            for i in range(int((currentchrLen-winwidth)/slideSize)):
                fillNA.append((0,0,[max(len(vcfMap1_ref[currentchrID]),len(vcfMap2[currentchrID])),'NA'],'NA'))
            FstMapByChrom[currentchrID]=fillNA
        return copy.deepcopy(FstMapByChrom)         

if __name__ == '__main__':
    dbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)

    if fsttype=='R' or fsttype=='r':
        fst=Fst(vcffileslist,dbtools,depthfilenames)
        fst.setMode(fsttype,options.outjoin_innerjoin, outputpath,windowWidth,slideSize,chromtable,minlength,tempdpip=Util.ip,tempdbusername=Util.username,tempdbpw=Util.password,tempdpname=Util.ghostdbname)
        fst.executeCaculate()
        
        
    elif fsttype == 'G' or fsttype == 'g':
        fst=Fst(vcffileslist,dbtools,depthfilenames)
        fst.setMode(fsttype, options.outjoin_innerjoin,outputpath, windowWidth, slideSize,chromtable,minlength)
        fst.executeCaculate()

    dbtools.disconnect()
    print("finsihed")
