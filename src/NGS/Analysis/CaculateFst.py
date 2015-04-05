# -*- coding: UTF-8 -*-
from NGS.BasicUtil import *
from itertools import combinations
import NGS.BasicUtil.Util
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
parser.add_option("-d", "--chromdbname", dest="chromdbname",# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="write report to FILE")
parser.add_option("-c", "--chromtable", dest="chromtable",# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="write report to FILE")
parser.add_option("-v","--vcffile",dest="vcffile",action="append", default=[],help="default infile1_infile2")
parser.add_option("-t","--fsttype",dest="fsttype",help="R(r)/G(g)")
parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
parser.add_option("-s","--slideSize",dest="slideSize",help="default infile2_infile1")#
parser.add_option("-m","--minlength",dest="minlength")
parser.add_option("-o","--outputpath",dest="outputpath")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
# if len(sys.argv) < 7:
#     print("python CaculateFst.py [vcf1] [vcf2] [vcf3]....[globe_Fst(G)/reletivepaire_Fsts(R)] [winwidth] [slidesize] [chromtable]")
#     exit(-1)
outputpath=options.outputpath.strip()
minlength=options.minlength
windowWidth=int(options.winwidth)
slideSize=int(options.slideSize)
chromdbname=options.chromdbname
chromtable = options.chromtable
fsttype=options.fsttype
primaryID = "chrID"
vcffileslist=options.vcffile
sql = "select * from " + chromtable

class Fst():
    def __init__(self):
        super().__init__()
#        self.doubleVcfMap = {}
        self.FstMapByChrom = {}  # {chr:[(first_snp_pos,last_snp_pos,fst),(),()],chr:[],chr:[]}
        self.distMap = {}
#     def alin2PopSnpPos(self,innerjoin_outjoin="i", *vcfMap):
#         """input:
#         two map fomart like this {chrNo:[(pos,REF,ALT,INFO,FORMAT,sample,...),(pos,REF,ALT,INFO,FORMAT,sample,...),,,,,],chrNo:[],,,,,,}
#         output:
#         one map like this {chrNo:[(pos,REF,ALT,(INFO,FORMAT,sample,...),(INFO,FORMAT,sample,...)),(,,,(),()),,,,,],chrNo:[],,,}
#                                                 from pop1                        from pop2
#         """
#         doubleVcfMap={}
#         multipleVcfMap={}
#         for currentChrom in vcfMap[0].keys():
# #             self.FstMapByChrom[currentChrom] = []
#             doubleVcfMap[currentChrom] = []
#             multipleVcfMap[currentChrom]=[]
#             for SNPrec in vcfMap[0][currentChrom]:
#                 posInPop1 = SNPrec[0]
#                 RefInPop1 = SNPrec[1]
#                 AltInPop1 = SNPrec[2]
#                 skipthisrec=False
#                 elementToAppend=[posInPop1,RefInPop1,AltInPop1,SNPrec[3:]]
#                 for vcfMap_obj_idx in range(1,len(vcfMap[:])):
#                     vcfMap_obj=vcfMap[vcfMap_obj_idx]
#                     if currentChrom not in vcfMap_obj:
#                         print("alin2PopSnpPos",currentChrom,"didn't find in vcfMap2")
#                         if innerjoin_outjoin=="i":
#                             skipthisrec=True
#                             break
#                         elif innerjoin_outjoin=="o":
#                             elementToAppend.append(None)
#                     low = 0
#                     high = len(vcfMap_obj[currentChrom]) - 1
#                     
#                     if re.search(r"[A-Za-z]+,[A-Za-z]+", AltInPop1) != None:  # multiple allels
#                         continue
#     #                dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", SNPrec[3])
#     #                 print(dp4.group(0))
#                     
#                     while low <= high:
#                         mid = (low + high)>>1
#                         if vcfMap_obj[currentChrom][mid][0]<posInPop1:
#                             low=mid+1
#                         elif vcfMap_obj[currentChrom][mid][0]>posInPop1:
#                             high=mid-1
#                         else:
#                             if AltInPop1 == vcfMap_obj[currentChrom][mid][2]:#same alt alle
#                                 if vcfMap_obj_idx!=len(vcfMap):
#                                     elementToAppend.append(vcfMap_obj[currentChrom][mid][3:])
#                                 elif vcfMap_obj_idx==len(vcfMap):
#                                     multipleVcfMap[currentChrom].append(elementToAppend)
#                             elif innerjoin_outjoin=="i":
#                                 skipthisrec=True
#                                 print(currentChrom,posInPop1,AltInPop1,vcfMap_obj[currentChrom][mid][2],"different alt allele,should skip this rec,but i have no time to improve this now")
#                             elif innerjoin_outjoin=="o":
#                                 if vcfMap_obj_idx!=len(vcfMap):
#                                     elementToAppend.append(None)
#                                 elif vcfMap_obj_idx==len(vcfMap):
#                                     multipleVcfMap[currentChrom].append(elementToAppend)                                
#                                 
#                             break
#                     else:
#                         if innerjoin_outjoin=="i" and skipthisrec:
#                             #ignore the rec
#                             break
#                         elif innerjoin_outjoin=="o":
#                             if vcfMap_obj_idx!=len(vcfMap):
#                                 elementToAppend.append(None)
#                             elif vcfMap_obj_idx==len(vcfMap):
#                                 multipleVcfMap[currentChrom].append(elementToAppend)                              
# #                     print("snp not found in vcfMap2",SNPrec)
# #                     self.doubleVcfMap[currentChrom].append(SNPrec+)
#         return multipleVcfMap

    def caculateFstAccordingdb(self,dbtools,chromstable,vcfNAME_POP1,vcfNAME_POP2,caculator,winwidth,slideSize,minlengthOfchrom):
        pop1 = VCFutil.VCF_Data(vcfNAME_POP1)  # new a class
        pop2 = VCFutil.VCF_Data(vcfNAME_POP2)  # new a class
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
                    pop1SeqOfAChr[currentchrID]=pop1.getVcfListByChrom(vcfNAME_POP1, currentchrID)
                    pop2SeqOfAChr[currentchrID]=pop2.getVcfListByChrom(vcfNAME_POP2, currentchrID)
                    self.caculateFst(pop1SeqOfAChr,pop2SeqOfAChr, fst_caculator,currentchrID,currentchrLen,winwidth,slideSize)
                else:#pop1 don't contation the current chromosome
                    fillNA=[(0,0,0,'NA')]
                    for i in range(int((currentchrLen-windowWidth)/slideSize)):
                        fillNA.append((0,0,0,'NA'))
                    self.FstMapByChrom[currentchrID]=fillNA
                                 
    def caculateFst(self, vcfMap1_ref, vcfMap2, caculator,currentchrID,currentchrLen, winwidth, slideSize):
        win = Util.Window()
        try:
#            self.doubleVcfMap={}
            doubleVcfMap = Util.alin2PopSnpPos([vcfMap1_ref, vcfMap2],"i")#produce self.doubleVcfMap{}
#            for currentChrom in self.doubleVcfMap.keys():
    #             self.FstMapByChrom[currentChrom]=[]
            win.winValueL = []
            print("after alin2PopSnpPos ,caculateFst value in "+currentchrID)
            
            win.slidWindowOverlap(doubleVcfMap[currentchrID], currentchrLen,winwidth, slideSize, caculator)
            self.FstMapByChrom[currentchrID] = copy.deepcopy(win.winValueL)
        except TypeError:#vcfMap2(pop2) don't contation the current chromosome
            print("caculateFst TypeError")
            fillNA=[(0,0,0,'NA')]
            for i in range(int((currentchrLen-windowWidth)/slideSize)):
                fillNA.append((0,0,0,'NA'))
            self.FstMapByChrom[currentchrID]=fillNA           

if __name__ == '__main__':
    dbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)

    if fsttype=='R' or fsttype=='r':
        
        allspeices=[]
        tableindextoarrayindex=[]
        treearrayprename=""
        for pathtoname in vcffileslist[:]:
            allspeices.append(re.search(r"[^/]*$",pathtoname).group(0).replace('.','_'))
            treearrayprename+=re.search(r"[^/]*$",pathtoname).group(0)[0]
        phyliparrayinfile=open(outputpath+treearrayprename+"phylip.arrayin"+str(windowWidth)+"_"+str(slideSize),'w')
        print("mysqltablename: "+treearrayprename+"treearray")
        arraytitle=""
        for name in allspeices:
            arraytitle+=(name+"\t")
        print("\t"+arraytitle+"\n")
        for namerow in allspeices:
            print(namerow[0:8]+"\n")        
            
        allkindofpaire = list(combinations(vcffileslist[:], 2))
        alldistMap={}
        tempdbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.ghostdbname)
        TABLES = {}
        TABLES[treearrayprename+"treearray"] = (
            "CREATE TABLE "+treearrayprename+"treearray ("
            " `chrID` varchar(128) NOT NULL ,"
            " `winNo` int(18) NOT NULL,"
            " PRIMARY KEY (`chrID`,`winNo`)"
            ")engine=innodb default charset=utf8"
            )
        tempdbtools.drop_table(treearrayprename+"treearray")
        time.sleep(SLEEP_FOR_NEXT_TRY)
        tempdbtools.create_table(TABLES)        
        for fstpaire in allkindofpaire:

            fstpaire1name = re.search(r"[^/]*$",fstpaire[0]).group(0).replace('.','_')
            fstpaire2name = re.search(r"[^/]*$", fstpaire[1]).group(0).replace('.','_')  # for linux
            tableindextoarrayindex.append((allspeices.index(fstpaire1name),allspeices.index(fstpaire2name)))
            
            outfile = open(outputpath+fstpaire1name + fstpaire2name + ".fst"+str(windowWidth)+"_"+str(slideSize), 'w')
            print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile)
            if re.search(r"indvd[^/]*$",fstpaire[0])!=None:
                MethodToSeqpop1="indvd"
            elif re.search(r"pool[^/]+",fstpaire[0])!=None:
                MethodToSeqpop1="pool"
            if re.search(r"indvd[^/]+",fstpaire[1])!=None:
                MethodToSeqpop2="indvd"
            elif re.search(r"pool[^/]+",fstpaire[1])!=None:
                MethodToSeqpop2="pool"
            fst_caculator = Caculators.Caculate_Fst(MethodToSeqpop1=MethodToSeqpop1, MethodToSeqpop2=MethodToSeqpop2)

            fst = Fst() 
            tempdbtools.operateDB("callproc", "mysql_sp_add_column", data=("ninglabvariantdata_tmp", treearrayprename+"treearray", (fstpaire1name[0:5]+fstpaire2name[0:5]), "text", "default null"))
                
            print("startcaculatefst:\n", fstpaire1name,fstpaire[0],'\n', fstpaire2name,fstpaire[1])
            fst.caculateFstAccordingdb(dbtools, chromtable, fstpaire[0], fstpaire[1], fst_caculator, windowWidth,slideSize,minlength)

            winCrossGenome = []
            for chrom in fst.FstMapByChrom.keys():
                for i in range(len(fst.FstMapByChrom[chrom])):
                    if fst.FstMapByChrom[chrom][i][3] != "NA":
                        winCrossGenome.append(fst.FstMapByChrom[chrom][i][3])
            exception = numpy.mean(winCrossGenome)
            std0 = numpy.std(winCrossGenome, ddof=0)
            std1 = numpy.std(winCrossGenome, ddof=1)
            del winCrossGenome
            
            totalChroms = dbtools.operateDB("select","select count(*) from "+chromtable+" where chrlength>="+minlength)[0][0]
            for i in range(0,totalChroms,20):
                currentsql=sql+" where chrlength>="+minlength+" order by "+primaryID+" limit "+str(i)+",20"
                result=dbtools.operateDB("select",currentsql)
                for row in result:
                    currentchrID=row[0]
                    currentchrLen=int(row[1])
                    if currentchrID in fst.FstMapByChrom:
                        for i in range(len(fst.FstMapByChrom[currentchrID])):
                            if fst.FstMapByChrom[currentchrID][i][3] != "NA":
                                zFst = (fst.FstMapByChrom[currentchrID][i][3] - exception) / std1
                                print(currentchrID + "\t" + str(i) + "\t" + str(fst.FstMapByChrom[currentchrID][i][0]) + "\t" + str(fst.FstMapByChrom[currentchrID][i][1])+ "\t" + str(fst.FstMapByChrom[currentchrID][i][2]) + "\t" + '%.15f'%(fst.FstMapByChrom[currentchrID][i][3]) + "\t" + '%.12f'%(zFst), file=outfile) 
                            else:
                                zFst = "NA"
                                print(currentchrID + "\t" + str(i) + "\t" + str(fst.FstMapByChrom[currentchrID][i][0]) + "\t" + str(fst.FstMapByChrom[currentchrID][i][1])+ "\t" + str(fst.FstMapByChrom[currentchrID][i][2]) + "\t" + str(fst.FstMapByChrom[currentchrID][i][3]) + "\t" + str(zFst), file=outfile)                    
#            for chrom in sorted(fst.FstMapByChrom.keys()):        
            sum = 0
            Number = 0
            for chrom in sorted(fst.FstMapByChrom.keys()):
                for i in range(len(fst.FstMapByChrom[chrom])):
                    if fst.FstMapByChrom[chrom][i][3] != 'NA':
                        Number += 1
                        sum += fst.FstMapByChrom[chrom][i][3]
                    tempdbtools.operateDB("insert","insert into "+treearrayprename+"treearray(chrID,winNo,"+fstpaire1name[0:5]+fstpaire2name[0:5]+") values(%s,%s,%s) on duplicate key update "+fstpaire1name[0:5]+fstpaire2name[0:5]+" = '"+str(fst.FstMapByChrom[chrom][i][3])+"'",data=(chrom,str(i),str(fst.FstMapByChrom[chrom][i][3])))
            alldistMap[fstpaire1name+fstpaire2name] = (sum / Number,allspeices.index(fstpaire1name))
            outfile.close()
        for n in alldistMap.keys():
            print(n + "\t" + str(alldistMap[n]), file=open(outputpath+"testdist.txt", 'a'))
        tatalwins = tempdbtools.operateDB("select", "select count(*) from "+treearrayprename+"treearray")[0][0]
        for j in range(0, tatalwins, 100):
            wins = tempdbtools.operateDB("select","select * from "+treearrayprename+"treearray order by chrID asc,winNo asc limit "+str(j) +",100")
            for win in wins:
                abandonthisWin=False
                tmparray=[[0 for x in range(len(allspeices))] for y in range(len(allspeices))]
                
#                print("\t"+arraytitle,file=phyliparrayinfile)
                for i in range(len(win[2:])):
                    tmparray[tableindextoarrayindex[i][0]][tableindextoarrayindex[i][1]]=str(win[i+2])
                    tmparray[tableindextoarrayindex[i][1]][tableindextoarrayindex[i][0]]=str(win[i+2])
                    if win[i+2]==None or win[i+2]=="NA" or win[i+2]=='NULL':
                        abandonthisWin=True
                if abandonthisWin:
                    continue
                print("    "+str(len(allspeices)),file=phyliparrayinfile)
                for i in range(len(allspeices)):
                    tmparray[i][i]='0'
                    try:
                        print(allspeices[i][0:8]+"  "+"\t".join(tmparray[i]),file=phyliparrayinfile)
                    except TypeError:
                        print(i,allspeices,tmparray)
                        print('Error:when making phyliparrayinfile')
                        exit(-1)
        tempdbtools.disconnect()
        phyliparrayinfile.close()
    elif fsttype == 'G' or fsttype == 'g':
        globalFstMapByChrom={}        
        specisnum=str(len(vcffileslist[:]))
        for majorpop in vcffileslist[:]:
            MethodToSeqpop1=None
            fstlist=[]
            vcfname=re.search(r"[^/]*$",majorpop).group(0)
            for othrpop in vcffileslist[:]:
                MethodToSeqpop2=None
                if majorpop == othrpop:
                    continue

                print("startcaculatefst", majorpop, othrpop)
                fstlist.append(Fst())
                if re.search(r"indvd[^/]+",majorpop)!=None:
                    MethodToSeqpop1="indvd"
                elif re.search(r"pool[^/]+",majorpop)!=None:
                    MethodToSeqpop1="pool"
                if re.search(r"indvd[^/]+",othrpop)!=None:
                    MethodToSeqpop2="indvd"
                elif re.search(r"pool[^/]+",othrpop)!=None:
                    MethodToSeqpop2="pool"
                fst_caculator = Caculators.Caculate_Fst(MethodToSeqpop1=MethodToSeqpop1, MethodToSeqpop2=MethodToSeqpop2)
                fstlist[-1].caculateFstAccordingdb(dbtools, chromtable, majorpop, othrpop, fst_caculator, windowWidth,slideSize,minlength)          
                vcfname+=("VS"+re.search(r"[^/]*$",othrpop).group(0)[0])
            outfile=open(outputpath+vcfname+'.gfst'+str(windowWidth)+"_"+str(slideSize)+"_"+specisnum,'w')
            print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile)
            if len(fstlist) != 0:
                for chrom in fstlist[0].FstMapByChrom.keys():
                    globalFstMapByChrom[chrom]=[]
                    for winNo in range(0,len(fstlist[0].FstMapByChrom[chrom])):
                        sumFstInAWin=0
                        Number=0
                        minNumberOfsnp=10000000000000000000000000000000000000000000
                        for i in range(0,len(fstlist)):
                            try:

                                if fstlist[i].FstMapByChrom[chrom][winNo][3]!= 'NA':
                                    Number+=1
                                    minNumberOfsnp=min(minNumberOfsnp,fstlist[i].FstMapByChrom[chrom][winNo][2])
                                    sumFstInAWin+=fstlist[i].FstMapByChrom[chrom][winNo][3]
                            except IndexError:
                                for j in range(0,len(fstlist)):
                                    print(str(i),str(j),vcffileslist[j],chrom,str(winNo),str(len(fstlist[j].FstMapByChrom[chrom])),fstlist[i].FstMapByChrom[chrom][winNo-1])
                                continue# always in the last position,and the value is caculate any way,so can't mispostion.
                        try:
                            gfst=sumFstInAWin/Number
                        except ZeroDivisionError:
                            gfst="NA"
                        globalFstMapByChrom[chrom].append((fstlist[0].FstMapByChrom[chrom][winNo][0],fstlist[0].FstMapByChrom[chrom][winNo][1],minNumberOfsnp,gfst))

                winCrossGenome = []
                for chrom in globalFstMapByChrom.keys():
                    for i in range(len(globalFstMapByChrom[chrom])):
                        if globalFstMapByChrom[chrom][i][3] != "NA":
                            winCrossGenome.append(globalFstMapByChrom[chrom][i][3])
                exception = numpy.mean(winCrossGenome)
                std0 = numpy.std(winCrossGenome, ddof=0)
                std1 = numpy.std(winCrossGenome, ddof=1)
                del winCrossGenome

                totalChroms = dbtools.operateDB("select","select count(*) from "+chromtable+" where chrlength>="+minlength)[0][0]
                for i in range(0,totalChroms,20):
                    currentsql=sql+" where chrlength>="+minlength+" order by "+primaryID+" limit "+str(i)+",20"
                    result=dbtools.operateDB("select",currentsql)
                    for row in result:
                        currentchrID=row[0]
                        currentchrLen=int(row[1])
                        if currentchrID in globalFstMapByChrom:                                
#                for chrom in sorted(globalFstMapByChrom.keys()):
                            for i in range(len(globalFstMapByChrom[currentchrID])):
                                if globalFstMapByChrom[currentchrID][i][3] != "NA":
                                    zgFst = (globalFstMapByChrom[currentchrID][i][3] - exception) / std1
#                                     print(globalFstMapByChrom[currentchrID][i][3])
                                    print(currentchrID + "\t" + str(i) + "\t" + str(globalFstMapByChrom[currentchrID][i][0]) + "\t" + str(globalFstMapByChrom[currentchrID][i][1]) +"\t" + str(globalFstMapByChrom[currentchrID][i][2]) + "\t" + '%.15f'%(globalFstMapByChrom[currentchrID][i][3]) + "\t" + '%.12f'%(zgFst), file=outfile)
                                else:
                                    zgFst = "NA"
                                    print(currentchrID + "\t" + str(i) + "\t" + str(globalFstMapByChrom[currentchrID][i][0]) + "\t" + str(globalFstMapByChrom[currentchrID][i][1]) +"\t" + str(globalFstMapByChrom[currentchrID][i][2]) + "\t" + globalFstMapByChrom[currentchrID][i][3] + "\t" + zgFst, file=outfile)
            outfile.close()
            if len(vcffileslist)==2:
                break
    dbtools.disconnect()
