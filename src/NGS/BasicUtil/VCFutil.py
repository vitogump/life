# -*- coding: UTF-8 -*-
'''
Created on 2013-6-30

@author: rui
'''

import os
import re, numpy, sys, pickle, copy

from NGS.BasicUtil import Util
import NGS.BasicUtil.DBManager as dbm


class AncestralAlleletabletools():
    def __init__(self, database="ninglabvariantdata", ip="10.2.48.140", usrname="root", pw="1234567"):
        super().__init__()
        self.dbtools = dbm.DBTools(ip, usrname, pw, database)
        self.dbname=database
        self.tablename=None

    def createtable(self, vcffilename="derived_alle_ref"):

        TABLES = {}
        tablename=re.search(r'[^/]*$',vcffilename).group(0)
        tablename=re.sub(r"[^\w^\d]","_",tablename)
        self.tablename=tablename
        TABLES[tablename] = (
            "CREATE TABLE " + tablename + " ("
            
            " `chrID` char(128) NOT NULL DEFAULT '',"
            " `snp_pos` bigint(20) NOT NULL DEFAULT '0',"
            " `snpID` char(128) NOT NULL,"
            " `ref_base` varchar(1000),"
            " alt_base varchar(1000),"
            
            " PRIMARY KEY (`chrID`,`snp_pos`) "
            ")ENGINE=InnoDB DEFAULT CHARSET=utf8"
            )
        signal=self.dbtools.create_table(TABLES)
        return signal
    def filldata(self,vcffilename,tablename):
        
        vcffile = open(vcffilename, 'r')
        vcfChromIndex = {}
        line = vcffile.readline()
        
        while re.search(r'^##', line) != None:
            line = vcffile.readline()
        
        if re.search(r'^#', line) != None:
            vcfChromIndex["title"] = re.split(r'\s+', line.strip())
        else:
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
            exit(-1) 
        vcffile.close()
        colslist=vcfChromIndex["title"][9:]
        for col in colslist:
            print("col name",col,"adding to mysql databases")
            self.dbtools.operateDB("callproc", "mysql_sp_add_column", data=(self.dbname, tablename, col, "char(128)", "default null"))
        a=os.system("""awk '$0!~/#/{OFS="\t";for(i=10;i<=NF;i++) printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$i"\t"FS;print ""}' """+vcffilename+">/tmp/tempstep1")
        if a!=0:
            print("error",""""awk '$0!~/#/{OFS="\t";for(i=10;i<=NF;i++) printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$i"\t"FS;print ""}' """+vcffilename+">/tmp/tempstep1")
        loaddatasql = "load data local infile '/tmp/tempstep1' into table " + tablename + " fields terminated by '\\t'"
        shellstatment = "mysql -uroot -p1234567 -D" + self.dbname.strip() + ' -e "' + loaddatasql + '"'
        
        a = os.system(shellstatment)
        if a!=0:
            print("error",shellstatment)
        os.system("rm tempstep1")
#     def filldata(self, vcfFileName, depthfileName, tablename="derived_alle_ref", posUniq=True,continuechrom=None,continuepos=None):
#         depthfile = Util.GATK_depthfile(depthfileName, depthfileName + ".index")
#         depth_linelist=None
#         vcffile = open(vcfFileName, 'r')
#         vcfline = vcffile.readline()
#         while re.search(r'^##', vcfline) != None:
#             vcfline = vcffile.readline()
#         
#         if re.search(r'^#', vcfline) != None:
#             poptitlelist = re.split(r'\s+', vcfline.strip())[9:]
#             print(poptitlelist)
#         else:
#             print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
#             exit(-1)   
#         for pop in poptitlelist:
#             self.dbtools.operateDB("callproc", "mysql_sp_add_column", data=("life_pilot", tablename, pop, "varchar(128)", "default null"))
#         popsdata = []#depth for ref or alt
#         if continuechrom!=None and continuepos!=None:
#             print("filldata",continuechrom,continuepos)
#             vcfpossearcher=VCF_Data(vcfFileName)
#             vcffile.seek(vcfpossearcher.VcfIndexMap[continuechrom])
#             vcfline= vcffile.readline()
#             while vcfline:
#                 vcflist = re.split(r'\s+', vcfline.strip())
#                 chrom = vcflist[0].strip()
#                 pos = int(vcflist[1].strip())
#                 print(chrom,pos)
#                 if chrom ==continuechrom and pos== continuepos:
#                     break
#                 vcfline =vcffile.readline()
#         else:
#             justiceGATKorSamtools = vcffile.readline()
#             vcflist = re.split(r'\s+', justiceGATKorSamtools.strip())
#             dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", vcflist[7])
#             refdep = 0;altalleledep = 0
#             if dp4 != None:#vcf from samtools 
#                 print("function for samtools vcf is still need to be finish")
#                 exit(-1)  
#             else:
#                 chrom = vcflist[0].strip()
#                 pos = int(vcflist[1].strip())
#                 snpID = vcflist[2].strip()
#                 REF = vcflist[3].strip()
#                 ALT = vcflist[4].strip()
#                 
#                 AD_idx = (re.split(":", vcflist[8])).index("AD")#gatk GT:AD:DP:GQ:PL
#                 sample_idx_in_vcf = 0
#                 for sample in vcflist[9:]:
#     
#                     samplename = poptitlelist[sample_idx_in_vcf]
#     
#                     sample_idx_in_vcf += 1
#                     species_idx = depthfile.title.index("Depth_for_" + samplename)
#                     if len(re.split(":", sample)) != len(re.split(":", vcflist[8])) and depth_linelist==None:# ./. when lack of variantion information,then consider the depthfile
#                         depth_linelist = depthfile.getdepthByPos(chrom, pos)
#     
#                         if int(depth_linelist[species_idx]) <= 1:
#                             popsdata.append('no covered')
#                         else:
#                             popsdata.append(depth_linelist[species_idx] + ",0")
#                         continue
#                     elif len(re.split(":", sample)) != len(re.split(":", vcflist[8])) and depth_linelist!=None:
#                         if int(depth_linelist[species_idx]) <= 1:
#                             popsdata.append('no covered')
#                         else:
#                             popsdata.append(depth_linelist[species_idx] + ",0")
#                         continue                        
# 
# 
#                     popsdata.append(re.split(":", sample)[AD_idx])
#                 depth_linelist=None
#                 print("insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poptitlelist[:-1]] + poptitlelist[-1:]) + ") select %s,%s,%s,%s,%s," + "%s,"*(len(poptitlelist) - 1) + "%s from dual where not exists( select * from "+tablename+" where "+tablename+".chrID='"+chrom+"' and "+tablename+".snp_pos="+str(pos)+")", (chrom, pos, snpID, REF, ALT) + tuple(popsdata))
#                 self.dbtools.operateDB("insert", "insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poptitlelist[:-1]] + poptitlelist[-1:]) + ") select %s,%s,%s,%s,%s," + "%s,"*(len(poptitlelist) - 1) + "%s from dual where not exists( select * from "+tablename+" where "+tablename+".chrID='"+chrom+"' and "+tablename+".snp_pos="+str(pos)+")", data=(chrom, pos, snpID, REF, ALT) + tuple(popsdata))
#                             
#             
#         for vcfline in vcffile:
#             
#             vcflist = re.split(r'\s+', vcfline.strip())
#             print(vcfline)
#             if posUniq and pos == int(vcflist[1].strip()):
#                 continue
#             chrom = vcflist[0].strip()
#             pos = int(vcflist[1].strip())
#             snpID = vcflist[2].strip()
#             REF = vcflist[3].strip()
#             ALT = vcflist[4].strip()
#             
#             AD_idx = (re.split(":", vcflist[8])).index("AD")#gatk GT:AD:DP:GQ:PL
#             sample_idx_in_vcf = 0
#             popsdata = []
#             for sample in vcflist[9:]:
#                 samplename = poptitlelist[sample_idx_in_vcf]
#                 sample_idx_in_vcf += 1
#                 species_idx = depthfile.title.index("Depth_for_" + samplename)
#                 if len(re.split(":", sample)) != len(re.split(":", vcflist[8])) and depth_linelist==None:# ./.
#                     depth_linelist = depthfile.getdepthByPos(chrom, pos)
#                     if int(depth_linelist[species_idx]) <= 1:
#                         popsdata.append('no covered')
#                     else:
#                         popsdata.append(depth_linelist[species_idx] + ",0")
#                     continue
#                 elif len(re.split(":", sample)) != len(re.split(":", vcflist[8])) and depth_linelist!=None:
#                     if int(depth_linelist[species_idx]) <= 1:
#                         popsdata.append('no covered')
#                     else:
#                         popsdata.append(depth_linelist[species_idx] + ",0")
#                     continue                    
# #                 AD_depth = re.split(",", re.split(":", sample)[AD_idx])
# 
#                 popsdata.append(re.split(":", sample)[AD_idx])
#             depth_linelist=None
#             print("insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poptitlelist[:-1]] + poptitlelist[-1:]) + ") select %s,%s,%s,%s,%s," + "%s,"*(len(poptitlelist) - 1) + "%s from dual where not exists( select * from "+tablename+" where "+tablename+".chrID='"+chrom+"' and "+tablename+".snp_pos="+str(pos)+")", (chrom, pos, snpID, REF, ALT) + tuple(popsdata))
#             self.dbtools.operateDB("insert","insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poptitlelist[:-1]] + poptitlelist[-1:]) + ") select %s,%s,%s,%s,%s," + "%s,"*(len(poptitlelist) - 1) + "%s from dual where not exists( select * from "+tablename+" where "+tablename+".chrID='"+chrom+"' and "+tablename+".snp_pos="+str(pos)+")", data=(chrom, pos, snpID, REF, ALT) + tuple(popsdata))
#         depthfile.closedepthfile()
#         vcffile.close()
    def getflankseqs(self, chrom,chromlen, snpstartpos, snpendpos, idxedreffilehandler, refindex, flanklen,outfile, tablename="derived_alle_ref"):

        testfile=open("testsnpfile.txt",'a')
        snps = self.dbtools.operateDB("select", "select * from " + tablename + " where chrID='" + chrom + "' and snp_pos>= " + str(snpstartpos) + " and snp_pos<=" + str(snpendpos))
        RefSeqMap = Util.getRefSeqBypos(idxedreffilehandler, refindex, chrom, snpstartpos-flanklen, snpendpos+flanklen,chromlen)
        
        for snp in snps:
            currentsnpPos = snp[1]
            if len(snp[3]) != 1 or len(snp[4]) != 1:
        #                        print(snp[4])
                continue# skip indel
            currentsnpID=chrom+"_"+str(snp[1])
            if currentsnpPos + 25 <= RefSeqMap[chrom][0] + len(RefSeqMap[chrom]) - 1 and currentsnpPos - 25 > RefSeqMap[chrom][0] :
                snpflankseq = ''.join(RefSeqMap[chrom][(currentsnpPos - 25 - RefSeqMap[chrom][0]):(currentsnpPos + 25 - RefSeqMap[chrom][0] + 1)])
                print(currentsnpID,snpflankseq[25],file=testfile)
                snpflankseq=snpflankseq[0:25]+'N'+snpflankseq[26:]
                
            elif currentsnpPos <= RefSeqMap[chrom][0] + len(RefSeqMap[chrom]) - 1 and currentsnpPos + 25 > RefSeqMap[chrom][0] + len(RefSeqMap[chrom]) - 1:
                snpflankseq = ''.join(RefSeqMap[chrom][(currentsnpPos - 25 - RefSeqMap[chrom][0]):(currentsnpPos - RefSeqMap[chrom][0] + 1)])
                print(currentsnpID,snpflankseq[25],file=testfile)
                snpflankseq=snpflankseq[0:25]+'N'
                
            elif currentsnpPos - 25 <= RefSeqMap[chrom][0]:
                snpflankseq = ''.join(RefSeqMap[chrom][(currentsnpPos - RefSeqMap[chrom][0]):(currentsnpPos + 25 - RefSeqMap[chrom][0] + 1)])
                print(currentsnpID,snpflankseq[0],file=testfile)
                snpflankseq = 'N'+snpflankseq[1:26]
                
            else:
                print("what's wrong with the func getflankseqs ?")
                exit(-1)
#            if currentsnpPos + 25 <= RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1 and currentsnpPos - 25 > RefSeqMap[lastchromNo][0] :
#            snpflankseq = ''.join(RefSeqMap[chrom][(currentsnpPos - 25 - RefSeqMap[chrom][0]):(currentsnpPos + 25 - RefSeqMap[chrom][0] + 1)])
#            print(currentsnpID, snpflankseq[25], file=testfile)
#             snpflankseq = snpflankseq[0:25] + 'N' + snpflankseq[26:]
            print(">" + currentsnpID + "\n" + snpflankseq, end='\n', file=outfile)
        testfile.close()
        #                    print("update "+finaltable+" set fafilepos="+str(filepos)+" where snpID='"+currentsnpID+"'")
    def callblast(self,pathtoblastn,pathtorefdb,queryfaFile,BlastOutFile):
        #outfmt chose 6 suggest by zhaoyiqiang
        shellstatment=pathtoblastn+" -query "+queryfaFile+" -task blastn -db "+pathtorefdb+" -out "+BlastOutFile +" -outfmt 7 -num_threads 4"
        print(shellstatment)
        a = os.system(shellstatment)
        if a != 0:
            print("DerivedalleleProcessor : callblast func os.system return not 0")
            exit(-1)
        print(shellstatment,a,"OK")
    def extarctAncestryAlleleFromBlastOut(self,BlastOutFile,ancestryrefFile,ancestryrefidx,tablename="derived_alle_ref",ancestralsnptable=None):
        ancestryreffile=open(ancestryrefFile,'r')
        ancestrysnpflank=open(tablename+"ancestrysnpflank.fa",'w')
        a = os.popen("awk '$1!~/^#/ && $5==1 && $4>26 && $6==0 {print $0}' " + BlastOutFile)
    #    hits=a.readlines()
    
        lastbasesAccur = {}
        onegroup=[]
        revcom=False
    #    initial
        hit = a.readline()
        hitlist = re.split(r"\s+", hit)
    
        sendpos = int(hitlist[9])
        sstartpos = int(hitlist[8])
        qstartpos = int(hitlist[6])
        blastlen=int(hitlist[3])
        snp_loc_s=sstartpos+26-qstartpos
        snpindex = 26 - qstartpos
        if sstartpos > sendpos:
            temp = sstartpos
            sstartpos = sendpos
            sendpos = temp
            revcom=True
        lastsnpID = hitlist[0]
        chrom= hitlist[1]
        RefSeqMap = Util.getRefSeqBypos(refFastahander=ancestryreffile, refindex=ancestryrefidx, currentChromNO=chrom, startpos=sstartpos, endpos=sendpos)
        if revcom:
            tempStr=RefSeqMap[chrom][1:]
            tempStr.reverse()
            RefSeqMap[chrom][1:]=Util.complementary(tempStr)
            revcom=False
            
        lastbasesAccur[RefSeqMap[chrom][snpindex + 1]] = [(chrom, sstartpos, sendpos)]
        onegroup.append((RefSeqMap[chrom][snpindex + 1],blastlen))
        for hit in a:
            print(hit)
            hitlist = re.split(r"\s+", hit)
            chrom = hitlist[1]
            sstartpos = int(hitlist[8])
            sendpos = int(hitlist[9])
            qstartpos = int(hitlist[6])
            blastlen=int(hitlist[3])
            snp_loc_s=sstartpos+26-qstartpos
            snpindex = 26 - qstartpos
            if sstartpos > sendpos:
                temp = sstartpos
                sstartpos = sendpos
                sendpos = temp
                revcom=True
            if lastsnpID == hitlist[0]:
                RefSeqMap = Util.getRefSeqBypos(refFastahander=ancestryreffile, refindex=ancestryrefidx, currentChromNO=chrom, startpos=sstartpos, endpos=sendpos)
                if revcom:
                    tempStr=RefSeqMap[chrom][1:]
                    tempStr.reverse()
                    RefSeqMap[chrom][1:]=Util.complementary(tempStr)
                    revcom=False            
                print(lastsnpID,RefSeqMap[chrom][snpindex + 1],str(snp_loc_s),"".join(RefSeqMap[chrom][1:]),file=ancestrysnpflank)
                if RefSeqMap[chrom][snpindex + 1] in lastbasesAccur:
                    lastbasesAccur[RefSeqMap[chrom][snpindex + 1]].append((chrom, sstartpos, sendpos))
                else:
                    lastbasesAccur[RefSeqMap[chrom][snpindex + 1]] = [(chrom, sstartpos, sendpos)]
                onegroup.append((RefSeqMap[chrom][snpindex + 1],blastlen))
            else:
#                出入数据库 按照不同的主键 即原来是snpid 现在换成别的

                snppos=re.search(r"_(\d+)",lastsnpID).group(1)
                snpChrom=re.search(r"(.+)_(\d+)",lastsnpID).group(1)
                onegroup.sort(key=lambda listRec:listRec[1])                           
                if len(onegroup)==1 or onegroup[0][1]-onegroup[1][1]>=15:#first , only one query id,second longest hit 15 bases greater than the second longest hit
                    if ancestralsnptable!=None and self.dbtools.operateDB("select","select count(*) from "+ancestralsnptable+" where chrID= '"+chrom+"' and snp_start_pos= "+str(snp_loc_s))[0][0]==0:
                        print("update " + tablename + " set ancestralallel='" + onegroup[0][0] + "' where chrID='" + snpChrom + "'and snp_pos="+snppos)
                        self.dbtools.operateDB("update", "update " + tablename + " set ancestralallel='" + onegroup[0][0] + "' where chrID='" + snpChrom + "'and snp_pos="+snppos)
                    else:
                        print("select count(*) from "+ancestralsnptable+" where chrID= '"+chrom+"' and snp_start_pos= "+str(snppos),self.dbtools.operateDB("select","select count(*) from "+ancestralsnptable+" where chrID= '"+chrom+"' and snp_start_pos= "+str(snppos)))
                elif (len(lastbasesAccur.keys()) == 1  and self.dbtools.operateDB("select","select count(*) from "+ancestralsnptable+" where chrID= '"+chrom+"' and snp_start_pos= "+str(snp_loc_s))[0][0]==0):
                    for bases in lastbasesAccur:#only once
                        print("update " + tablename + " set ancestralallel='" + bases + "' where chrID='" + snpChrom + "' and snp_pos="+snppos)
                        self.dbtools.operateDB("update", "update " + tablename + " set ancestralallel='" + bases + "' where chrID='" + snpChrom + "' and snp_pos="+snppos)
                elif len(lastbasesAccur.keys()) == 0:
                    print(" len(lastbasesAccur.keys()) == 0")
                    exit(-1)
                RefSeqMap = Util.getRefSeqBypos(refFastahander=ancestryreffile, refindex=ancestryrefidx, currentChromNO=chrom, startpos=sstartpos, endpos=sendpos)
                if revcom:
                    tempStr=RefSeqMap[chrom][1:]
                    tempStr.reverse()
                    RefSeqMap[chrom][1:]=Util.complementary(tempStr)
                    revcom=False            
                print(hitlist[0],RefSeqMap[chrom][snpindex + 1],str(snp_loc_s),"".join(RefSeqMap[chrom][1:]),file=ancestrysnpflank)
    #            dbtools.operateDB("update", "update " + finaltable + " set chicken='" + RefSeqMap[chrom][snpindex + 1] + "' where snpID='" + hitlist[0] + "'")
                lastsnpID = hitlist[0]
                
                lastbasesAccur.clear()
                lastbasesAccur[RefSeqMap[chrom][snpindex + 1]] = [(chrom, sstartpos, sendpos)]
        print("finish")
        ancestryreffile.close()
    def fillarchicpop(self,archicpopVcfFile,depthFile,chromtable,archicpopNameindepthFile,tablename="derived_alle_ref",archicpopfieldNameintable="archicpop"):
        """
        abandon the snps which exist in archicpopVcfFile but absence in all others pop snp sets 
        """
        depthfile = Util.GATK_depthfile(depthFile, depthFile + ".index")
        species_idx = depthfile.title.index("Depth_for_" + archicpopNameindepthFile)
        archicpop = VCF_Data(archicpopVcfFile)
        totalChroms = self.dbtools.operateDB("select","select count(*) from "+chromtable)[0][0]
        for i in range(0,totalChroms,20):
            currentsql="select * from " + chromtable+" order by chrlength desc limit "+str(i)+",20"
            result=self.dbtools.operateDB("select",currentsql)
            for row in result:
                
                currentchrID=row[0]
                print(currentchrID+":",end="")
                currentchrLen=int(row[2])
                archicpopSeqOfAChr={}
                archicpopSeqOfAChr[currentchrID]=archicpop.getVcfListByChrom(archicpopVcfFile, currentchrID)
                allsnpsInAchr=self.dbtools.operateDB("select","select snp_pos,alt_base from "+tablename+" where chrID='"+currentchrID+"'")
                for snp in allsnpsInAchr:
                    snp_pos=int(snp[0])
                    ALT=snp[1]
                    low=0
                    high=len(archicpopSeqOfAChr[currentchrID])-1
                    while low <=high:
                        mid=(low+high)>>1
                        if archicpopSeqOfAChr[currentchrID][mid][0]< snp_pos:
                            low=mid+1
                        elif archicpopSeqOfAChr[currentchrID][mid][0]> snp_pos:
                            high=mid-1
                        else:#find the pos
                            pos, REF, ALT, INFO,FORMAT,samples = archicpopSeqOfAChr[currentchrID][mid]
                            dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", INFO)
                            refdep=0;altalleledep=0
                            if dp4!=None:#vcf from samtools 
                                refdep = int(dp4.group(1)) + int(dp4.group(2))
                                altalleledep = int(dp4.group(3)) + int(dp4.group(4))    
                            else:
                                AD_idx=(re.split(":",FORMAT)).index("AD")#gatk GT:AD:DP:GQ:PL
                                for sample in samples:
                                    if len(re.split(":",sample))==1:# ./.
                                        continue
                                    AD_depth=re.split(",",re.split(":",sample)[AD_idx])
                                    try :
                                        refdep+=int(AD_depth[0])
                                        altalleledep+=int(AD_depth[1])
                                    except ValueError:
                                        print(sample,end="")
                            popsdata=ALT+":"+str(refdep)+","+str(altalleledep)
                            break
                    else:
                        depth_linelist = depthfile.getdepthByPos(currentchrID, snp_pos)
                        if int(depth_linelist[species_idx]) <= 1:
                            popsdata="no covered"
                        else:
                            popsdata=ALT+":"+depth_linelist[species_idx] + ",0"
#                     print(snp[0],end="\t")     
                    self.dbtools.operateDB("update", "update " + tablename + " set "+archicpopfieldNameintable+" = '" + popsdata+"' where chrID="+"'"+currentchrID+"' and snp_pos="+str(snp[0]))

     
class VCF_Data():
    def __init__(self, vcffileName):
        super().__init__()
        self.VcfMap_AllChrom = {}
        self.VcfIndexMap = {}
        try:
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
        except:
            VCF_Data.indexVCF(VCFName=vcffileName, indexFileName=(vcffileName + ".myindex"))
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
    @staticmethod
    def indexVCF(VCFName, indexFileName):
        """
        {chrom:position_in_file_of_first_SNP_of_this_chrom,chrom:position,,,,,,}
        """
        vcffile = open(VCFName, 'r')
        vcfChromIndex = {}
        line = vcffile.readline()
        
        while re.search(r'^##', line) != None:
            line = vcffile.readline()
        
        if re.search(r'^#', line) != None:
            vcfChromIndex["title"] = re.split(r'\s+', line)
        else:
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
            exit(-1)        
        currentChrom = None
        lastPosition = vcffile.tell()
        print(line)
        line = vcffile.readline()
        print(line)
        while line:      
            linelist = re.split(r"\s+", line)
            if currentChrom != linelist[0]:
                currentChrom = linelist[0]
                vcfChromIndex[currentChrom] = lastPosition
            lastPosition = vcffile.tell()
    
            line = vcffile.readline()
        pickle.dump(vcfChromIndex, open(indexFileName, 'wb'))
        vcffile.close()
    def getVcfListByChrom(self, vcfFileName, chrom,posUniq=True,considerINDEL=False):
        """
            return a list that contain all vcf record of a chrom
        """
        VcfList_A_Chrom = []
        vcfFile = open(vcfFileName, 'r')
        try:
            print("getVcfListByChrom", self.VcfIndexMap[chrom], chrom)            
            vcfFile.seek(self.VcfIndexMap[chrom])
            line = vcfFile.readline().strip()
        except KeyError:
            print(chrom + "didn't find in " + vcfFileName)
            return []
        while line and (re.split(r'\s+', line))[0] == chrom:
            linelist = re.split(r'\s+', line)
            samples=linelist[9:len(linelist)]
            chrom = linelist[0].strip()
            pos = int(linelist[1].strip())
            REF = linelist[3].strip()
            ALT = linelist[4].strip()
            if considerINDEL and len(REF)>1 and len(ALT)>1:
                continue
            INFO = linelist[7]
            FORMAT=linelist[8]
            line = vcfFile.readline()
            if posUniq and VcfList_A_Chrom and pos==VcfList_A_Chrom[0]:
                print("VCFutil unique the vcf pos",line,VcfList_A_Chrom[-1])
                continue
            VcfList_A_Chrom.append((pos, REF, ALT, INFO,FORMAT,samples))
            
        return copy.deepcopy(VcfList_A_Chrom)
        vcfFile.close()

    def getVcfMap(self, vcfFileName):
        """
        this func is from bio\test\posAroundGene\func.py ,and did some improvement,that is add  INFO = collist[7],and add INFO into
        read the vcffile into a map which keys are chrom,values are a list of tuple
        {chrNo:[(pos,REF,ALT,INFO),(pos,REF,ALT,INFO),,,,,],chrNo:[],,,,,,},the order of the tuples in the list,is according pos,
        you we can search a record by  binary chop search
        no matter self.VcfMap_AllChrom has or has not value,the value will be clean
        """
        vcfMap = {}
        vcfFile = open(vcfFileName, 'r')
        
        line = vcfFile.readline()
        while re.search(r'^##', line) != None:
    #        print(line)
            line = vcfFile.readline()
        if re.search(r'^#', line) != None:
            lineslist = vcfFile.readlines()
        else:
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'\n" + line)
            exit(-1)
    #    print("pass")
        currentLine = 0
        totalRecs = len(lineslist)
        while currentLine != totalRecs:
    #        print(currentLine)
            collist = re.split(r'\s+', lineslist[currentLine])
            samples=collist[9:len(collist)]
            chrom = collist[0].strip()
            pos = int(collist[1].strip())
            REF = collist[3].strip()
            ALT = collist[4].strip()
            INFO = collist[7]
            FORMAT = collist[8]
            if chrom in vcfMap:
                vcfMap[chrom].append((pos, REF, ALT, INFO,FORMAT,samples))
            else:
                vcfMap[chrom] = [(pos, REF, ALT, INFO,FORMAT,samples)]
            currentLine += 1
        vcfFile.close()
        self.VcfMap_AllChrom = vcfMap
#         for line in self.VcfMap["scaffold8"]:
#             print(line,file =open("vcfMapdata.txt",'a'))
