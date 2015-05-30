'''
Created on 2014-4-24

@author: liurui
'''
from NGS.BasicUtil import Util, VCFutil,Caculators
import NGS.BasicUtil.DBManager as dbm
import os,re


primaryID = "chrID"
class Dstatistics():
    def __init__(self,database="life_pilot",ip="10.2.48.140",usrname="root",pw="1234567",allpopssnptable="derived_alle_ref"):
        super().__init__()
        self.dbtools=dbm.DBTools(ip,usrname,pw,database)
        self.DMapByChrom={}#{chrom1:[(first_snp_pos,last_snp_pos,fst),(),(),...],chrom2:[],....}
        self.allpopssnptable=allpopssnptable
#     def 
    def caculateDstatisticsAccordingdb(self,chromstable,p1name,p2name,p3name,caculator,winwidth,minlengthOfchrom='0',mindeptojudgefix=20):
        totalChroms = self.dbtools.operateDB("select","select count(*) from "+chromstable+" where chrlength>="+minlengthOfchrom)[0][0]
        for i in range(0,totalChroms,20):
            currentsql="select * from " + chromstable+" where chrlength>="+minlengthOfchrom+" order by "+primaryID+" limit "+str(i)+",20"
            chrs=self.dbtools.operateDB("select",currentsql)
            for chrinfo in chrs:
                currentchrID=chrinfo[0]
                currentchrLen=int(chrinfo[2])
                snpsInAChr=self.dbtools.operateDB("select","select snp_pos,ref_base,alt_base,ancestralallel,archicpop,"+p1name+","+p2name+","+p3name+" from "+self.allpopssnptable+" where chrID='"+currentchrID+"'")
                all4popallposInAChr=[]
                if len(snpsInAChr)==0:
                    self.DMapByChrom[currentchrID]=[(0,0,'NA',0,'NA')]
                    continue
                for snp in snpsInAChr:
                    if snp[4]==None or snp[4].strip() == "no covered" or snp[5].strip()=="no covered" or snp[6].strip()=="no covered" or snp[7].strip()=="no covered":
                        continue
                    snp_pos=snp[0]
                    ref_base=snp[1].strip()
                    alt_base=snp[2].strip()
#                     ancestralallel=snp[3].strip().upper()
                    if snp[4]==None or re.search(r'[\w\W]+[,][\w\W]+:\d+,\d+',snp[4])!=None:
                        continue
                    archicpop=re.search(r'([ATCGatcg]+):(\d+),(\d+)',snp[4])
                    archic_base=archicpop.group(1).strip().upper()
                    if len(re.split(r',',alt_base))!=1:
                        continue
#                     if snp[3]!=None and archic_base!= snp[3].strip().upper():
#                         continue
                    if archic_base!=alt_base:
                        continue                    
                    if archicpop.group(2).strip()!='0' and archicpop.group(3).strip()!='0':
                        continue
                    if int(archicpop.group(2).strip())+int(archicpop.group(3).strip())<=mindeptojudgefix:
                        continue
                    elif archicpop.group(2).strip()=='0':
                        if snp[3]!=None or archic_base!= snp[3].strip().upper():
                            continue
                        A_base_idx=1#alt_allele is the ancestral allele
                    elif archicpop.group(3).strip()=='0':
                        if snp[3]!=None or snp[3].strip().upper()!=ref_base:
                            continue
                        A_base_idx=0#ref_allele is the ancestral allele
                    all4popallposInAChr.append((snp_pos,snp[5],snp[6],snp[7],A_base_idx))
                if winwidth==None:
                    for site in all4popallposInAChr:
                        caculator.process(site)
                    else:
                        ABBAcount,BABAcount,D_fixed,D_snp,noofsnps=caculator.getResult()
                    self.DMapByChrom[currentchrID]=[(ABBAcount,BABAcount,D_fixed,noofsnps,D_snp)]
                else:
                    print("ing....")#call self.caculateDstatistics() to caculate D in every windows
                    pass
                        
    def caculateDstatistics(self,p1,p2,p3,p4,caculator,currentchrID,currentchrLen, winwidth=None):
        win = Util.Window()

class MakeDerivedAlleletable():
    def __init__(self, database="life_pilot", ip="10.2.48.140", usrname="root", pw="1234567"):
        super().__init__()
        self.dbtools = dbm.DBTools(ip, usrname, pw, database)

    def createtable(self, tablename="derived_alle_ref"):

        TABLES = {}
        TABLES[tablename] = (
            "CREATE TABLE " + tablename + " ("
            
            " `chrID` varchar(128) NOT NULL DEFAULT '',"
            " `snp_pos` bigint(20) NOT NULL DEFAULT '0',"
            " `snpID` varchar(128) NOT NULL,"
            " `ref_base` tinytext,"
            " alt_base tinytext,"
            " `ancestralallel` tinytext,"
            " `archicpop` varchar(128),"
            
            " PRIMARY KEY (`chrID`,`snp_pos`) "
            ")ENGINE=InnoDB DEFAULT CHARSET=utf8"
            )
        self.dbtools.drop_table(tablename)
        self.dbtools.create_table(TABLES)

    def filldata(self, vcfFileName, depthfileName, tablename="derived_alle_ref", posUniq=True,continuechrom=None,continuepos=None):
        depthfile = Util.GATK_depthfile(depthfileName, depthfileName + ".index")
        depth_linelist=None
        vcffile = open(vcfFileName, 'r')
        vcfline = vcffile.readline()
        while re.search(r'^##', vcfline) != None:
            vcfline = vcffile.readline()
        
        if re.search(r'^#', vcfline) != None:
            poptitlelist = re.split(r'\s+', vcfline.strip())[9:]
            print(poptitlelist)
        else:
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
            exit(-1)   
        for pop in poptitlelist:
            self.dbtools.operateDB("callproc", "mysql_sp_add_column", data=("life_pilot", tablename, pop, "varchar(128)", "default null"))
        popsdata = []#depth for ref or alt
        if continuechrom!=None and continuepos!=None:
            print("filldata",continuechrom,continuepos)
            vcfpossearcher=VCFutil.VCF_Data(vcfFileName)
            vcffile.seek(vcfpossearcher.VcfIndexMap[continuechrom])
            vcfline= vcffile.readline()
            while vcfline:
                vcflist = re.split(r'\s+', vcfline.strip())
                chrom = vcflist[0].strip()
                pos = int(vcflist[1].strip())
                print(chrom,pos)
                if chrom ==continuechrom and pos== continuepos:
                    break
                vcfline =vcffile.readline()
        else:
            justiceGATKorSamtools = vcffile.readline()
            vcflist = re.split(r'\s+', justiceGATKorSamtools.strip())
            dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", vcflist[7])
            refdep = 0;altalleledep = 0
            if dp4 != None:#vcf from samtools 
                print("function for samtools vcf is still need to be finish")
                exit(-1)  
            else:
                chrom = vcflist[0].strip()
                pos = int(vcflist[1].strip())
                snpID = vcflist[2].strip()
                REF = vcflist[3].strip()
                ALT = vcflist[4].strip()
                
                AD_idx = (re.split(":", vcflist[8])).index("AD")#gatk GT:AD:DP:GQ:PL
                sample_idx_in_vcf = 0
                for sample in vcflist[9:]:
    
                    samplename = poptitlelist[sample_idx_in_vcf]
    
                    sample_idx_in_vcf += 1
                    species_idx = depthfile.title.index("Depth_for_" + samplename)
                    if len(re.split(":", sample)) != len(re.split(":", vcflist[8])) and depth_linelist==None:# ./. when lack of variantion information,then consider the depthfile
                        depth_linelist = depthfile.getdepthByPos(chrom, pos)
    
                        if int(depth_linelist[species_idx]) <= 1:
                            popsdata.append('no covered')
                        else:
                            popsdata.append(depth_linelist[species_idx] + ",0")
                        continue
                    elif len(re.split(":", sample)) != len(re.split(":", vcflist[8])) and depth_linelist!=None:
                        if int(depth_linelist[species_idx]) <= 1:
                            popsdata.append('no covered')
                        else:
                            popsdata.append(depth_linelist[species_idx] + ",0")
                        continue                        


                    popsdata.append(re.split(":", sample)[AD_idx])
                depth_linelist=None
                print("insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poptitlelist[:-1]] + poptitlelist[-1:]) + ") select %s,%s,%s,%s,%s," + "%s,"*(len(poptitlelist) - 1) + "%s from dual where not exists( select * from "+tablename+" where "+tablename+".chrID='"+chrom+"' and "+tablename+".snp_pos="+str(pos)+")", (chrom, pos, snpID, REF, ALT) + tuple(popsdata))
                self.dbtools.operateDB("insert", "insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poptitlelist[:-1]] + poptitlelist[-1:]) + ") select %s,%s,%s,%s,%s," + "%s,"*(len(poptitlelist) - 1) + "%s from dual where not exists( select * from "+tablename+" where "+tablename+".chrID='"+chrom+"' and "+tablename+".snp_pos="+str(pos)+")", data=(chrom, pos, snpID, REF, ALT) + tuple(popsdata))
                            
            
        for vcfline in vcffile:
            
            vcflist = re.split(r'\s+', vcfline.strip())
            print(vcfline)
            if posUniq and pos == int(vcflist[1].strip()):
                continue
            chrom = vcflist[0].strip()
            pos = int(vcflist[1].strip())
            snpID = vcflist[2].strip()
            REF = vcflist[3].strip()
            ALT = vcflist[4].strip()
            
            AD_idx = (re.split(":", vcflist[8])).index("AD")#gatk GT:AD:DP:GQ:PL
            sample_idx_in_vcf = 0
            popsdata = []
            for sample in vcflist[9:]:
                samplename = poptitlelist[sample_idx_in_vcf]
                sample_idx_in_vcf += 1
                species_idx = depthfile.title.index("Depth_for_" + samplename)
                if len(re.split(":", sample)) != len(re.split(":", vcflist[8])) and depth_linelist==None:# ./.
                    depth_linelist = depthfile.getdepthByPos(chrom, pos)
                    if int(depth_linelist[species_idx]) <= 1:
                        popsdata.append('no covered')
                    else:
                        popsdata.append(depth_linelist[species_idx] + ",0")
                    continue
                elif len(re.split(":", sample)) != len(re.split(":", vcflist[8])) and depth_linelist!=None:
                    if int(depth_linelist[species_idx]) <= 1:
                        popsdata.append('no covered')
                    else:
                        popsdata.append(depth_linelist[species_idx] + ",0")
                    continue                    
#                 AD_depth = re.split(",", re.split(":", sample)[AD_idx])

                popsdata.append(re.split(":", sample)[AD_idx])
            depth_linelist=None
            print("insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poptitlelist[:-1]] + poptitlelist[-1:]) + ") select %s,%s,%s,%s,%s," + "%s,"*(len(poptitlelist) - 1) + "%s from dual where not exists( select * from "+tablename+" where "+tablename+".chrID='"+chrom+"' and "+tablename+".snp_pos="+str(pos)+")", (chrom, pos, snpID, REF, ALT) + tuple(popsdata))
            self.dbtools.operateDB("insert","insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poptitlelist[:-1]] + poptitlelist[-1:]) + ") select %s,%s,%s,%s,%s," + "%s,"*(len(poptitlelist) - 1) + "%s from dual where not exists( select * from "+tablename+" where "+tablename+".chrID='"+chrom+"' and "+tablename+".snp_pos="+str(pos)+")", data=(chrom, pos, snpID, REF, ALT) + tuple(popsdata))
        depthfile.closedepthfile()
        vcffile.close()
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
        shellstatment=pathtoblastn+" -query "+queryfaFile+" -task blastn -db "+pathtorefdb+" -out "+BlastOutFile +" -outfmt 7 -num_threads 6"
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
        archicpop = VCFutil.VCF_Data(archicpopVcfFile)
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

        
