'''
Created on 2014-4-24

@author: liurui
'''
from NGS.BasicUtil import Util
import NGS.BasicUtil.DBManager as dbm
import os
import re
class MakeDerivedAlleletable():
    def __init__(self, database="life_pilot", ip="localhost", usrname="root", pw="1234567"):
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

    def filldata(self, vcfFileName, depthfileName, tablename="derived_alle_ref", posUniq=True):
        depthfile = Util.GATK_depthfile(depthfileName, depthfileName + ".index")
        vcffile = open(vcfFileName, 'r')
        vcfline = vcffile.readline()
        while re.search(r'^##', vcfline) != None:
            vcfline = vcffile.readline()
        
        if re.search(r'^#', vcfline) != None:
            poplist = re.split(r'\s+', vcfline.strip())[9:]
            print(poplist)
        else:
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
            exit(-1)   
        for pop in poplist:
            self.dbtools.operateDB("callproc", "mysql_sp_add_column", data=("life_pilot", tablename, pop, "varchar(128)", "default null"))
        popsdata = []#depth for ref or alt
        
        justiceGATKorSamtools = vcffile.readline()
        vcflist = re.split(r'\s+', justiceGATKorSamtools.strip())
        dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", vcflist[7])
        refdep = 0;altalleledep = 0
        if dp4 != None:#vcf from samtools 
            print("function for samtools vcf is still need to be finish")
            exit(-1)
            refdep = int(dp4.group(1)) + int(dp4.group(2))
            altalleledep = int(dp4.group(3)) + int(dp4.group(4))    
        else:
            chrom = vcflist[0].strip()
            pos = int(vcflist[1].strip())
            snpID = vcflist[2].strip()
            REF = vcflist[3].strip()
            ALT = vcflist[4].strip()
            
            AD_idx = (re.split(":", vcflist[8])).index("AD")#gatk GT:AD:DP:GQ:PL
            sample_idx_in_vcf = 0
            for sample in vcflist[9:]:

                samplename = poplist[sample_idx_in_vcf]

                sample_idx_in_vcf += 1
                species_idx = depthfile.title.index("Depth_for_" + samplename)
                if len(re.split(":", sample)) != len(re.split(":", vcflist[8])):# ./. when lack of variantion information,then consider the depthfile
                    depth_linelist = depthfile.getdepthByPos(chrom, pos)

                    if int(depth_linelist[species_idx]) <= 1:
                        popsdata.append('no covered')
                    else:
                        popsdata.append(depth_linelist[species_idx] + ",0")
                    continue
                AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                refdep += int(AD_depth[0])
                altalleledep += int(AD_depth[1])
                popsdata.append(re.split(":", sample)[AD_idx])
            print("insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poplist[:-1]] + poplist[-1:]) + ") values(%s,%s,%s,%s,%s," + "%s,"*(len(poplist) - 1) + "%s)", (chrom, pos, snpID, REF, ALT) + tuple(popsdata))
            self.dbtools.operateDB("insert", "insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poplist[:-1]] + poplist[-1:]) + ") values(%s,%s,%s,%s,%s," + "%s,"*(len(poplist) - 1) + "%s)", data=(chrom, pos, snpID, REF, ALT) + tuple(popsdata))
                            
            
        for vcfline in vcffile:
            
            vcflist = re.split(r'\s+', vcfline.strip())
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
                samplename = poplist[sample_idx_in_vcf]
                sample_idx_in_vcf += 1
                species_idx = depthfile.title.index("Depth_for_" + samplename)
                if len(re.split(":", sample)) != len(re.split(":", vcflist[8])):# ./.
                    depth_linelist = depthfile.getdepthByPos(chrom, pos)
                    if int(depth_linelist[species_idx]) <= 1:
                        popsdata.append('no covered')
                    else:
                        popsdata.append(depth_linelist[species_idx] + ",0")
                    continue
                AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                refdep += int(AD_depth[0])
                altalleledep += int(AD_depth[1])
                popsdata.append(re.split(":", sample)[AD_idx])
            print("insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poplist[:-1]] + poplist[-1:]) + ") values(%s,%s,%s,%s,%s," + "%s,"*(len(poplist) - 1) + "%s)", (chrom, pos, snpID, REF, ALT) + tuple(popsdata))
            self.dbtools.operateDB("insert", "insert into " + tablename + "(chrID,snp_pos,snpID,ref_base,alt_base," + "".join([e + "," for e in poplist[:-1]] + poplist[-1:]) + ") values(%s,%s,%s,%s,%s," + "%s,"*(len(poplist) - 1) + "%s)", data=(chrom, pos, snpID, REF, ALT) + tuple(popsdata))
        depthfile.closedepthfile()
        vcffile.close()
    def getflankseqs(self, chrom,chromlen, snpstartpos, snpendpos, idxedreffilehandler, refindex, flanklen,outfile, tablename="derived_alle_ref"):

        testfile=open("testsnpfile.txt",'a')
        snps = self.dbtools.operateDB("select", "select * from " + tablename + " where chrID='" + chrom + "' and snp_pos>= " + str(snpstartpos) + " and snp_pos<=" + str(snpendpos))
        RefSeqMap = Util.getRefSeqBypos(idxedreffilehandler, refindex, chrom,chromlen, snpstartpos-flanklen, snpendpos+flanklen)
        
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
            snpflankseq = snpflankseq[0:25] + 'N' + snpflankseq[26:]
            print(">" + currentsnpID + "\n" + snpflankseq, end='\n', file=outfile)
        testfile.close()
        #                    print("update "+finaltable+" set fafilepos="+str(filepos)+" where snpID='"+currentsnpID+"'")
    def callblast(self,pathtoblastn,pathtorefdb,queryfaFile,BlastOutFile):
        shellstatment=pathtoblastn+" -query "+queryfaFile+" -task blastn -db "+pathtorefdb+" -out "+BlastOutFile +" -outfmt 7 -num_threads 4"
        print(shellstatment)
        a = os.system(shellstatment)
        if a != 0:
            print("DerivedalleleProcessor : callblast func os.system return not 0")
            exit(-1)
        print(shellstatment,a)
    def extarctAncestryAlleleFromBlastOut(self,BlastOutFile,ancestryrefFile,ancestryrefidx,tablename="derived_alle_ref",ancestralsnptable=None):
        ancestryreffile=open("ancestryrefFile",'r')
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
        snp_loc_s=sstartpos-26+qstartpos
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
            snp_loc_s=sstartpos-26+qstartpos
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
                print(lastsnpID,RefSeqMap[chrom][snpindex + 1],str(snp_loc_s),"".join(RefSeqMap[chrom][1:]),file=open(tablename+"ancestrysnpflank.fa",'a'))
                if RefSeqMap[chrom][snpindex + 1] in lastbasesAccur:
                    lastbasesAccur[RefSeqMap[chrom][snpindex + 1]].append((chrom, sstartpos, sendpos))
                else:
                    lastbasesAccur[RefSeqMap[chrom][snpindex + 1]] = [(chrom, sstartpos, sendpos)]
                onegroup.append((RefSeqMap[chrom][snpindex + 1],blastlen))
            else:
#                出入数据库 按照不同的主键 即原来是snpid 现在换成别的
                snppos=re.search(chrom+r"_(\d+)",lastsnpID).group(1)
                onegroup.sort(key=lambda listRec:listRec[1])
                if len(onegroup)==1 or onegroup[0][1]-onegroup[0][2]>=15:#first , only one query id,second longest hit 15 bases greater than the second longest hit
                    if ancestralsnptable!=None and self.dbtools.operateDB("select","select count(*) from "+ancestralsnptable+" where chrID= '"+chrom+"' and snp_start_pos= '"+str(snppos)+"'")[0][0]==0:
                        self.dbtools.operateDB("update", "update " + tablename + " set ancestralallel='" + onegroup[0][0] + "' where chrID='" + chrom + "' snp_pos='"+snppos+"'")
                elif len(lastbasesAccur.keys()) == 1:
                    for bases in lastbasesAccur:#only once
                        self.dbtools.operateDB("update", "update " + tablename + " set ancestralallel='" + bases + "' where chrID='" + chrom + "' snp_pos='"+snppos+"'")
                elif len(lastbasesAccur.keys()) == 0:
                    exit(-1)
                RefSeqMap = Util.getRefSeqBypos(refFastahander=ancestryreffile, refindex=ancestryrefidx, currentChromNO=chrom, startpos=sstartpos, endpos=sendpos)
                if revcom:
                    tempStr=RefSeqMap[chrom][1:]
                    tempStr.reverse()
                    RefSeqMap[chrom][1:]=Util.complementary(tempStr)
                    revcom=False            
                print(hitlist[0],RefSeqMap[chrom][snpindex + 1],str(snp_loc_s),"".join(RefSeqMap[chrom][1:]),file=open(tablename+"ancestrysnpflank.fa",'a'))
    #            dbtools.operateDB("update", "update " + finaltable + " set chicken='" + RefSeqMap[chrom][snpindex + 1] + "' where snpID='" + hitlist[0] + "'")
                lastsnpID = hitlist[0]
                
                lastbasesAccur.clear()
                lastbasesAccur[RefSeqMap[chrom][snpindex + 1]] = [(chrom, sstartpos, sendpos)]

                ancestryreffile.close()
        
        
        
        
