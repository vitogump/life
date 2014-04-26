'''
Created on 2014-4-24

@author: liurui
'''
from NGS.BasicUtil import Util
import NGS.BasicUtil.DBManager as dbm
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
    def getflankseqs(self, chrom, snpstartpos, snpendpos, tablename="derived_alle_ref", idxedreffilehandler, refindex, flanklen):
        snps = self.dbtools.operateDB("select", "select * from " + tablename + " where chrID='" + chrom + "' and snp_pos>= " + str(snpstartpos) + " and snp_pos<=" + str(snpendpos))
        RefSeqMap = Util.getRefSeqBypos(idxedreffilehandler, refindex, chrom, snpstartpos, snpendpos+flanklen)
        for snp in snps:
            currentsnpPos = int(snp[1])
            if len(snp[3]) != 1 or len(snp[4]) != 1:
        #                        print(snp[4])
                continue# skip indel
            if currentsnpPos + 25 <= RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1 and currentsnpPos - 25 > RefSeqMap[lastchromNo][0] :
                snpflankseq = ''.join(RefSeqMap[currentsnpChrId][(currentsnpPos - 25 - RefSeqMap[currentsnpChrId][0]):(currentsnpPos + 25 - RefSeqMap[currentsnpChrId][0] + 1)])
                print(currentsnpID, snpflankseq[25], file=testfile)
                snpflankseq = snpflankseq[0:25] + 'N' + snpflankseq[26:]
                
            elif currentsnpPos <= RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1 and currentsnpPos + 25 > RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1:
                snpflankseq = ''.join(RefSeqMap[currentsnpChrId][(currentsnpPos - 25 - RefSeqMap[currentsnpChrId][0]):(currentsnpPos - RefSeqMap[currentsnpChrId][0] + 1)])
                print(currentsnpID, snpflankseq[25], file=testfile)
                snpflankseq = snpflankseq[0:25] + 'N'
                
            elif currentsnpPos - 25 <= RefSeqMap[lastchromNo][0]:
                snpflankseq = ''.join(RefSeqMap[currentsnpChrId][(currentsnpPos - RefSeqMap[currentsnpChrId][0]):(currentsnpPos + 25 - RefSeqMap[currentsnpChrId][0] + 1)])
                print(currentsnpID, snpflankseq[0], file=testfile)
                snpflankseq = 'N' + snpflankseq[1:26]
                
            elif currentsnpPos > RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1:
                RefSeqMap, lastchromNo = Util.getRefSeqMap(duckreffile, currentChromNO=currentchrID, preBaseTotal=RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1)
                snpflankseq = ''.join(RefSeqMap[currentsnpChrId][(currentsnpPos - RefSeqMap[currentsnpChrId][0]):(currentsnpPos + 25 - RefSeqMap[currentsnpChrId][0] + 1)])
                print(currentsnpID, snpflankseq[0], file=testfile)
                snpflankseq = 'N' + snpflankseq[1:26]
                
            else:
                print("what's wrong?")
            print(">" + currentsnpID + "\n" + snpflankseq, end='\n', file=outfile)
        #                    print("update "+finaltable+" set fafilepos="+str(filepos)+" where snpID='"+currentsnpID+"'")
            dbtools.operateDB("update", "update " + finaltable + " set fafilepos=" + str(filepos) + " where snpID='" + currentsnpID + "'")
            filepos = int(outfile.tell())
        
        
        
        
