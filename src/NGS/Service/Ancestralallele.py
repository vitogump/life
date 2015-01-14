'''
Created on 2014-11-30

@author: liurui
'''
import os, numpy, sys, re

from NGS.BasicUtil import Util
import NGS.BasicUtil.DBManager as dbm
from NGS.BasicUtil import *


class AncestralAlleletabletools():
    def __init__(self, database="ninglabvariantdata", ip="10.2.48.140", usrname="root", pw="1234567",dbgenome="genomebasicinfo"):
        super().__init__()
        self.dbvariant = dbm.DBTools(ip, usrname, pw, database)
        self.dbgenome=dbm.DBTools(ip, usrname, pw, dbgenome)
        
        #dbtmp means never use the table in the software,you can delete the table anytime without check dependency
        self.dbtmp=dbm.DBTools(ip, usrname, pw, "ninglabvariantdata_tmp")
        self.dbtmpname="ninglabvariantdata_tmp"#
        
        self.dbvariant_name=database
        self.dbgenomename=dbgenome

    def createtable(self, vcffilename="derived_alle_ref"):

        TABLES = {}
        tablename=re.search(r'[^/]*$',vcffilename).group(0)
        tablename=re.sub(r"[^\w^\d]","_",tablename)
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
        signal=self.dbvariant.create_table(TABLES)
        return tablename
    def filldata(self,vcffilename,tablename):
        """
        createtable for every vcf  file,and filldata
        """
        
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
            self.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(self.dbvariant_name, tablename, col, "char(128)", "default null"))
        a=os.system("""awk '$0!~/#/{OFS="\t";printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t";for(i=10;i<=NF;i++) printf $i"\t"FS;print ""}' """+vcffilename+">"+vcffilename+"tempstep1")
        if a!=0:
            print("error",""""awk '$0!~/#/{OFS="\t";printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t";for(i=10;i<=NF;i++) printf $i"\t"FS;print ""}'"""+vcffilename+">"+vcffilename+"tempstep1")
        loaddatasql = "load data local infile '"+vcffilename+"tempstep1"+"' into table " + tablename + " fields terminated by '\\t'"
        shellstatment = "mysql -uroot -p1234567 -D" + self.dbvariant_name.strip() + ' -e "' + loaddatasql + '"'
        
        a = os.system(shellstatment)
        if a!=0:
            print("error",shellstatment)
        os.system("rm "+vcffilename+"tempstep1")
    def getflankseqs(self, chrom,chromlen, snpstartpos, snpendpos, idxedreffilehandler, refindex, flanklen,outfile, tablename="derived_alle_ref"):

        testfile=open("testsnpfile.txt",'a')
        snps = self.dbvariant.operateDB("select", "select * from " + tablename + " where chrID='" + chrom + "' and snp_pos>= " + str(snpstartpos) + " and snp_pos<=" + str(snpendpos))
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
            print("Ancestralallele.py : callblast func os.system return not 0")
            exit(-1)
        print(shellstatment,a,"OK")
    def extarctAncestryAlleleFromBlastOut(self,BlastOutFile,ancestryrefFile,ancestralgenomename,ancestryrefidx,tablename="derived_alle_ref",ancestralsnptable=None):
        ancestryreffile=open(ancestryrefFile,'r')
        ancestrysnpflank=open(tablename+"ancestrysnpflank.fa",'w')
        print(" query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score")
        a = os.popen("awk '$1!~/^#/ && $5==1 && $4>40 && $6==0 {print $0}' " + BlastOutFile)
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
                    if ancestralsnptable!=None and self.dbvariant.operateDB("select","select count(*) from "+ancestralsnptable+" where chrID= '"+chrom+"' and snp_start_pos= "+str(snp_loc_s))[0][0]==0:
                        print("update " + tablename + " set "+ancestralgenomename+" ='" + onegroup[0][0] + "' where chrID='" + snpChrom + "'and snp_pos="+snppos)
                        self.dbvariant.operateDB("update", "update " + tablename + " set "+ancestralgenomename+" ='" + onegroup[0][0] + "' where chrID='" + snpChrom + "'and snp_pos="+snppos)
                    else:
                        print("select count(*) from "+ancestralsnptable+" where chrID= '"+chrom+"' and snp_start_pos= "+str(snppos),self.dbvariant.operateDB("select","select count(*) from "+ancestralsnptable+" where chrID= '"+chrom+"' and snp_start_pos= "+str(snppos)))
                elif (len(lastbasesAccur.keys()) == 1  and self.dbvariant.operateDB("select","select count(*) from "+ancestralsnptable+" where chrID= '"+chrom+"' and snp_start_pos= "+str(snp_loc_s))[0][0]==0):
                    for bases in lastbasesAccur:#only once
                        print("update " + tablename + " set "+ancestralgenomename+" ='" + bases + "' where chrID='" + snpChrom + "' and snp_pos="+snppos)
                        self.dbvariant.operateDB("update", "update " + tablename + " set "+ancestralgenomename+" ='" + bases + "' where chrID='" + snpChrom + "' and snp_pos="+snppos)
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
    #            dbvariant.operateDB("update", "update " + finaltable + " set chicken='" + RefSeqMap[chrom][snpindex + 1] + "' where snpID='" + hitlist[0] + "'")
                lastsnpID = hitlist[0]
                
                lastbasesAccur.clear()
                lastbasesAccur[RefSeqMap[chrom][snpindex + 1]] = [(chrom, sstartpos, sendpos)]
        print("finish")
        ancestryreffile.close()
    def fillAncestral(self,archicpopVcfFile,depthFile,archicpopNameindepthFile,chromtable,toplevelsnptablename="ducksnp_toplevel"):
        """
        abandon the snps which exist in archicpopVcfFile but absence in all others pop snp sets 
        """
        depthfile = Util.GATK_depthfile(depthFile, depthFile + ".index")
        species_idx = depthfile.title.index("Depth_for_" + archicpopNameindepthFile)
        
        archicpop_colname=re.search(r'[^/]*$',archicpopVcfFile).group(0)
        archicpop_colname=re.sub(r"[^\w^\d]","_",archicpop_colname)
        self.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(self.dbvariant_name, toplevelsnptablename, archicpop_colname, "char(128)", "default null"))       
        archicpop = VCFutil.VCF_Data(archicpopVcfFile)
        totalChroms = self.dbgenome.operateDB("select","select count(*) from "+chromtable)[0][0]
        for i in range(0,totalChroms,20):
            currentsql="select * from " + chromtable+" order by chrlength desc limit "+str(i)+",20"
            result=self.dbgenome.operateDB("select",currentsql)
            for row in result:
                
                currentchrID=row[0]
                print(currentchrID+":",end="")
                currentchrLen=int(row[1])
                archicpopSeqOfAChr={}
                archicpopSeqOfAChr[currentchrID]=archicpop.getVcfListByChrom(archicpopVcfFile, currentchrID)
                allsnpsInAchr=self.dbvariant.operateDB("select","select snp_pos,alt_base from "+toplevelsnptablename+" where chrID='"+currentchrID+"'")
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
                                        print("Ancestralallele.fillAncestral except ValueError",sample,end="")
                            popsdata=ALT+":"+str(refdep)+","+str(altalleledep)
                            break
                    else:
                        depth_linelist = depthfile.getdepthByPos(currentchrID, snp_pos)
                        if int(depth_linelist[species_idx]) <= 1:
                            popsdata="no covered"
                        else:
                            popsdata=ALT+":"+depth_linelist[species_idx] + ",0"
                    #change to insert if exist skip
                    print("update",popsdata)
                    self.dbvariant.operateDB("update", "update " + toplevelsnptablename + " set "+archicpop_colname+" = '" + popsdata+"' where chrID="+"'"+currentchrID+"' and snp_pos="+str(snp[0]))
    def leftjoinSelectedTables(self,chromtable,outtable_file_Name,vcftables=[],toplevelsnptable="ducksnp_toplevel",FORMAT="GT:AD:DP:GQ:PL"):
        outfile=open(outtable_file_Name,'w')
        outtable_Name=re.search(r'[^/]*$', outtable_file_Name).group(0)
        outtable_Name=re.sub(r"[^\w^\d]","_",outtable_Name)
        totalChroms = self.dbgenome.operateDB("select","select count(*) from "+chromtable)[0][0]
        outtable_titlelist=[a[0].strip() for a in self.dbvariant.operateDB("select","select column_name  from information_schema.columns where table_schema='"+self.dbvariant_name+"' and table_name='"+toplevelsnptable+"'")]  
        toplevellen=len(outtable_titlelist)
        outtable_titlelist=outtable_titlelist+vcftables
        print(*outtable_titlelist,sep="\t",file=outfile)
        for i in range(0,totalChroms,20):
            currentsql="select * from " + chromtable+" order by chrlength desc limit "+str(i)+",20"
            result=self.dbgenome.operateDB("select",currentsql)
            for row in result:
                sqlselectstatementpart="select t.*"
                sqlfromstatementpart=" from "+toplevelsnptable+" as t "
                currentchrID=row[0]
                print(currentchrID+":",end="")
                #sql statement produce part1
                for vcftable in vcftables:
                    
                    titlelist=[a[0].strip() for a in self.dbvariant.operateDB("select","select column_name  from information_schema.columns where table_schema='"+self.dbvariant_name+"' and table_name='"+vcftable+"'")]              
                    indvdnameslist=titlelist[5:]
                    print(vcftable,indvdnameslist)
                    sqlselectstatementpart=sqlselectstatementpart+","+vcftable.strip()+".alt_base as "+vcftable+"_alt_base"
                    for indvdname in indvdnameslist:
                        sqlselectstatementpart=sqlselectstatementpart+","+vcftable.strip()+"."+indvdname.strip()
                print(sqlselectstatementpart)
                #sql statement produce part2
                for vcftable in vcftables:
                    sqlfromstatementpart=sqlfromstatementpart+" left join "+vcftable.strip()+" using(chrID,snp_pos)"
                #sql where statement append
                sqlstatement=sqlselectstatementpart+sqlfromstatementpart+" where chrID='"+currentchrID+"'"
                print(sqlstatement)
                allsnpOfJoinTableinAchr=self.dbvariant.operateDB("select",sqlstatement)
                #process value,merge into one col
                NumOfColOftoplevel_fix=len(self.dbvariant.operateDB("select","select column_name  from information_schema.columns where table_schema='"+self.dbvariant_name+"' and table_name='"+toplevelsnptable+"'"))
                for rec in allsnpOfJoinTableinAchr:
                    NumOfColOftoplevel=NumOfColOftoplevel_fix
                    recToPrint=list(rec[0:NumOfColOftoplevel])
                    print("recToPrintpre",recToPrint,"rec",rec)
                    for vcftable in vcftables:
                        titlelist=[a[0].strip() for a in self.dbvariant.operateDB("select","select column_name  from information_schema.columns where table_schema='"+self.dbvariant_name+"' and table_name='"+vcftable+"'")]
                        indvdnameslist=titlelist[5:]
                        refdep=0;altalleledep=0
                        print("lastpart",rec[NumOfColOftoplevel:NumOfColOftoplevel+1+len(indvdnameslist)])
                        if rec[NumOfColOftoplevel]==None:
                            recToPrint=recToPrint+["unknow"]
                            NumOfColOftoplevel=NumOfColOftoplevel+1+len(indvdnameslist)
#                             continue
                        else:
                            ALT=rec[NumOfColOftoplevel]
                            AD_idx=(re.split(":",FORMAT)).index("AD")
                            for sample in rec[NumOfColOftoplevel+1:NumOfColOftoplevel+1+len(indvdnameslist)]:
                                if len(re.split(":",sample))==1:
                                    continue
                                AD_depth=re.split(",",re.split(":",sample)[AD_idx])
                                try:
                                    refdep+=int(AD_depth[0])
                                    altalleledep+=int(AD_depth[1])
                                except ValueError:
                                    print("ValueError",sample,end="")
                            popsdata=ALT+":"+str(refdep)+","+str(altalleledep)
                            recToPrint=recToPrint+[popsdata]
                            NumOfColOftoplevel=NumOfColOftoplevel+1+len(indvdnameslist)
                    print(*recToPrint,sep="\t",file=outfile)
        self.dbtmp.operateDB("copytableschema","create table "+outtable_Name+" like "+self.dbvariant_name.strip()+"."+toplevelsnptable.strip())
        for colname in outtable_titlelist[toplevellen-1:]:
            self.dbtmp.operateDB("callproc", "mysql_sp_add_column",data=(self.dbvariant_name, outtable_Name, colname, "char(128)", "default null"))
        loaddatasql="load data local infile '"+outtable_file_Name+"' into table "+outtable_Name+" fields terminated by '\\t'"
        a=os.system("mysql -uroot -p1234567 -D" + self.dbtmpname.strip() + ' -e "' + loaddatasql + '"')
        if a!=0:
            print(a,"maybe error")
        else:
            print("finish")

                        
                