'''
Created on 2014-11-30

@author: liurui
'''
import copy, time, pysam
import os, numpy, sys, re
import pickle, configparser

from NGS.BasicUtil import Util, VCFutil
import NGS.BasicUtil.DBManager as dbm
import web.dba as dba










currentpath=os.path.realpath(__file__)
currentpath[:currentpath.find("life/src")]+"life/com/config.properties"
cfparser = configparser.ConfigParser()
cfparser.read(currentpath[:currentpath.find("life/src")]+"life/com/config.properties")

password=cfparser.get("mysqldatabase","password")
mindepthforJudgefixref=10
class dynamicInsertUpdateAncestralContext():#here use the fast edition
    def __init__(self,dbvariantstools,refseqfa,toplevelsnptablename="mspsgjlksy10pop_toplevel_pekingduckref_new",insertauthority=True,changetableauthority=True):
        self.refseqfahandler=open(refseqfa,'r')
        self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE={}#{vcfname:[vcfobj,pysam.samfile1,pysam.samfile2,,,,],,,,}
        self.mapofSNPrecsForeachVCFpop_mapBYchrom={}
        try:
            self.refseqfafasteridx=pickle.load(open(refseqfa+".myfasteridx",'rb'))
        except IOError:
            Util.generateFasterRefIndex(refseqfa, refseqfa+".myfasteridx")
            self.refseqfafasteridx=pickle.load(open(refseqfa+".myfasteridx",'rb'))
#         self.vcfnameKEY_vcfobjVALUE={}#{vcfname:vcfobj,,,,}
        self.currentchrLen=None
        self.toplevelsnptablename=toplevelsnptablename
        self.dbvariantstools=dbvariantstools
        self.toplevelsnptable_titlelist=[a[0].strip() for a in self.dbvariantstools.operateDB("select", "select column_name  from information_schema.columns where table_schema='" + Util.vcfdbname + "' and table_name='" + toplevelsnptablename + "'")]
        outgroupBAMconfigs=re.split(r";",Util.outgroupVCFBAMconfig_beijingref)
        for configfile in outgroupBAMconfigs:
            fp=open(configfile,'r')
            for line in fp:
                vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                if vcffilename_obj!=None:
                    vcfname=vcffilename_obj.group(1).strip()
                    archicpop_colname=re.search(r'[^/]*$',vcfname).group(0)
                    archicpop_colname=re.sub(r"[^\w^\d]","_",archicpop_colname)
                    print("init dynamicInsertUpdateAncestralContext",vcfname,archicpop_colname,end="\t")
                    for outgroupname_ALT in self.toplevelsnptable_titlelist[6::2]:
                        if archicpop_colname+"_alt"==outgroupname_ALT:
                            break
                    else:
                        print(archicpop_colname+"_alt does not exist in topleveltable" )
                        if changetableauthority:
                            self.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(Util.vcfdbname, toplevelsnptablename, archicpop_colname+"_alt", "char(128)", "default null"))
                            self.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(Util.vcfdbname, toplevelsnptablename, archicpop_colname+"_dep", "char(128)", "default null"))
                        else:
                            pass
                    #always a
                    self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                    self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcffilename_obj.group(1).strip()))
                elif line.split():
                    print("dynamicInsertUpdateAncestralContext",line.strip())
                    self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
                    
            fp.close()
    def __del__(self):
        for vcfname in self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE.keys():
            for pysamobj in self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][1:]:
                pysamobj.close()
    def getRECsforCHR(self,chrom,currentchrLen=None):
        self.currentchrLen=currentchrLen
        self.mapofSNPrecsForeachVCFpop_mapBYchrom={}
        self.mapofSNPrecsForeachVCFpop_mapBYchrom[chrom]={}
        for vcfname_tmp in self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE.keys():
            self.mapofSNPrecsForeachVCFpop_mapBYchrom[chrom][vcfname_tmp]=self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname_tmp][0].getVcfListByChrom(chrom,MQfilter=None)
    def insertorUpdatetopleveltable(self,PopSnpAligned,flankseqfafile,SNP_flanklen):
        """
            first element of a snp is the snp position
            PopSnpAligned are required to be ordered
            func is fill context or/and ancestral (using outgroupVCFBAMconfig_beijingref in config.properties) 
        """
        for chrom in PopSnpAligned.keys():
            #pre process for context
            RefSeqMap=Util.getRefSeqBypos_faster(self.refseqfahandler, self.refseqfafasteridx, chrom, PopSnpAligned[chrom][0][0]-SNP_flanklen, PopSnpAligned[chrom][-1][0]+SNP_flanklen, self.currentchrLen)
            #pre process for context end
            SNPrec_of_one_chrom_invcf=[]
            insertsql_statement_list=[];updatesql_statement_list=[]
            insertsql_data_list=[];updatesql_date_list=[]
            for snp in PopSnpAligned[chrom]:
                print("insertorUpdatetopleveltable",snp)
                snp_pos=snp[0]
                snp_recINtopleveltable=self.dbvariantstools.operateDB("select","select * from "+self.toplevelsnptablename+" where chrID='"+chrom+"' and snp_pos="+str(snp_pos))
                if not snp_recINtopleveltable or snp_recINtopleveltable==0:#empty
                    insertsql_statement="insert into "+self.toplevelsnptablename + " (chrID,snp_pos,snpID,ref_base,alt_base,context,"
                    insertsql_date=[chrom,snp_pos,".",snp[1],snp[2]]
                    updatesql_statement=None
                else:
                    """
                    snp_recINtopleveltable [('KB743962.1', 60011, '.', 'A', 'T', 'TAA', None, None, 'T', '39,0', None, None, 'T', '21,0')]

                    """
                    print("snp_recINtopleveltable",snp_recINtopleveltable)
                    insertsql_statement=None
                    updatesql_date=[]
                    updatevcflist=[]
                    alt_idx=6
                    for outgroupname_ALT in snp_recINtopleveltable[0][6::2]:
                        if None==outgroupname_ALT:
                            updatevcflist.append(self.toplevelsnptable_titlelist[alt_idx][:-4])#know which field need to be filled
                        alt_idx+=2
                        
                    updatesql_statement="update "+ self.toplevelsnptablename+" set "
                for vcfname in self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE.keys():
                    archicpop_colname=re.search(r'[^/]*$',vcfname).group(0)
                    archicpop_colname=re.sub(r"[^\w^\d]","_",archicpop_colname)
                    if insertsql_statement==None and archicpop_colname in updatevcflist:
                        updatesql_statement+=(archicpop_colname+"_alt =%s,"+archicpop_colname+"_dep =%s,")
                    elif updatesql_statement!=None and insertsql_statement==None:
                        print("warning! "+archicpop_colname+" not in the ",updatevcflist)
                        print("updatesql_statement is ",updatesql_statement)
                        continue# this vcf is filled
                    elif updatesql_statement==None and insertsql_statement!=None:
                        insertsql_statement+=(archicpop_colname+"_alt,"+archicpop_colname+"_dep,")
                    
                    if chrom in  self.mapofSNPrecsForeachVCFpop_mapBYchrom:
                        pass
                    else:
                        self.getRECsforCHR(chrom)
                    SNPrec_of_one_chrom_invcf=self.mapofSNPrecsForeachVCFpop_mapBYchrom[chrom][vcfname]
                    """state: updatesql_statement= update toplevelsnptablename set processed_vcfname1_alt=%s,processed_vcfname1_dep=%s,
                            insertsql_statement=insert into toplevelsnptablename (chrID,snp_pos,snpID,ref_base,alt_base,context,processed_vcfname_m_alt,processed_vcfname_m_dep,
                    """
                    low=0;ALT=snp[2];high=len(SNPrec_of_one_chrom_invcf)-1
                    while low<=high:
                        mid=(low+high)>>1
                        if SNPrec_of_one_chrom_invcf[mid][0]<snp_pos:
                            low=mid+1
                        elif SNPrec_of_one_chrom_invcf[mid][0]>snp_pos:
                            high=mid-1
                        else:#find the pos
#                             print("here is a  bug still not finished",mid)
                            pos, REF, archicpop_ALT, INFO,FORMAT,samples = SNPrec_of_one_chrom_invcf[mid]
                            dp4=re.search(r"DP4=(\d*)", INFO)
                            AF=re.search(r"AF=([\d\.e-]+)", INFO).group(1)
                            refdep=0;altalleledep=0
                            if dp4!=None and AF!=None: #samtools
                                refdep = int(dp4.group(1))*(1- float(AF))
                                altalleledep = int(dp4.group(1)) * float(AF)
                            else:#vcf from indvd 
                                if "AD" not in FORMAT:
                                    print("AD not in ",FORMAT )
                                    continue
                                AD_idx=(re.split(":",FORMAT)).index("AD")#gatk GT:AD:DP:GQ:PL 
                                for sample in samples:
                                    if len(re.split(r":",sample))==1:# ./.
                                        continue
                                    AD_depth=re.split(r",",re.split(":",sample)[AD_idx])
                                    try:
                                        refdep+=int(AD_depth[0])
                                        altalleledep+=int(AD_depth[1])
                                    except ValueError:
                                        print("Ancestralallele.fillAncestral except ValueError",sample,end="")
                            popsdata_alt=archicpop_ALT
                            popsdata_dep=str(refdep)+","+str(altalleledep)
                            break
                    else:
                        sum_depth=0
                        popsdata_alt=ALT
                        for samfile in self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][1:]:
                            ACGTdep=samfile.count_coverage(chrom,snp_pos-1,snp_pos)
                            for dep in ACGTdep:
                                sum_depth+=dep[0]
                        if sum_depth<4:
                            popsdata_dep="no covered"
                        else:
                            popsdata_dep=str(sum_depth) + ",0"
                    #continnue sql_statement
                    """state: updatesql_date=[],insertsql_date=[chrom,snp_pos,".",snp[1],snp[2]]
                    updatesql_statement= update toplevelsnptablename set processed_vcfname1_alt=%s,processed_vcfname1_dep=%s,
                    insertsql_statement=insert into toplevelsnptablename (chrID,snp_pos,snpID,ref_base,alt_base,context,processed_vcfname_m_alt,processed_vcfname_m_dep,
                    """
                    if updatesql_statement!=None and insertsql_statement==None:
                        updatesql_date+=[popsdata_alt,popsdata_dep]
                        print("updatesql_date",updatesql_date)
                    elif  updatesql_statement==None and insertsql_statement!=None:
                        insertsql_date+=[popsdata_alt,popsdata_dep]
                    else:
                        print("error: updatesql_statement!=None or insertsql_statement!=None")
                #for one snp
                #process context 
                if (updatesql_statement!=None and  snp_recINtopleveltable[0][5]==None) or insertsql_statement!=None:
                    snpID=chrom+"_"+str(snp[0])
                    if snp_pos + SNP_flanklen<=RefSeqMap[chrom][0]+len(RefSeqMap[chrom])-1 and snp_pos- SNP_flanklen>RefSeqMap[chrom][0]:
                        snpflankseq="".join(RefSeqMap[chrom][(snp_pos-SNP_flanklen-RefSeqMap[chrom][0]):(snp_pos+SNP_flanklen-RefSeqMap[chrom][0]+1)])
                        if insertsql_statement!=None:
                            insertsql_date.insert(5, snpflankseq[SNP_flanklen-1:SNP_flanklen+2])                        
                        elif snp_recINtopleveltable[0][5]==None:
                            updatesql_statement+="context=%s"
                            updatesql_date.append(snpflankseq[SNP_flanklen-1:SNP_flanklen+2])
                            print("updatesql_date",updatesql_date)
                        currentsnpID=chrom+"_"+str(snp_pos)+snpflankseq[SNP_flanklen]+":"+snp[1]
                        snpflankseq=snpflankseq[0:SNP_flanklen]+'N'+snpflankseq[SNP_flanklen+1:]

                    elif snp_pos <=RefSeqMap[chrom][0]+len(RefSeqMap[chrom])-1 and snp_pos-SNP_flanklen>RefSeqMap[chrom][0]:
                        snpflankseq="".join(RefSeqMap[chrom][(snp_pos-SNP_flanklen-RefSeqMap[chrom][0]):(snp_pos-RefSeqMap[chrom][0]+1)])
                        if insertsql_statement!=None:
                            insertsql_date.insert(5,snpflankseq[SNP_flanklen-1:SNP_flanklen+1]+"N")                        
                        elif snp_recINtopleveltable[0][5]==None:
                            updatesql_statement+="context=%s"
                            updatesql_date.append(snpflankseq[SNP_flanklen-1:SNP_flanklen+1]+"N")
                            print("updatesql_date",updatesql_date)
                        currentsnpID=chrom+"_"+str(snp_pos)+snpflankseq[SNP_flanklen]+":"+snp[1]
                        snpflankseq=snpflankseq[0:SNP_flanklen]+'N'
                        
                    elif snp_pos-SNP_flanklen<=RefSeqMap[chrom][0] and snp_pos + SNP_flanklen<=RefSeqMap[chrom][0]+len(RefSeqMap[chrom])-1:
                        snpflankseq="".join(RefSeqMap[chrom][(snp_pos - RefSeqMap[chrom][0]):(snp_pos + SNP_flanklen - RefSeqMap[chrom][0] + 1)])
                        if insertsql_statement!=None:
                            insertsql_date.insert(5,"N"+snpflankseq[0:2])                        
                        elif snp_recINtopleveltable[0][5]==None:
                            updatesql_statement+="context=%s"
                            updatesql_date.append("N"+snpflankseq[0:2])
                            print("updatesql_date",updatesql_date)
                        currentsnpID=chrom+"_"+str(snp_pos)+snpflankseq[0]+":"+snp[1]
                        snpflankseq = 'N'+snpflankseq[1:SNP_flanklen+1]
                        
                    else:
                        print(len(RefSeqMap[chrom]))
                        currentsnpID=">error"
                        snpflankseq="error"
                        print("Error! what's wrong with the func insertorUpdatetopleveltable?")
                    print(">" + currentsnpID + "\n" + snpflankseq, end='\n', file=flankseqfafile)
                #context process end ,append into multiple statements    
                if updatesql_statement!=None:
                    if snp_recINtopleveltable[0][5]!=None:
                        updatesql_statement=updatesql_statement[:-1]
                    updatesql_statement+=" where snp_pos="+str(snp_pos)+" and chrID='"+chrom+"'"
                    updatesql_statement_list.append(updatesql_statement)
                    print("updatesql_date_list",updatesql_date_list)
                    updatesql_date_list.append(tuple(updatesql_date))
                    print("updatesql_date_list",updatesql_date_list)
                elif insertsql_statement!=None:
                    insertsql_statement=insertsql_statement[:-1]+")VALUES ("+"%s,"*(len(insertsql_date))
                    insertsql_statement=insertsql_statement[:-1]+")"
                    insertsql_statement_list.append(insertsql_statement)
                    insertsql_data_list.append(tuple(insertsql_date))
#                     self.dbvariantstools.operateDB("",insertsql_statement,data=tuple(insertsql_date))
            if insertsql_statement_list !=[]:
#                 print(insertsql_statement_list)
                self.dbvariantstools.operateDB("insert",*insertsql_statement_list,data=insertsql_data_list)
#                 snp=self.dbvariantstools.operateDB("select","select * from "+self.toplevelsnptablename+" where chrID='"+insertsql_data_list[0][0]+"' and snp_pos='"+str(insertsql_data_list[0][1])+"'")
#                 print(snp)
            elif updatesql_statement_list !=[] and len(updatesql_date_list)==len(updatesql_statement_list):
                print(updatesql_statement_list,updatesql_date_list)
                self.dbvariantstools.operateDB("update",*updatesql_statement_list,data=updatesql_date_list)
                
#             print("insertsql_data_list",insertsql_data_list)
#             print("updatesql_date_list",updatesql_date_list)
#             snp=self.dbvariantstools.operateDB("select","select * from "+self.topleveltablejudgeancestralname+" where chrID='"+insertsql_data_list[0][0]+"' and snp_pos='"+str(insertsql_data_list[0][1])+"'")

class AncestralAlleletabletools():
    def __init__(self, database="ninglabvariantdata", ip="10.2.48.140", usrname="root", pw="1234567",dbgenome="genomebasicinfo"):
        super().__init__()
        self.dbvariant = dbm.DBTools(ip, usrname, pw, database)
        self.dbgenome=dbm.DBTools(ip, usrname, pw, dbgenome)
#         self.forchenli=open("chenliidentifysheep",'w')
        #dbtmp means never use the table in the software,you can delete the table anytime without check dependency
        self.dbtmp=dbm.DBTools(ip, usrname, pw, Util.ghostdbname)
        self.dbtmpname=Util.ghostdbname#
        self.session=dba.getVarSession()
        self.dbvariant_name=database
        self.dbgenomename=dbgenome

    def createtable(self, vcffilename="derived_alle_ref",drop=False):

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

        signal,key=self.dbvariant.create_table(TABLES)
        if signal!="OK":
            print("signal",signal,drop)
            while signal=="already exist" and drop:
                self.dbvariant.drop_table(key)
                signal,key=self.dbvariant.create_table(TABLES)
                print(signal)
            else:
                if signal=="already exist" and not drop:
                    exit(-1)
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
        firstSNPrec=re.split(r"\s+",vcffile.readline().strip())
        fields_OF_INFO=re.split(r";",firstSNPrec[7])
 
        
        vcffile.close()
#         colslist=vcfChromIndex["title"][9:]
        if re.search(r"indvd[^/]+",vcffilename)!=None:
            print("indvd")
            colslist=["DP","AF","AN"]
            colidxmap={}
            for e in fields_OF_INFO:
                for col in colslist:
                    if (col+"=" in e) and re.search(r"^"+col+"=",e)!=None :
                        colidxmap[col]=fields_OF_INFO.index(e)
            for col in colslist:
                print("col name",col,"adding to mysql databases")
                self.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(self.dbvariant_name, tablename, col, "char(128)", "default null"))
            a=os.system(""" awk '$0!~/#/&&length($5)==1{OFS="\t";print $1,$2,$3,$4,$5,$8}' """+vcffilename+""" |awk '{OFS="\t";split($6,myarr,";");print $1,$2,$3,$4,$5,myarr["""+str(colidxmap[colslist[0]])+"""],myarr["""+str(colidxmap[colslist[1]])+"""],myarr["""+str(colidxmap[colslist[2]])+"""]}' |sed 's/DP=//g'|sed 's/AF=//g'|sed 's/AN=//g'|awk '{OFS="\t";if($7==1){$7=1};print $0}' > """+vcffilename+"tempstep1")
            if a!=0:
                print("error",""" awk '$0!~/#/&&length($5)==1{OFS="\t";print $1,$2,$3,$4,$5,$8}' """+vcffilename+""" |awk '{OFS="\t";split($6,myarr,";");print $1,$2,$3,$4,$5,myarr["""+str(colidxmap[colslist[0]])+"""],myarr["""+str(colidxmap[colslist[1]])+"""],myarr["""+str(colidxmap[colslist[2]])+"""]}' |sed 's/DP=//g'|sed 's/AF=//g'|sed 's/AN=//g'|awk '{OFS="\t";if($7==1){$7=1};print $0}' > """+vcffilename+"tempstep1")
        
        elif re.search(r"pool[^/]+",vcffilename)!=None:
            print("pool if NF>10 need recode the programm")
#             colslist=vcfChromIndex["title"][9:]
            colslist=["DP","AF"]
            colidxmap={}
            for e in fields_OF_INFO:
                for col in colslist:
                    if (col+"=" in e) and re.search(r"^"+col+"=",e)!=None :
                        colidxmap[col]=fields_OF_INFO.index(e)
            for col in colslist:
                print("col name",col,"adding to mysql databases")
                self.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(self.dbvariant_name, tablename, col, "char(128)", "default null"))
            a=os.system("""awk '$0!~/#/&&length($5)==1{OFS="\t";print $1,$2,$3,$4,$5,$10}' """+vcffilename+""" |awk '{OFS="\t";split($6,myarr,":");split(myarr[2],mydep,",");print $1,$2,$3,$4,$5,(myarr["""+str(colidxmap[1])+"""]),mydep[2]/(mydep[1]+mydep[2])}'>"""+vcffilename+"tempstep1")
            if a!=0:
                print("""awk '$0!~/#/&&length($5)==1{OFS="\t";print $1,$2,$3,$4,$5,$10}' """+vcffilename+""" |awk '{OFS="\t";split($6,myarr,":");split(myarr[2],mydep,",");print $1,$2,$3,$4,$5,mydep[2]/(mydep[1]+mydep[2]),(myarr[3])}'>"""+vcffilename+"tempstep1")    
        else:
            print("skip",vcffilename)
            return
        loaddatasql = "load data local infile '"+vcffilename+"tempstep1"+"' into table " + tablename + " fields terminated by '\\t'"
        shellstatment = "mysql -uroot -p"+password+" -D" + self.dbvariant_name.strip() + ' -e "' + loaddatasql + '"'
        
        a = os.system(shellstatment)
        if a!=0:
            print("error",shellstatment)
        os.system("rm "+vcffilename+"tempstep1")
    def getflankseqstooutfile(self, chrom,chromlen, startpostocollecteSNP, endpostocollectSNP, idxedreffilehandler,ancestralgenomenameaddtotable, refindex, flanklen,outfile,vfilelinelists,info=None, tablename="derived_alle_ref"):
                            
        RefSeqMap = Util.getRefSeqBypos_faster(idxedreffilehandler, refindex, chrom, startpostocollecteSNP-flanklen, endpostocollectSNP+flanklen,chromlen)
        for snp in vfilelinelists:
            """
            snp=
            [chrom,  pos, snpID,  refbase, altbase,qual,fiter,info,,,]
            """
            
            currentsnpPos = int(snp[1])
            if currentsnpPos>=chromlen:
                print("warnning,this loci may not match this fasta file")
                continue
            if False and len(snp[3]) != 1 :
        #                        print(snp[4])
                """not in use
                """
                continue# skip indel
            currentsnpID=chrom+"_"+str(snp[1])+"_" +snp[7][:70] if info==None else chrom+"_"+str(snp[1])+"_" +re.search(r"",snp[7]).group(1)
            if currentsnpPos + flanklen < RefSeqMap[chrom][0] + len(RefSeqMap[chrom]) - 1 and currentsnpPos - flanklen > RefSeqMap[chrom][0] :
                snpflankseq = ''.join(RefSeqMap[chrom][(currentsnpPos - flanklen - RefSeqMap[chrom][0]):(currentsnpPos + flanklen - RefSeqMap[chrom][0] + 1)])
#                 self.dbvariant.operateDB("update","update "+tablename+" set context='"+snpflankseq[flanklen-1:flanklen+2]+"' where chrID='"+ chrom + "' and snp_pos= "+str(currentsnpPos))
                currentsnpID+="_"+snpflankseq[flanklen]+":"+snp[3]+snp[4]
                snpflankseq=snpflankseq[0:flanklen]+'N'+snpflankseq[flanklen+1:]
                
            elif currentsnpPos <= RefSeqMap[chrom][0] + len(RefSeqMap[chrom]) - 1 and currentsnpPos - flanklen > RefSeqMap[chrom][0]:
                snpflankseq = ''.join(RefSeqMap[chrom][(currentsnpPos - flanklen - RefSeqMap[chrom][0]):(currentsnpPos - RefSeqMap[chrom][0] + 1)])
#                 print(currentsnpID,snpflankseq[flanklen],file=testfile)
#                 self.dbvariant.operateDB("update","update "+tablename+" set context='"+snpflankseq[flanklen-1:flanklen+1]+"N' where chrID='"+ chrom + "' and snp_pos= "+str(currentsnpPos))
                currentsnpID+="_"+snpflankseq[flanklen]+":"+snp[3]+snp[4]
                snpflankseq=snpflankseq[0:flanklen]+'N'
                
            elif currentsnpPos - flanklen <= RefSeqMap[chrom][0] and currentsnpPos + flanklen<=RefSeqMap[chrom][0]+len(RefSeqMap[chrom])-1:
                snpflankseq = ''.join(RefSeqMap[chrom][(currentsnpPos - RefSeqMap[chrom][0]):(currentsnpPos + flanklen - RefSeqMap[chrom][0] + 1)])
#                 print(currentsnpID,snpflankseq[0],file=testfile)
#                 self.dbvariant.operateDB("update","update "+tablename+" set context='N"+snpflankseq[0:2]+"' where chrID='"+ chrom + "' and snp_pos= "+str(currentsnpPos))
                currentsnpID+="_"+snpflankseq[0]+":"+snp[3]+snp[4]
                snpflankseq = 'N'+snpflankseq[1:flanklen+1]
                
            else:
                print(snp,currentsnpPos,RefSeqMap[chrom][0],len(RefSeqMap[chrom])-1)
                print("what's wrong with the func getflankseqs ?")
                exit(-1)

            print(">" + currentsnpID + "\n" + snpflankseq, end='\n', file=outfile)
    def getflankseqs(self, chrom,chromlen, startpostocollecteSNP, endpostocollectSNP, idxedreffilehandler,ancestralgenomenameaddtotable, refindex, flanklen,outfile, tablename="derived_alle_ref"):
#         testfile=open("testsnpfile.txt",'a')
#         self.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(self.dbvariant_name, tablename, "context", "char(3)", "default null"))
        #add a temp function code block to process the snp before startpos of collecteSNP
        tempsnps=self.dbvariant.operateDB("select","select * from " + tablename + " where chrID='" + chrom + "' and snp_pos>= 2 and snp_pos<=" + str(startpostocollecteSNP))
        RefSeqMaptemp = Util.getRefSeqBypos_faster(idxedreffilehandler, refindex, chrom, 1, startpostocollecteSNP+1,chromlen)
        for snp in tempsnps:
            self.dbvariant.operateDB("update","update "+tablename+" set context='"+''.join(RefSeqMaptemp[chrom][snp[1]-1:snp[1]+2])+"' where chrID='"+ chrom + "' and snp_pos= "+str(snp[1]))
            l=copy.deepcopy(RefSeqMaptemp[chrom][snp[1]-1:])
            l[1]="N"
            print(">"+chrom+"_"+str(snp[1])+"\n"+"".join(l), file=outfile)
            
        #temp function code block end
        snps = self.dbvariant.operateDB("select", "select * from " + tablename + " where chrID='" + chrom + "' and snp_pos>= " + str(startpostocollecteSNP) + " and snp_pos<=" + str(endpostocollectSNP))
        RefSeqMap = Util.getRefSeqBypos_faster(idxedreffilehandler, refindex, chrom, startpostocollecteSNP-flanklen, endpostocollectSNP+flanklen,chromlen)
        
        for snp in snps:
            currentsnpPos = snp[1]
            if len(snp[3]) != 1 :
        #                        print(snp[4])
                continue# skip indel
            currentsnpID=chrom+"_"+str(snp[1])
            if currentsnpPos + flanklen <= RefSeqMap[chrom][0] + len(RefSeqMap[chrom]) - 1 and currentsnpPos - flanklen > RefSeqMap[chrom][0] :
                snpflankseq = ''.join(RefSeqMap[chrom][(currentsnpPos - flanklen - RefSeqMap[chrom][0]):(currentsnpPos + flanklen - RefSeqMap[chrom][0] + 1)])
                self.dbvariant.operateDB("update","update "+tablename+" set context='"+snpflankseq[flanklen-1:flanklen+2]+"' where chrID='"+ chrom + "' and snp_pos= "+str(currentsnpPos))
                snpflankseq=snpflankseq[0:flanklen]+'N'+snpflankseq[flanklen+1:]
            elif currentsnpPos <= RefSeqMap[chrom][0] + len(RefSeqMap[chrom]) - 1 and currentsnpPos - flanklen > RefSeqMap[chrom][0]:
                snpflankseq = ''.join(RefSeqMap[chrom][(currentsnpPos - flanklen - RefSeqMap[chrom][0]):(currentsnpPos - RefSeqMap[chrom][0] + 1)])
#                 print(currentsnpID,snpflankseq[flanklen],file=testfile)
                self.dbvariant.operateDB("update","update "+tablename+" set context='"+snpflankseq[flanklen-1:flanklen+1]+"N' where chrID='"+ chrom + "' and snp_pos= "+str(currentsnpPos))
                snpflankseq=snpflankseq[0:flanklen]+'N'
            elif currentsnpPos - flanklen <= RefSeqMap[chrom][0]:
                snpflankseq = ''.join(RefSeqMap[chrom][(currentsnpPos - RefSeqMap[chrom][0]):(currentsnpPos + flanklen - RefSeqMap[chrom][0] + 1)])
#                 print(currentsnpID,snpflankseq[0],file=testfile)
                self.dbvariant.operateDB("update","update "+tablename+" set context='N"+snpflankseq[0:2]+"' where chrID='"+ chrom + "' and snp_pos= "+str(currentsnpPos))
                snpflankseq = 'N'+snpflankseq[1:flanklen+1]
                
            else:
                print("what's wrong with the func getflankseqs ?")
                exit(-1)
#            if currentsnpPos + 25 <= RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1 and currentsnpPos - 25 > RefSeqMap[lastchromNo][0] :
#            snpflankseq = ''.join(RefSeqMap[chrom][(currentsnpPos - 25 - RefSeqMap[chrom][0]):(currentsnpPos + 25 - RefSeqMap[chrom][0] + 1)])
#            print(currentsnpID, snpflankseq[25], file=testfile)
#             snpflankseq = snpflankseq[0:25] + 'N' + snpflankseq[26:]
#             print(">" + currentsnpID + "\n" + snpflankseq, end='\n', file=outfile)
#         testfile.close()
        #                    print("update "+finaltable+" set fafilepos="+str(filepos)+" where snpID='"+currentsnpID+"'")
    def getregionseqstooutfile(self, chrom, chromlen, startpostocollecteREGIONs, endpostocollectREGIONs, idxedreffilehandler, ancestralgenomenameaddtotable, refindex, minRegionLEN, outfile, regionsOfOneChrom, tablename):
        RefSeqMap = Util.getRefSeqBypos_faster(idxedreffilehandler, refindex, chrom, startpostocollecteREGIONs-minRegionLEN, endpostocollectREGIONs+minRegionLEN,chromlen)
        for region in regionsOfOneChrom:
            regionStart=int(region[1])
            regionEnd=int(region[2])
            if minRegionLEN<(regionEnd-regionStart):
                c=minRegionLEN-(regionEnd-regionStart)
                if regionStart -c>RefSeqMap[chrom][0] :
                    regionStart=regionStart -c
                elif regionEnd+c<RefSeqMap[chrom][0] + len(RefSeqMap[chrom]) - 1 :
                    regionEnd=regionEnd+c
            if regionEnd <= RefSeqMap[chrom][0] + len(RefSeqMap[chrom]) - 1 and regionStart > RefSeqMap[chrom][0] :
                snpflankseq = ''.join(RefSeqMap[chrom][(regionStart - RefSeqMap[chrom][0]):(regionEnd - RefSeqMap[chrom][0] + 1)])
                currentsnpID=chrom+"_"+str(regionStart)+"_"+str(regionEnd)
                print(">" + currentsnpID + "\n" + snpflankseq, end='\n', file=outfile)
    def callblast(self,pathtoblastn,pathtorefdb,queryfaFile,BlastOutFile):
        #outfmt chose 6 suggest by zhaoyiqiang
        shellstatment=pathtoblastn+" -query "+queryfaFile+" -task blastn -db "+pathtorefdb+" -out "+BlastOutFile +" -outfmt 7 -num_threads 8"
        print(shellstatment)
        a = os.system(shellstatment)
        if a != 0:
            print("Ancestralallele.py : callblast func os.system return not 0")
            exit(-1)
        print(shellstatment,a,"OK")
    def extarctBlastOut(self,BlastOutFile,flanklen=25):
        print(" query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score")
        a = os.popen("awk '$1!~/^#/ && $4==121 && $5==1  {print $0}' " + BlastOutFile)
        posfile=open(BlastOutFile+".pos",'w')
    
        lastbasesAccur = {}
        onegroup=[]
        revcom="+"
    #    initial
        hit = a.readline()
        hitlist = re.split(r"\s+", hit)
    
        sendpos = int(hitlist[9])
        sstartpos = int(hitlist[8])
        qstartpos = int(hitlist[6])
        blastlen=int(hitlist[3])
        
        snpindex = flanklen+1 - qstartpos
        if sstartpos > sendpos:
            temp = sstartpos
            sstartpos = sendpos
            sendpos = temp
            revcom="-"
        snp_loc_s=sstartpos+flanklen
        lastsnpID = hitlist[0]
        chrom= hitlist[1]
#         RefSeqMap = Util.getRefSeqBypos(refFastahander=ancestryreffile, refindex=ancestryrefidx, currentChromNO=chrom, startpos=sstartpos, endpos=sendpos)
#         if revcom:
#             tempStr=RefSeqMap[chrom][1:]
#             tempStr.reverse()
#             RefSeqMap[chrom][1:]=Util.complementary(tempStr)
#             revcom=False
            
#         lastbasesAccur[RefSeqMap[chrom][snpindex + 1]] = [(chrom, sstartpos, sendpos)]
        onegroup.append((lastsnpID,blastlen,chrom, sstartpos, sendpos))
        for hit in a:
            print(hit)
            hitlist = re.split(r"\s+", hit)
            chrom = hitlist[1]
            sstartpos = int(hitlist[8])
            sendpos = int(hitlist[9])
            qstartpos = int(hitlist[6])
            blastlen=int(hitlist[3])
            snp_loc_s=sstartpos+flanklen
            snpindex = flanklen+1 - qstartpos
            if sstartpos > sendpos:
                temp = sstartpos
                sstartpos = sendpos
                sendpos = temp
                revcom="-"
            else:
                revcom="+"
            if lastsnpID == hitlist[0]:

                onegroup.append((lastsnpID,blastlen,chrom, sstartpos, sendpos))
            else:
                snppos=re.search(r"_(\d+)",lastsnpID).group(1)
                snpChrom=re.search(r"(.+)_(\d+)",lastsnpID).group(1)
                onegroup.sort(key=lambda listRec:listRec[1])                           
                if len(onegroup)==1 or onegroup[0][1]-onegroup[1][1]>=15:#first , only one query id,second longest hit 15 bases greater than the second longest hit
                    chr_pos_info_base=re.split(r"_",onegroup[0][0])
                    print(*onegroup[0][2:],chr_pos_info_base[0],chr_pos_info_base[1],revcom,onegroup[0][1],chr_pos_info_base[2],file=posfile)
#                 RefSeqMap = Util.getRefSeqBypos(refFastahander=ancestryreffile, refindex=ancestryrefidx, currentChromNO=chrom, startpos=sstartpos, endpos=sendpos)
#                 if revcom:
#                     tempStr=RefSeqMap[chrom][1:]
#                     tempStr.reverse()
#                     RefSeqMap[chrom][1:]=Util.complementary(tempStr)
#                     revcom=False            
#                 print(hitlist[0],RefSeqMap[chrom][snpindex + 1],str(snp_loc_s),"".join(RefSeqMap[chrom][1:]),file=ancestrysnpflank)
    #            dbvariant.operateDB("update", "update " + finaltable + " set chicken='" + RefSeqMap[chrom][snpindex + 1] + "' where snpID='" + hitlist[0] + "'")
                lastsnpID = hitlist[0]
                onegroup=[(lastsnpID,blastlen,chrom, sstartpos, sendpos)]
        posfile.close()
        print("finish")

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
                onegroup=[(RefSeqMap[chrom][snpindex + 1],blastlen)]
                lastbasesAccur.clear()
                lastbasesAccur[RefSeqMap[chrom][snpindex + 1]] = [(chrom, sstartpos, sendpos)]
        print("finish")
        ancestryreffile.close()
    def fillAncestral(self,archicpopVcfFile,chromlist,toplevelsnptablename="ducksnp_toplevel"):
        """
        abandon the snps which exist in archicpopVcfFile but absence in all others pop snp sets 
        """
        

        for vcffilename in archicpopVcfFile.keys():
#             archicpop = VCFutil.VCF_Data(archicpopVcfFile)
            archicpop=archicpopVcfFile[vcffilename][0]
            archicpop_colname=re.search(r'[^/]*$',vcffilename).group(0)
            archicpop_colname=re.sub(r"[^\w^\d]","_",archicpop_colname)            
            print(archicpop_colname+"_alt",archicpop_colname+"_dep")
            self.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(self.dbvariant_name, toplevelsnptablename, archicpop_colname+"_alt", "char(128)", "default null"))
            self.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(self.dbvariant_name, toplevelsnptablename, archicpop_colname+"_dep", "char(128)", "default null"))  
            print("callproc", "mysql_sp_add_column","done",time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
#             depthfile = Util.GATK_depthfile(depthFile, depthFile + ".index")
#             print("Util.GATK_depthfile done",time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
#             species_idx_list=[]
#     #         species_idx = depthfile.title.index("Depth_for_" + archicpopNameindepthFile)
#             for i in range(0,len(depthfile.title)):
#                 if re.search(r""+archicpopNameindepthFile,depthfile.title[i])!=None:
#                     species_idx_list.append(i)        
#     后面程序需要调一下，考虑一下indvd 的 AF DP仿照Caculators里面的东西
    #         totalChroms = self.dbgenome.operateDB("select","select count(*) from "+chromtable)[0][0]
    #         for i in range(0,totalChroms,20):
    #             currentsql="select * from " + chromtable+" order by chrlength desc limit "+str(i)+",20"
    #             result=self.dbgenome.operateDB("select",currentsql)
            lastposofdepthfilefp=0
            for currentchrID in chromlist:
                print(currentchrID+":",end="")
            
                archicpopSeqOfAChr={}
                archicpopSeqOfAChr[currentchrID]=archicpop.getVcfListByChrom(currentchrID)
                allsnpsInAchr=self.dbvariant.operateDB("select","select snp_pos,alt_base from "+toplevelsnptablename+" where chrID='"+currentchrID+"'")
                updatesql_statments=[]
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
                            pos, REF, archicpop_ALT, INFO,FORMAT,samples = archicpopSeqOfAChr[currentchrID][mid]
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
                            popsdata_alt=archicpop_ALT
                            popsdata_dep=str(refdep)+","+str(altalleledep)
                            break
                    else:
                        sum_depth=0
                        for samfile in archicpopVcfFile[vcffilename][1:]:
                            ACGTdep=samfile.count_coverage(currentchrID,snp_pos-1,snp_pos)
                        for dep in ACGTdep:
                            sum_depth+=dep[0]
                        if sum_depth<4:
                            popsdata_alt=ALT
                            popsdata_dep="no covered"
                        else:
                            popsdata_alt=ALT
                            popsdata_dep=str(sum_depth) + ",0"

                    updatesql_statments.append("update " + toplevelsnptablename + " set "+archicpop_colname+"_alt = '" + popsdata_alt+"',"+archicpop_colname+"_dep= '"+popsdata_dep+"' where chrID="+"'"+currentchrID+"' and snp_pos="+str(snp[0]))
                
                if  updatesql_statments:
                    print("updatesql_statments",len(updatesql_statments))
                    for update_statments in updatesql_statments:
                        self.session.execute(update_statments)
                        self.session.commit()
#                     self.dbvariant.operateDB("update", *updatesql_statments)
    def leftjoinSelectedTables(self,chromlist,outtable_file_Name,depthfilenames,vcftables=[],toplevelsnptable="ducksnp_toplevel",drop=False):
        depthobjmap={};lastposofdepthfilefp={}#
        for vcftablename in depthfilenames.keys():
            if depthfilenames[vcftablename]!=None:
                lastposofdepthfilefp[vcftablename]=0
                depthobjmap[vcftablename]=Util.GATK_depthfile(depthfilenames[vcftablename][0], depthfilenames[vcftablename][0] + ".index")
            else:
                depthobjmap[vcftablename]=None

        species_idx_map={}
        
        for vcftable in vcftables:
            species_idx_map[vcftable]=[]
            if depthfilenames[vcftablename]!=None:
                for name in depthfilenames[vcftable][1:]:
                    for i in range(0,len(depthobjmap[vcftable].title)):
                        if re.search(r""+name,depthobjmap[vcftable].title[i])!=None:
                            species_idx_map[vcftable].append(i)
#             popnamesubStringInDepthTitle=re.search(r"^([^_]+)_.*",vcftable).group(1)
#             species_idx_map[popnamesubStringInDepthTitle]=[]
#             for depthfilehanlder in depthfilehandlerlist:
#                 for i in range(0,len(depthfilehanlder.title)):
#                     if re.search(r""+popnamesubStringInDepthTitle,depthfilehanlder.title[i])!=None:
#                         species_idx_map[popnamesubStringInDepthTitle].append(i)   
        
        outfile=open(outtable_file_Name,'w')
        outtable_Name=re.search(r'[^/]*$', outtable_file_Name).group(0)
        outtable_Name=re.sub(r"[^\w^\d]","_",outtable_Name)
        outtable_titlelist=[a[0].strip() for a in self.dbvariant.operateDB("select","select column_name  from information_schema.columns where table_schema='"+self.dbvariant_name+"' and table_name='"+toplevelsnptable+"'")]  
        toplevellen=len(outtable_titlelist)
        outtable_titlelist=outtable_titlelist+vcftables
        print(*outtable_titlelist,sep="\t",file=outfile)
        for currentchrID in chromlist: 
            sqlselectstatementpart="select t.*"
            sqlfromstatementpart=" from "+toplevelsnptable+" as t "
            #sql statement produce part1
            for vcftable in vcftables:
                sqlselectstatementpart=sqlselectstatementpart+","+vcftable.strip()+".alt_base as "+vcftable+"_alt_base,"+vcftable.strip()+".AF as "+vcftable+"_AF"
            
            #sql statement produce part2
            for vcftable in vcftables:
                sqlfromstatementpart=sqlfromstatementpart+" left join "+vcftable.strip()+" using(chrID,snp_pos)"
            #sql where statement append
            sqlstatement=sqlselectstatementpart+sqlfromstatementpart+" where chrID='"+currentchrID+"'"
            print(sqlstatement)
            allsnpOfJoinTableinAchr=self.dbvariant.operateDB("select",sqlstatement)
            #process value,merge into one col
            NumOfColOftoplevel_fix=len(self.dbvariant.operateDB("select","select column_name  from information_schema.columns where table_schema='"+self.dbvariant_name+"' and table_name='"+toplevelsnptable+"'"))
            #               print to outfile and then  fix the depth to None value
            for rec in allsnpOfJoinTableinAchr:
                snp_pos=rec[1]
                NumOfColOftoplevel=NumOfColOftoplevel_fix
                recToPrint=list(rec[0:NumOfColOftoplevel])
                for vcftable in vcftables:
                    if rec[NumOfColOftoplevel]==None:
                        if depthfilenames[vcftable]==None:
                            popsdata="unknow" #may be fixed as ref ,check the depth file ..... need to be done
                        else:
                            sum_depth=0
                            #find the depth file which the title contain the popnamesubStringInDepthTitle
                            depth_linelist = depthobjmap[vcftable].getdepthByPos_optimized(currentchrID, snp_pos)
#                             lastposofdepthfilefp[vcftable]=depthobjmap[vcftable].depthfilefp.tell()
#                             popnamesubStringInDepthTitle=re.search(r"^([^_]+)_.*",vcftable).group(1).strip()
                            
                            # accumulate  the depth
                            for idx in species_idx_map[vcftable]:
                                sum_depth+=int(depth_linelist[idx])
                            if sum_depth<20:
                                popsdata="no covered"
                            else:
                                ALT=rec[4]
                                popsdata=ALT+":0"
                    else:
                        ALT=rec[NumOfColOftoplevel]
                        AF=rec[NumOfColOftoplevel+1]
                        popsdata=ALT+":"+str(AF)
            
                    recToPrint=recToPrint+[popsdata]
                    NumOfColOftoplevel=NumOfColOftoplevel+2   #alt_base and AF col
                print(*recToPrint,sep="\t",file=outfile)
        if drop:
            try:
                self.dbtmp.drop_table(outtable_Name)
            except:
                print("except in drop table", outtable_Name)
        self.dbtmp.operateDB("copytableschema","create table "+outtable_Name+" like "+self.dbvariant_name.strip()+"."+toplevelsnptable.strip())
        for colname in outtable_titlelist[toplevellen-1:]:
            self.dbtmp.operateDB("callproc", "mysql_sp_add_column",data=(self.dbvariant_name, outtable_Name, colname, "char(128)", "default null"))
        loaddatasql="load data local infile '"+outtable_file_Name+"' into table "+outtable_Name+" fields terminated by '\\t'"
        print(loaddatasql)
        outfile.close()
#         a=os.system("mysql -uroot -p1234567 -D" + self.dbtmpname.strip() + ' -e "' + loaddatasql + '"')
#         if a!=0:
#             print(a,"maybe error")
#         else:
#             print("finish")

                        
                