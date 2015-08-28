# -*- coding: UTF-8 -*-

'''
Created on 2014-6-24

@author: liurui
'''

import src.NGS.BasicUtil.DBManager as dbm
from optparse import OptionParser
import pickle
import re

from NGS.BasicUtil import Util, VCFutil

tempdbname="temp"
parser = OptionParser()
parser.add_option("-c", "--topleveltable", dest="topleveltable",# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="write report to FILE")
parser.add_option("-f", "--farsurebutfew", dest="farsurebutfew",help="far sure but with few locs")
parser.add_option("-o","--outfileprename",dest="outfilepreName",help="outfilepreName with path")
parser.add_option("-a","--ancestralspeciescolname",dest="ancestralspeciescolname",help="ancestralspecisname")
parser.add_option("-m","--mindepth",dest="mindepth",help="mindepth for both archicpop and ancestralallel")#
# parser.add_option("-e", "--bedfile", action="append",dest="bedfile",default=[],help="measure region of bedfile")
parser.add_option("-s","--speciesesName",action="append",dest="speciesesName",default=[],help="speciesName in table")#
                                                                                                                                                          
(options, args) = parser.parse_args()
allpop_with_derived_alletable=options.topleveltable
ancestralspeciescolname=options.ancestralspeciesname.strip()
farsurebutfew=options.farsurebutfew.strip()
mindepth=int(options.mindepth)
if __name__ == '__main__':
    dbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    tableprename=""
    TABLES = {}
    for bedfileName in args[:]:

        tableprename+=re.search(r"[^/]*$",bedfileName).group(0)[0]
    tablename=tableprename+Util.random_str()
    TABLES[tablename] = (
        "CREATE TABLE " + tablename + " ("
        " `snpID` varchar(128) NOT NULL ,"
        " `region` varchar(128) NOT NULL ,"
        " `DAF` double default 100 ,"
        " `MAF` double default 100 ,"
        " PRIMARY KEY (`snpID`,`region`)"
        ")"
        )
    
    tempdbtools=dbm.DBTools(Util.ip, Util.username, Util.password, Util.ghostdbname)
    tempdbtools.create_table(TABLES)
    titlelist=[a[0].strip() for a in dbtools.operateDB("select","select column_name  from information_schema.columns where table_schema='"+Util.vcfdbname+"' and table_name='"+options.topleveltable+"'")]
    ancestralspeciesidx=titlelist.index(ancestralspeciescolname)
    farsurebutfewidx=titlelist.index(farsurebutfew)
    print(titlelist)
    for speciesName in options.speciesesName:
        outfile=open(options.outfilepreName+speciesName,'w')
        print("snpID\tregion\tDAF\tMAF\tgroup",file=outfile)
        speciesidx=titlelist.index(speciesName)
        
        for bedfileName in args[:]:
            
            regionmap=Util.bedfiletools(bedfileName)
            regionName=re.search(r"[^/]*$",bedfileName).group(0).replace('.','_')[0:10]+str(args.index(bedfileName)+1)
#             tempdbtools.operateDB("callproc", "mysql_sp_add_column", data=(options.dbname, tablename, regionName+"DAF", "double", "default 100"))
#             tempdbtools.operateDB("callproc", "mysql_sp_add_column", data=(options.dbname, tablename, regionName+"MAF", "double", "default 100"))

            for chrom in regionmap:
                for startpos,endpos,optionalfields in regionmap[chrom]:
                    print("select * from "+allpop_with_derived_alletable + " where chrID='"+chrom+"' and  snp_pos>="+str(startpos)+" and snp_pos<="+str(endpos))
                    selectedsnps = dbtools.operateDB("select","select * from "+allpop_with_derived_alletable + " where chrID='"+chrom+"' and snp_pos>="+str(startpos)+" and snp_pos<="+str(endpos))
                    if selectedsnps==None:
                        continue
                    for snp in selectedsnps:
                        if snp[speciesidx]==None or snp[speciesidx].strip()=="no covered" or snp[speciesidx]=="no covered" or re.search(r'[\w\W]+[,][\w\W]+:\d+,\d+',snp[4])!=None:
                            continue
                        snpid=snp[0].strip()+"_"+str(snp[1])
                        archicpop=re.search(r'([ATCGatcg]+):(\d+),(\d+)',snp[speciesidx])
                        archic_base=archicpop.group(1).strip().upper()
#                         archicsnp = re.search(r"(\w*):(\d*),(\d*)", snp[speciesidx])
#                         archicaltalle=archicsnp.group(1).strip().upper()
                        archicrefcount=int(archicpop.group(2).strip())
                        archicaltcount=int(archicpop.group(3).strip())
#                         ancestralallel=snp[farsurebutfewidx].strip().upper()
                        ref_base=snp[3].strip().upper()
#                         if snp[farsurebutfewidx]!=None:
#                             if snp[farsurebutfewidx].strip().upper()!=archicaltalle:
#                                 print("skip:",snp)
#                                 continue
                        if (archicrefcount!=0 and archicaltcount!=0) or (archicrefcount+archicaltcount)<mindepth:
                            print("skip:",snp)
                            continue
                        popsnpcount=re.split(r',',snp[speciesidx])
                        if int(popsnpcount[0].strip())+int(popsnpcount[1].strip())<mindepth:
                            continue
                        if archicpop.group(2).strip()=='0':
                            if snp[farsurebutfewidx]!=None and archic_base!= snp[farsurebutfewidx].strip().upper():
                                continue
                            DAF=int(popsnpcount[1])/(int(popsnpcount[0])+int(popsnpcount[1]))#alt_allele is the ancestral allele
                        elif archicpop.group(3).strip()=='0':
                            if snp[farsurebutfewidx]!=None and snp[farsurebutfewidx].strip().upper()!=ref_base:
                                continue
                            DAF=int(popsnpcount[0])/(int(popsnpcount[0])+int(popsnpcount[1]))#ref_allele is the ancestral allele
#                         if archicaltalle==snp[3].strip().upper():#ancestral allel is the same with ref allele type
#                             DAF=int(popsnpcount[0])/(int(popsnpcount[0])+int(popsnpcount[1]))
#                         elif archicaltalle==snp[4].strip().upper():#ancestral allel is the same with alt allele type
#                             DAF=int(popsnpcount[1])/(int(popsnpcount[0])+int(popsnpcount[1]))
#                         else:
#                             continue
                        MAF=min(int(popsnpcount[0]),int(popsnpcount[1]))/(int(popsnpcount[0])+int(popsnpcount[1]))
                        print("insert into "+tablename+"(snpID,region,DAF,MAF) values(%s,%s,%s,%s)")
                        tempdbtools.operateDB("insert","insert into "+tablename+"(snpID,region,DAF,MAF) values(%s,%s,%s,%s)",data=(snpid,regionName,DAF,MAF))
                    else:
                        print("no snp in " ,bedfileName,chrom,regionmap[chrom])
        totalsnps = tempdbtools.operateDB("select", "select count(*) from " + tablename)[0][0]
        if totalsnps==0:
            break

        for i in range(0,totalsnps,1000):
            snps = tempdbtools.operateDB("select","select * from "+tablename+" limit "+str(i)+",1000")
            for snp in snps:   
                print(snp[0],str(snp[2]),str(snp[3]),snp[1],file=outfile)
                
#         tempdbtools.drop_table(tablename)
        print(args[:])
        outfile.close()

