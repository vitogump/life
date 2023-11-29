# -*- coding: UTF-8 -*-
'''
Created on 2015-5-7

@author: liurui
'''
import copy
from optparse import OptionParser
import random
import re

from src.NGS.BasicUtil import Util
import src.NGS.BasicUtil.DBManager as dbm


parser = OptionParser()
parser.add_option("-t","--toplevelsnptable",dest="toplevelsnptable",default="ducksnp_toplevel",help="depth of the folder to output")
parser.add_option("-m","--minlength",dest="minlength",help="require least chrom length")
parser.add_option("-A","--minAN",dest="minAN")
parser.add_option("-d","--snpperkb",dest="snpperkb")
parser.add_option("-o","--outputfilename",dest="outputfilename")
parser.add_option("-v", "--vcftablelist", dest="vcftablelist",action="append",default=[],help="")
(options, args) = parser.parse_args()
minlength=options.minlength;toplevelsnptable=options.toplevelsnptable;snpperkb=int(options.snpperkb);vcftableslist=options.vcftablelist;minAN=options.minAN
dadisnpfile=open(options.outputfilename,'w')
outgroupidx_in_topleveltable=6;minoutgroupdepth=30
if __name__ == '__main__':
    genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    dbvariantstools=dbm.DBTools(Util.ip, Util.username,Util.password, Util.vcfdbname)
    toplevelsnptable_titlelist=[a[0].strip() for a in dbvariantstools.operateDB("select", "select column_name  from information_schema.columns where table_schema='" + "ninglabvariantdata" + "' and table_name='" + toplevelsnptable + "'")]
    selectedchroms=genomedbtools.operateDB("select","select * from "+Util.pekingduckchromtable+" where chrlength>="+minlength)
    ######################## title print ##############################
    print(Util.pekingduckchromtable[:9],toplevelsnptable_titlelist[outgroupidx_in_topleveltable],"Allele1",sep="\t",end="\t",file=dadisnpfile)
    for vcftable_name in vcftableslist:
        popName=re.split(r'_',vcftable_name)[0]
        print(popName,end="\t",file=dadisnpfile)
    print("Allele2",end="\t",file=dadisnpfile)
    for vcftable_name in vcftableslist:
        popName=re.split(r'_',vcftable_name)[0]
        print(popName,end="\t",file=dadisnpfile)
    print("Gene\tPosition",file=dadisnpfile)
    ############               finish title print ##################################
    for row in selectedchroms:
        currentchrID=row[0]
        currentchrLen=int(row[1])
        sqlselectstatementpart="select t.*"
#         sqlfromstatementpart=" from "+toplevelsnptable+" as t "
        for vcftable in vcftableslist:
            sqlselectstatementpart=sqlselectstatementpart+","+vcftable.strip()+".AC as "+vcftable+"_AC,"+vcftable.strip()+".AN as "+vcftable+"_AN"
        sqlselectstatementpart+=" from "+toplevelsnptable+" as t "
        for vcftable in vcftableslist:
            sqlselectstatementpart=sqlselectstatementpart+" inner join "+vcftable.strip()+" using(chrID,snp_pos)"
        sqlselectstatementpart=sqlselectstatementpart+" where chrID='"+currentchrID+"' and t.context is not null and length(t.ref_base)=length(t.alt_base) and (t."+toplevelsnptable_titlelist[outgroupidx_in_topleveltable] +" regexp ':0,' or t."+toplevelsnptable_titlelist[outgroupidx_in_topleveltable] +" regexp ',0$')"
        for vcftable in vcftableslist:
            sqlselectstatementpart=sqlselectstatementpart+" and "+vcftable+".AN >="+minAN
        print(sqlselectstatementpart)
        allsnpOfJoinTableinAchr=dbvariantstools.operateDB("select",sqlselectstatementpart)
        if len(allsnpOfJoinTableinAchr)>int(currentchrLen/snpperkb)*2:
            allsnpOfJoinTableinAchr_sampled_idxlist=random.sample([j for j in range(len(allsnpOfJoinTableinAchr))],int(currentchrLen/snpperkb))
            allsnpOfJoinTableinAchr_sampled_idxlist.sorted()
        else:
            allsnpOfJoinTableinAchr_sampled_idxlist=[j for j in range(len(allsnpOfJoinTableinAchr))]
        for sampled_idx in allsnpOfJoinTableinAchr_sampled_idxlist:
            contextwithinspeces=allsnpOfJoinTableinAchr[sampled_idx][5].upper()
            contextoutgroup=copy.copy(contextwithinspeces)
            outgroup=allsnpOfJoinTableinAchr[sampled_idx][6];outgroupBase=re.split(r':',outgroup);outgroupdepth=re.split(r",",outgroupBase[1]);outgroupBase=outgroupBase[0].upper()
            ALT=allsnpOfJoinTableinAchr[sampled_idx][4].upper().strip()
            sampled_idx_find_satisfied=sampled_idx
            #find the nearest correct rec
            while (int(outgroupdepth[0])<=minoutgroupdepth and int(outgroupdepth[1])<=minoutgroupdepth) or outgroupBase!=ALT:
                sampled_idx_find_satisfied+=1
                if sampled_idx_find_satisfied==allsnpOfJoinTableinAchr_sampled_idxlist[-1] or sampled_idx_find_satisfied==len(allsnpOfJoinTableinAchr) or sampled_idx_find_satisfied==allsnpOfJoinTableinAchr_sampled_idxlist[allsnpOfJoinTableinAchr_sampled_idxlist.index(sampled_idx)+1] :
                    break
                
                outgroup=allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][6];outgroupBase=re.split(r':',outgroup);outgroupdepth=re.split(r",",outgroupBase[1]);outgroupBase=outgroupBase[0].upper()
                ALT=allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][4].upper().strip()
            else:#find
                contextwithinspeces=allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][5].upper()
                contextoutgroup=copy.copy(contextwithinspeces)
                if int(outgroupdepth[0])==0:
                    contextoutgroup=contextoutgroup[0]+ALT+contextoutgroup[2]
                REF=allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][3].upper()
                postion=allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][1]
                print(contextwithinspeces,contextoutgroup,REF,sep="\t",end="\t",file=dadisnpfile)
                for vcftable_idx in range(len(vcftableslist)):
                    print(str(int(allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][7+vcftable_idx*2+1])-int(allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][7+vcftable_idx*2])),end="\t",file=dadisnpfile)
                print(ALT,end="\t",file=dadisnpfile)
                for vcftable_idx in range(len(vcftableslist)):
                    print(allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][7+vcftable_idx*2],end="\t",sep="\t",file=dadisnpfile)
                print(currentchrID.replace(".","_"),postion,sep="\t",file=dadisnpfile)


    dadisnpfile.close()
    print("finsih")