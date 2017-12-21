'''
Created on 2017-12-19

@author: Dr.liu
'''
from optparse import OptionParser
import random
import re

from src.NGS.BasicUtil import Util, VCFutil
import src.NGS.BasicUtil.DBManager as dbm
from src.NGS.Service import Ancestralallele


parser = OptionParser()
parser.add_option("-t","--toplevelsnptable",dest="toplevelsnptable",default="ducksnp_toplevel",help="depth of the folder to output")
parser.add_option("-m","--minlength",dest="minlength",help="require least chrom length")
parser.add_option("-q","--quantizing",dest="quantizing",action="append",help="select from some first param of -v ")
parser.add_option("-n","--noofindvds2quantizing",dest="noofindvds2quantizing",help="together with -q ")
parser.add_option("-C","--countsnpnumberfromvcf",dest="countsnpnumberfromvcf",default=None)
parser.add_option("-d","--snpperkb",dest="snpperkb")
parser.add_option("-c","--chromlist",dest="chromlist")
parser.add_option("-o","--outputfilename",dest="outputfilename")
parser.add_option("-v", "--vcffile", dest="vcffile",action="append",default=[],nargs=2,help="")
(options, args) = parser.parse_args()

minlength=options.minlength;toplevelsnptable=options.toplevelsnptable;snpperkb=float(options.snpperkb);vcffilelist=options.vcffile
outgroupidx_in_topleveltable=[6,8];minoutgroupdepth=30

noofindvds2quantizing=int(options.noofindvds2quantizing)
dadisnpfile=open(options.outputfilename+"dilutetodensity"+options.snpperkb.strip(),'w')
dbvariantstools=dbm.DBTools(Util.ip, Util.username,Util.password, Util.vcfdbname)
dynamicIU_toptable_obj=Ancestralallele.dynamicInsertUpdateAncestralContext(dbvariantstools,Util.beijingreffa,options.topleveltablejudgeancestral)

flankseqfafile=open(options.outputfilename+re.search(r"[^/]*$",options.chromlist).group(0)+".fa","w")
if __name__ == '__main__':
    chromlistfile=open(options.chromlistfilename,"r")
    selectedchroms=[]
    for chrrow in chromlistfile:
        chrrowlist=re.split(r'\s+',chrrow.strip())
        selectedchroms.append((chrrowlist[0].strip(),int(chrrowlist[1].strip())))
    chromlistfile.close()
    ###########title print   and produce vcfobjlist############
    vcfobjmap={}
    
    toplevelsnptable_titlelist=[a[0].strip() for a in dbvariantstools.operateDB("select", "select column_name  from information_schema.columns where table_schema='" + "ninglabvariantdata" + "' and table_name='" + toplevelsnptable + "'")]
    
    print(Util.pekingduckchromtable[:9],toplevelsnptable_titlelist[outgroupidx_in_topleveltable[0]][:8],"Allele1",sep="\t",end="\t",file=dadisnpfile)
    vcftablesidxlist_toquantizing=[]
    idx=0
    for vcftfile_name,minAN in vcffilelist:
        vcfobjmap[vcftfile_name]=VCFutil.VCF_Data(vcftfile_name)
        if options.quantizing!=None:
            for quantizingvcftable in options.quantizing:
                if quantizingvcftable.strip()==vcftfile_name.strip():
                    vcftablesidxlist_toquantizing.append(idx)
        idx+=1
        popName=re.split(r'_',vcftfile_name)[0]
        print(popName,end="\t",file=dadisnpfile)
    print("vcftablesidxlist_toquantizing",vcftablesidxlist_toquantizing)
    if options.quantizing!=None:
        print("quantizpool",end="\t",file=dadisnpfile)
    print("Allele2",end="\t",file=dadisnpfile)
    for vcftable_name,minAN in vcffilelist:
        popName=re.split(r'_',vcftable_name)[0]
        print(popName,end="\t",file=dadisnpfile)
    if options.quantizing!=None:
        print("quantizpool",end="\t",file=dadisnpfile)
    print("Gene\tPosition",file=dadisnpfile)
    ############ #########              finish title print                       ##################################
    totalsnp=0;totallength=0;totallengthduilt=0;
    for currentchrID,currentchrLen in selectedchroms:
        dynamicIU_toptable_obj.currentchrLen=currentchrLen
        listOfpopvcfRecsMapByChr=[]
        ####################        produce input to Util.alinmultPopSnpPos  ##########################################
        for vcfname,minAN in vcffilelist:
            vcflistOfAchr=vcfobjmap[vcfname].getVcfListByChrom(currentchrID,MQfilter=None)
            listOfpopvcfRecsMapByChr.append({currentchrID:vcflistOfAchr})
            ###########################     produce input end        ##################################################
        fulloutjoinSNPs=Util.alinmultPopSnpPos(listOfpopvcfRecsMapByChr, "o")
        totalsnpforAchr=len(fulloutjoinSNPs[currentchrID])
        dilute = snpperkb*currentchrLen / (1000 * totalsnpforAchr)
        #use currentchrLen because almost all region has been sequenced
        if totalsnpforAchr>=snpperkb*currentchrLen / 1000:
            totallength+=currentchrLen
            totalsnp+=totalsnpforAchr
            totallengthduilt+=dilute*totalsnpforAchr
            sample_idxlistOfaJoinTable=random.sample([j for j in range(totalsnpforAchr)],int(dilute*totalsnpforAchr)+1)
            sample_idxlistOfaJoinTable.sort()
        else:
            print("skip this chrom",currentchrID)
            continue
        ############        filter MQ , minAN and ancestral info #######################################################################
        for sampled_idx in sample_idxlistOfaJoinTable:
            snp_pos=fulloutjoinSNPs[currentchrID][sampled_idx][0];REF=fulloutjoinSNPs[currentchrID][sampled_idx][1];ALT=fulloutjoinSNPs[currentchrID][sampled_idx][2]
            continuesearch=-1;sampled_idx_find_satisfied=sampled_idx;direction=1
            while continuesearch==-1:# -1 continuesearch; 1 OUTgroup1 passed; 2 secondgroupbase
                snprec_in_toplevel=dbvariantstools.operateDB("select","select * from "+toplevelsnptable+" where chrID='"+currentchrID+"' and snp_pos='"+str(fulloutjoinSNPs[currentchrID][sampled_idx_find_satisfied][0])+"'")
                if sampled_idx_find_satisfied==len(fulloutjoinSNPs[currentchrID]) or (sampled_idx!=sample_idxlistOfaJoinTable[-1] and sampled_idx_find_satisfied==sample_idxlistOfaJoinTable[sample_idxlistOfaJoinTable.index(sampled_idx)+1] ) or sampled_idx_find_satisfied==-1 or (sampled_idx!=sample_idxlistOfaJoinTable[0] and sampled_idx_find_satisfied==sample_idxlistOfaJoinTable[sample_idxlistOfaJoinTable.index(sampled_idx)-1]):
                    if direction==-1:
                        print("search snp out of range ",currentchrID,snp_pos)
                        totalsnp-=1
                        break
                    direction=-1
                    print("direction changed",direction)
                    sampled_idx_find_satisfied=sampled_idx#start again, but the opposit deriction
                elif (not snprec_in_toplevel or snprec_in_toplevel[0][outgroupidx_in_topleveltable[0]]==None):
                    dynamicIU_toptable_obj.insertorUpdatetopleveltable({currentchrID:[(snp_pos,REF,ALT)]}, flankseqfafile, 5)
                    continue#shoud excute only once for one position
                elif snprec_in_toplevel and snprec_in_toplevel[0][outgroupidx_in_topleveltable[0]]!=None:
                    OUTgroup1=re.split(r",",fulloutjoinSNPs[currentchrID][sampled_idx_find_satisfied][outgroupidx_in_topleveltable[0]+1])
                else:
                    print("warning ! this may be a except situation")
                    sampled_idx_find_satisfied+=direction
                snprec_in_toplevel=dbvariantstools.operateDB("select","select * from "+toplevelsnptable+" where chrID='"+currentchrID+"' and snp_pos='"+str(fulloutjoinSNPs[currentchrID][sampled_idx_find_satisfied][0])+"'")
                
                    
                    continue
                elif snprec_in_toplevel and (snprec_in_toplevel[0][outgroupidx_in_topleveltable[0]]==None or len(snprec_in_toplevel[0][5]))!=3:
                    sampled_idx_find_satisfied+=direction
                for pop in 
                OUTgroup1=re.split(r",",