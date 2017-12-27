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


indvdNo_Of_POOL=10

parser = OptionParser()
parser.add_option("-t","--toplevelsnptable",dest="toplevelsnptable",default="ducksnp_toplevel",help="depth of the folder to output")
parser.add_option("-m","--minlength",dest="minlength",help="require least chrom length")
parser.add_option("-q","--quantizing",dest="quantizing",action="append",help="select from some first param of -v ")
parser.add_option("-n","--noofindvds2quantizing",dest="noofindvds2quantizing",help="together with -q ")

parser.add_option("-d","--snpperkb",dest="snpperkb")
parser.add_option("-c","--chromlist",dest="chromlist")
parser.add_option("-o","--outputfilename",dest="outputfilename")
parser.add_option("-v", "--vcffile", dest="vcffile",action="append",default=[],nargs=2,help="vcffile minAN")
(options, args) = parser.parse_args()

toplevelsnptable=options.toplevelsnptable;snpperkb=float(options.snpperkb);vcffilelist=options.vcffile#;minlength=options.minlength
outgroupidx_in_topleveltable=[6,8];minoutgroupdepth=30

noofindvds2quantizing=int(options.noofindvds2quantizing)
dadisnpfile=open(options.outputfilename+"dilutetodensity"+options.snpperkb.strip(),'w')
dbvariantstools=dbm.DBTools(Util.ip, Util.username,Util.password, Util.vcfdbname)
dynamicIU_toptable_obj=Ancestralallele.dynamicInsertUpdateAncestralContext(dbvariantstools,Util.beijingreffa,options.toplevelsnptable)

flankseqfafile=open(options.outputfilename+re.search(r"[^/]*$",options.chromlist).group(0)+".fa","a")
# recf=open("recf","w")
if __name__ == '__main__':
    chromlistfile=open(options.chromlist,"r")
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
        popName=re.split(r'\.',re.search(r"[^/]*$",vcftfile_name).group(0))[0]
        print(popName,end="\t",file=dadisnpfile)
    print("vcftablesidxlist_toquantizing",vcftablesidxlist_toquantizing)
    if options.quantizing!=None:
        print("quantizpool",end="\t",file=dadisnpfile)
    print("Allele2",end="\t",file=dadisnpfile)
    for vcftfile_name,minAN in vcffilelist:
        popName=re.split(r'\.',re.search(r"[^/]*$",vcftfile_name).group(0))[0]
        print(popName,end="\t",file=dadisnpfile)
    if options.quantizing!=None:
        print("quantizpool",end="\t",file=dadisnpfile)
    print("Gene\tPosition",file=dadisnpfile)
    ############ #########              finish title print                       ##################################
    totalsnp=0;totallength=0;totallengthduilt=0;totalduiltsnp=0
    for currentchrID,currentchrLen in selectedchroms:
        dynamicIU_toptable_obj.getRECsforCHR(currentchrID,currentchrLen)
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
            totallengthduilt+=dilute*currentchrLen
            sample_idxlistOfaJoinTable=random.sample([j for j in range(totalsnpforAchr)],int(dilute*totalsnpforAchr)+1)
            sample_idxlistOfaJoinTable.sort()
        else:
            print("skip this chrom",currentchrID)
            continue
        ############        filter MQ , minAN and ancestral info #######################################################################
        for sampled_idx in sample_idxlistOfaJoinTable:
            continuesearch=-1;sampled_idx_find_satisfied=sampled_idx;direction=1;A_base_idx=-1
            while continuesearch==-1:# -1 continuesearch; 1 OUTgroup1 passed; 2 secondgroupbase

                if sampled_idx_find_satisfied==len(fulloutjoinSNPs[currentchrID]) or (sampled_idx!=sample_idxlistOfaJoinTable[-1] and sampled_idx_find_satisfied==sample_idxlistOfaJoinTable[sample_idxlistOfaJoinTable.index(sampled_idx)+1] ) or sampled_idx_find_satisfied==-1 or (sampled_idx!=sample_idxlistOfaJoinTable[0] and sampled_idx_find_satisfied==sample_idxlistOfaJoinTable[sample_idxlistOfaJoinTable.index(sampled_idx)-1]):
                    if direction==-1:
                        print("search snp out of range around snppos",currentchrID,fulloutjoinSNPs[currentchrID][sampled_idx][0])
                        break
                    direction=-1
                    print("direction changed",direction)
                    sampled_idx_find_satisfied=sampled_idx+direction#start again, but the opposit deriction
                    continue
                snp_pos=fulloutjoinSNPs[currentchrID][sampled_idx_find_satisfied][0];REF=fulloutjoinSNPs[currentchrID][sampled_idx_find_satisfied][1];ALT=fulloutjoinSNPs[currentchrID][sampled_idx_find_satisfied][2]
                snprec_in_toplevel=dbvariantstools.operateDB("select","select * from "+toplevelsnptable+" where chrID='"+currentchrID+"' and snp_pos='"+str(fulloutjoinSNPs[currentchrID][sampled_idx_find_satisfied][0])+"'")                
                if (not snprec_in_toplevel or snprec_in_toplevel[0][outgroupidx_in_topleveltable[0]]==None):
                    dynamicIU_toptable_obj.insertorUpdatetopleveltable({currentchrID:[(snp_pos,REF,ALT)]}, flankseqfafile, 5)
                    continue#shoud excute only once for one position , to add a snp recs in toplevel table
                elif snprec_in_toplevel and snprec_in_toplevel[0][outgroupidx_in_topleveltable[0]]!=None:
                    #snp and outgroup base exist in topleveltable , extract ancestral allele
                    OUTgroupAlt=snprec_in_toplevel[0][outgroupidx_in_topleveltable[0]].upper().strip()
                    OUTgroup1=re.split(r",",snprec_in_toplevel[0][outgroupidx_in_topleveltable[0]+1].strip())
                    if len(OUTgroup1)==2 and OUTgroup1[0].strip()=='0' and int(OUTgroup1[1])>=minoutgroupdepth and ALT.upper()==OUTgroupAlt:
                        A_base_idx=1
                    elif len(OUTgroup1)==2 and OUTgroup1[1].strip()=='0' and int(OUTgroup1[0])>=minoutgroupdepth:
                        A_base_idx=0
                    else:
                        print("skip this snp,outgroup coverage not sufficent or different direction")
                        sampled_idx_find_satisfied+=direction
                        continue
                    #check minAN and MQ
                    NoOfAllele1Obsed=[];NoOfAllele2Obsed=[];NoAl1=0;NoAl2=0
                    pop_idx=-1
                    ################## one pop by one pop
                    for vcftable_name,minAN in vcffilelist:
                        AN=-1;AC=-1;MQvalue=-1;pop_idx+=1
#                         print(fulloutjoinSNPs[currentchrID][sampled_idx_find_satisfied][3+pop_idx],file=recf)
                        if fulloutjoinSNPs[currentchrID][sampled_idx_find_satisfied][3+pop_idx]==None:
                            AC=0;MQvalue=28;AF=0
                            AN=int(minAN)
#                             print(vcftable_name,"\nAN",AN,"AC",AC)
                        else:
                            INFO,FORMAT,INDVDlist=fulloutjoinSNPs[currentchrID][sampled_idx_find_satisfied][3+pop_idx]
                            try:#try collect MQ,AC AN
                                MQvalue=float(re.search(r"MQ=([\d\.]+);",INFO).group(1))
                                if re.search(r"indvd[^/]+",vcftable_name)!=None:
                                    AF = float(re.search(r"AF=([\d\.e-]+)[;,]", INFO).group(1))
                                    AC = int(re.search(r"AC=(\d+)[;,]", INFO).group(1))
                                    AN = int(re.search(r"AN=(\d+)[;,]", INFO).group(1))
#                                     print("indvd",vcftable_name,"\nAN",AN,"AC",AC)
                                elif re.search(r"pool[^/]+",vcftable_name)!=None:
                                    refdep = 0;altalleledep = 0
                                    AD_idx = (re.split(":", FORMAT)).index("AD")
                                    for indsample in INDVDlist:
                                        if len(re.split(":",indsample))==1:#./.
                                            continue
                                        AD_depth=re.split(",", re.split(":", indsample)[AD_idx])
                                        refdep += int(AD_depth[0])
                                        altalleledep += int(AD_depth[1])
                                    DP=altalleledep+refdep
                                    AF=altalleledep/DP
                                    AC=int(indvdNo_Of_POOL*2*AF)
#                                     print("pool",vcftable_name,AC)
                                    if DP/2>=int(minAN):
                                        AN=indvdNo_Of_POOL*2
                                    else:
#                                         print("DP/2<minAN",file=recf)
                                        break
                                    #for filter only,sometimes it will greater than indvdNo_Of_POOL
                            except:
#                                 print("exception",file=recf)
                                break
                        if MQvalue>=28 and  AN>=int(minAN) :
                            print("pop passed thersold",str(snp_pos),vcftable_name)
                            NoOfAllele1Obsed.append(AN-AC)
                            NoOfAllele2Obsed.append(AC)
                            if options.quantizing!=None and pop_idx in vcftablesidxlist_toquantizing:
                                NoAl1+=(1-AF)*(noofindvds2quantizing/len(options.quantizing))
                                NoAl2+=AF*(noofindvds2quantizing/len(options.quantizing))
                        else:
#                             print("break threshold",file=recf)
                            break
                    else:######## normal finished means position need to be print into dadiinputfile    ####################################
                        if [0 for e in  NoOfAllele1Obsed]==NoOfAllele1Obsed or [0 for e in NoOfAllele2Obsed]==NoOfAllele2Obsed:
                            print("skip pos that fixed in all pop")
                            sampled_idx_find_satisfied+=direction
                        else:
                            totalduiltsnp+=1
                            print("recode passed thersold",str(snp_pos))
                            continuesearch=1
                        continue
                    ############################### finish pop loop    #######################################
                
                else:
                    print("warning ! this may be a except situation")
                print("continue search")
                sampled_idx_find_satisfied+=direction
            else:#find
                print("find")
                contextwithinspeces=snprec_in_toplevel[0][5].upper()
                if A_base_idx==1:
                    contextoutgroup=contextwithinspeces[0]+ALT+contextwithinspeces[2]
                elif A_base_idx==0:
                    contextoutgroup=contextwithinspeces
                print(contextwithinspeces,contextoutgroup,REF,sep="\t",end="\t",file=dadisnpfile)
                print(*NoOfAllele1Obsed,sep="\t",end="\t",file=dadisnpfile)
                print(str(round(NoAl1)),ALT,sep="\t",end="\t",file=dadisnpfile)
                print(*NoOfAllele2Obsed,sep="\t",end="\t",file=dadisnpfile)
                print(str(round(NoAl2)),currentchrID.replace(".","_"),str(snp_pos),sep="\t",file=dadisnpfile)
    realdilutetime=totalduiltsnp/totalsnp
    realduiltelength=totallength*realdilutetime
    print("totallength",totallength)
    print("totalsnp",totalsnp,"totalduiltsnp(shoud equal to No. of SNP in dadi inputfile)",totalduiltsnp,"realduiltelength",realduiltelength,"totallengthduilt",totallengthduilt)
    dadisnpfile.close()
    flankseqfafile.close()
#     recf.close()
#                 snprec_in_toplevel=dbvariantstools.operateDB("select","select * from "+toplevelsnptable+" where chrID='"+currentchrID+"' and snp_pos='"+str(fulloutjoinSNPs[currentchrID][sampled_idx_find_satisfied][0])+"'")
                
                    
