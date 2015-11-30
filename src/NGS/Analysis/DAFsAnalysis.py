# -*- coding: UTF-8 -*-
'''
Created on 2015-8-27

@author: liurui
'''
import src.NGS.BasicUtil.DBManager as dbm
from optparse import OptionParser
import os,re,copy
from NGS.BasicUtil import Util, VCFutil

parser = OptionParser()
parser.add_option("-i", "--interval", dest="interval", nargs=3,
                  help="minvalue maxvalue breaks", metavar="FILE")
# parser.add_option("-v", "--vcffilename", dest="vcffilename",action="append",nargs=2,help="vcffilename multiple")
parser.add_option("-o","--outfileprename",dest="outfileprename",help="outfilepreName with path")
parser.add_option("-b","--bedlikefile",dest="bedlikefile",action="append",nargs=3,help="bedlikefile vcffilename minANforindvd_or_minDPforpool")#
parser.add_option("-C","--cdsfile",action="append",dest="cdsfiles",nargs=2,default=[],help="file depthfile")#
parser.add_option("-U","--utrfile",action="append",dest="utrfiles",nargs=2,default=[],help="file depthfile")#
parser.add_option("-I","--intronfile",action="append",dest="intronfiles",nargs=2,default=[],help="file depthfile")#
parser.add_option("-G","--interGenic",action="append",dest="interGenics",nargs=2,default=[],help="file depthfile")#
parser.add_option("-t","--topleveltablename",dest="topleveltablename",help="speciesName in table")                                                                                                                                                         
(options, args) = parser.parse_args()
pathtovcftools="/pub/tool/vcftools_0.1.12b/bin/vcftools "
dbvariantstools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.vcfdbname)
topleveltablename=options.topleveltablename
mindeptojudgefix_forancestral=30
#####################
VCFobj={};vcfnameKEY_depthfilename_titlenameVALUE_tojudgeancestrall={};vcfnameKEY_depthobjVALUE_tojudgeancestral={}
VCFobj["wigeon"]=VCFutil.VCF_Data("/home/bioinfo/liurui/data/vcffiles/uniqmap/taihudomesticgoose/taihudomesticgoose.pool.withindel.vcf")
VCFobj["fanya"]=VCFutil.VCF_Data("/home/bioinfo/liurui/data/vcffiles/uniqmap/fanya/fanya._pool.withindel.vcf")
vcfnameKEY_depthfilename_titlenameVALUE_tojudgeancestrall["wigeon"]=Util.GATK_depthfile("/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth","/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth.index")#here is a temp trick not a error
vcfnameKEY_depthfilename_titlenameVALUE_tojudgeancestrall["fanya"]=Util.GATK_depthfile("/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth","/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth.index")
vcfnameKEY_depthobjVALUE_tojudgeancestral["wigeon"]=["/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth",9]#here is a temp trick not a error
vcfnameKEY_depthobjVALUE_tojudgeancestral["fanya"]=["/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth",3]
####################################
minvalue = float(options.interval[0])
maxvalue = float(options.interval[1])
breaks = int(options.interval[2])
d_increase = (maxvalue - minvalue) / breaks
if __name__ == '__main__':
    ########## prepare ##########
    print("depthfile is not used currently")
    intervalFileName=options.outfileprename+".interval"
    intervalfile=open(intervalFileName,"w")
    while minvalue + d_increase<= maxvalue-d_increase :
        print(str(minvalue),str(minvalue+d_increase),sep="\t",file=intervalfile)
        minvalue+=d_increase
    else:
        if minvalue<maxvalue-d_increase:
            print("never get here",minvalue,maxvalue)
        print(str(minvalue),str(1.0),sep="\t",file=intervalfile)
        intervalfile.close()
    intervalfile=open(intervalFileName,'r')
    DAFintervalMap_SNPcounts_template={}
    MAFintervalMap_SNPcounts_template={}
    AFintervalMap_SNPcounts_template={}
    DAFsMapByFileName={}
    MAFsMapByFileName={}
    AFsMapByFileName={}
    for line in intervalfile:
        linelist=re.split(r'\s+',line.strip())
        DAFintervalMap_SNPcounts_template[float(linelist[0]),float(linelist[1])]=0
        MAFintervalMap_SNPcounts_template[float(linelist[0]),float(linelist[1])]=0
        AFintervalMap_SNPcounts_template[float(linelist[0]),float(linelist[1])]=0

    ################################bed##############
    print("start caculate bedfiles")
    for bedfilename,vcffilename,minAN in options.bedlikefile:#each col
        DAFintervalMap_SNPcounts=copy.deepcopy(DAFintervalMap_SNPcounts_template)
        MAFintervalMap_SNPcounts=copy.deepcopy(MAFintervalMap_SNPcounts_template)
        AFintervalMap_SNPcounts=copy.deepcopy(AFintervalMap_SNPcounts_template)
        bedfilenamewithoutpath=re.search(r"[^/]*$",bedfilename).group(0)
        vcffilenamewithoutpath=re.search(r"[^/]*$",vcffilename).group(0)
        os.system("rm "+vcffilename+bedfilenamewithoutpath+".recode.vcf"+" "+vcffilename+bedfilenamewithoutpath+".recode.vcf.myindex")

        if re.search(r"indvd[^/]+",vcffilename)!=None:
            a=os.system(pathtovcftools+" --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO AF --recode-INFO AN --bed "+bedfilename+" --vcf "+vcffilename+" --out "+vcffilename+bedfilenamewithoutpath)
        elif re.search(r"pool[^/]+",vcffilename)!=None:
            a=os.system(pathtovcftools+" --remove-indels --min-alleles 2 --max-alleles 2 --recode  --bed "+bedfilename+" --vcf "+vcffilename+" --out "+vcffilename+bedfilenamewithoutpath)
    #         filteredvcffilenamelist.append(vcffilename+bedfilenamewithoutpath+".recode.vcf")
        if a!=0:
            print("Polypolidy")
            print("""awk '{OFS="\\t";print $1,$2,$3,$4,$5,$6,$7,$8,"GT:AD:DP:GQ:PL","."}' """+vcffilename+" > "+vcffilename+"cut.vcf")
            b=os.system("""awk '{OFS="\\t";if($1!~/^#/){print $1,$2,$3,$4,$5,$6,$7,$8,"GT:AD:DP:GQ:PL","."}else if($0~/#CHROM/){print "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tpopname"}else{print $0}}' """+vcffilename+" > "+vcffilename+"cut.vcf")
            if b!=0:
                print("error b")
                exit(-1)
            treatPolypolidyasindvd=True
            a=os.system(pathtovcftools+" --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO AF --recode-INFO AN --bed "+bedfilename+" --vcf "+vcffilename+"cut.vcf --out "+vcffilename+bedfilenamewithoutpath)
        if a!=0:
            print(pathtovcftools+" --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO AF --recode-INFO AN --bed "+bedfilename+" --vcf "+vcffilename+"cut.vcf --out "+vcffilename+bedfilenamewithoutpath)
            print("error pathtovcftools")
            exit(-1)
            os.system("rm "+vcffilename+"cut.vcf")

        chromset=[]
#         VCFlistmapBycurchr=[]
#         for filterdvcf in filteredvcffilenamelist:#produce myindex file,vcfobj
#             vcfobjlist.append(VCFutil.VCF_Data(filterdvcf))
#             chromset+=vcfobjlist[-1].chromOrder
#             VCFlistmapBycurchr.append({})
        vcfobj=VCFutil.VCF_Data(vcffilename+bedfilenamewithoutpath+".recode.vcf")
        chromset=vcfobj.chromOrder
#         chromset=list(set(chromset))
        for currentchrID in chromset:
            allsnprecinAchr_mapbyvcfname={}
#             for vcfobj_idx in range(len(vcfobjlist)):
#                 VCFlistmapBycurchr[vcfobj_idx]={}
#                 VCFlistmapBycurchr[vcfobj_idx][currentchrID]=vcfobjlist[vcfobj_idx].getVcfListByChrom(currentchrID)
# #                 VCFlistmapBycurchr.append({currentchrID:vcfobj})
#             print(VCFlistmapBycurchr,file=open("aaaaa.txt",'a'))
#             alignedSNP_onechr=Util.alinmultPopSnpPos(VCFlistmapBycurchr, "o")
            alignedSNP_onechr=vcfobj.getVcfListByChrom(currentchrID)
#             countedAF=0
            for snp_aligned in alignedSNP_onechr:
                curpos=int(snp_aligned[0])
                ######caculate AF###############
                if re.search(r"indvd[^/]+",vcffilename+bedfilenamewithoutpath+".recode.vcf")!=None or treatPolypolidyasindvd:
                    AF = float(re.search(r"AF=([\d\.]+)", snp_aligned[3]).group(1))
                    AN = float(re.search(r"AN=([\d]+)", snp_aligned[3]).group(1))
                    if AN<int(minAN):
                        continue
                elif (not treatPolypolidyasindvd) and re.search(r"pool[^/]+",vcffilename+bedfilenamewithoutpath+".recode.vcf")!=None:
                    refdep = 0;altalleledep = 0
                    DP_idx=(re.split(":", snp_aligned[4])).index("DP")
                    AD_idx = (re.split(":", snp_aligned[4])).index("AD")  # gatk GT:AD:DP:GQ:PL   
                    for sample in snp_aligned[5]:
                        if len(re.split(":", sample)) == 1:  # ./.
                            continue                        
                        AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                        try :
                            refdep += int(AD_depth[0])
                            altalleledep += int(AD_depth[1])
                        except ValueError:
                            print(sample, end="|")
                    if (refdep==altalleledep and altalleledep==0) or altalleledep+ refdep<int(minAN):
                        print("no sample available in this pop")
                        continue
                    AF=altalleledep/(altalleledep+refdep)
# +++++++++++++++++++++++++++++++++++++++++++++==

                if AF>=0.5:
                    MAF=AF
                else:
                    MAF=1-AF
                for a,b in MAFintervalMap_SNPcounts.keys():
                    if MAF>=a and MAF<b and MAF!=0:
                        MAFintervalMap_SNPcounts[a,b]+=1
                        break
                snp=dbvariantstools.operateDB("select","select * from "+topleveltablename+" where chrID='"+currentchrID+"' and snp_pos='"+str(curpos)+"'")
                if not snp:
                    if allsnprecinAchr_mapbyvcfname=={}:
#                         allsnprecinAchr_mapbyvcfname["wigeon"]=VCFobj["wigeon"].getVcfListByChrom(currentchrID)
                        allsnprecinAchr_mapbyvcfname["fanya"]=VCFobj["fanya"].getVcfListByChrom(currentchrID)
                    low0=0
                    high0=len(allsnprecinAchr_mapbyvcfname["fanya"])-1
                    while low0<=high0:
                        mid0=(low0+high0)>>1
                        if allsnprecinAchr_mapbyvcfname["fanya"][mid0][0]<curpos:
                            low0=mid0+1
                        elif allsnprecinAchr_mapbyvcfname["fanya"][mid0][0]>curpos:
                            high0=mid0-1
                        else:#find the pos
                            pos, REF, archicpop_ALT, INFO,FORMAT,samples = allsnprecinAchr_mapbyvcfname["fanya"][mid0]
                            refdep=0;altalleledep=0
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
                            if refdep!=0  or altalleledep<mindeptojudgefix_forancestral:
                                A_base_idx=-1#alt_allele is the ancestral allele
                                low0=high0+1
                            elif refdep==0 and altalleledep>=mindeptojudgefix_forancestral:
                                A_base_idx=1
                                break
                            
                    else:
                        depth_linelist=vcfnameKEY_depthfilename_titlenameVALUE_tojudgeancestrall["fanya"].getdepthByPos_optimized(currentchrID,curpos)
                        if int(depth_linelist[vcfnameKEY_depthobjVALUE_tojudgeancestral["fanya"][1]])>=mindeptojudgefix_forancestral:
                            A_base_idx=0
                        else:
                            A_base_idx=-1

                    if A_base_idx==-1:
                        print(currentchrID,curpos,"snp not find in db and not sufficent in fanya,skip")
                        continue
                else:
                    wigeondepthlist1=re.split(r",",snp[0][13])
                    fanyadepthlist=re.split(r",",snp[0][9])
                    if len(wigeondepthlist1)==2 and len(fanyadepthlist)==2 and (int(wigeondepthlist1[0]) + int(wigeondepthlist1[1])>=mindeptojudgefix_forancestral or int(fanyadepthlist[0]) + int(fanyadepthlist[1])>=mindeptojudgefix_forancestral) and ((wigeondepthlist1[0].strip()=="0" and fanyadepthlist[0].strip()=="0") or (wigeondepthlist1[1].strip()=="0" and fanyadepthlist[1].strip()=="0") ):
                        if wigeondepthlist1[0].strip()=="0" and fanyadepthlist[0].strip()=="0":
                            A_base_idx=1
                        elif wigeondepthlist1[1].strip()=="0" and fanyadepthlist[1].strip()=="0":
                            A_base_idx=0
                        else:
                            print("never get here")
                    elif (len(fanyadepthlist)==2  and int(fanyadepthlist[0]) + int(fanyadepthlist[1])>=mindeptojudgefix_forancestral and (fanyadepthlist[0].strip()=="0" or fanyadepthlist[1].strip()=="0" )):#   or (snp[0][7]=="no covered" and len(depthlist2)==2 and int(depthlist2[0]) + int(depthlist2[1])>=mindeptojudgefix and (depthlist2[1].strip()=="0" or depthlist2[0].strip()=="0")):
                        if  fanyadepthlist[0].strip()=="0":
                            A_base_idx=1
                        elif fanyadepthlist[1].strip()=="0":
                            A_base_idx=0
                    else:
                        print("skip snp",snp)
                        A_base_idx=-1
                        continue                    
                                     



                if A_base_idx==-1:
                    print(currentchrID,curpos,"snp not find in db and not sufficent in fanya,skip")
                    continue

                #######ancestrall allele judged ,average DAF ################

                if A_base_idx==0:
                    DAF=1-AF
                else:
                    DAF=AF
                ##### count ########
                for a,b in AFintervalMap_SNPcounts.keys():
                    if AF>=a and AF<b and AF!=0:#
                        AFintervalMap_SNPcounts[a,b]+=1
                        break
                for a,b in DAFintervalMap_SNPcounts.keys():
                    if DAF>=a and DAF<b and DAF!=0:
                        DAFintervalMap_SNPcounts[a,b]+=1
                        break
        AFsMapByFileName[bedfilenamewithoutpath+vcffilenamewithoutpath]=copy.deepcopy(AFintervalMap_SNPcounts)
        DAFsMapByFileName[bedfilenamewithoutpath+vcffilenamewithoutpath]=copy.deepcopy(DAFintervalMap_SNPcounts)
        MAFsMapByFileName[bedfilenamewithoutpath+vcffilenamewithoutpath]=copy.deepcopy(MAFintervalMap_SNPcounts)
    

    os.system("rm "+vcffilename+bedfilenamewithoutpath+".recode.vcf"+" "+vcffilename+bedfilenamewithoutpath+".recode.vcf.myindex")
    ######start caculate #########################
    catalogfilslistlist=[]
    if options.intronfiles!=[]:
        catalogfilslistlist.append(options.intronfiles+["intron"])
    if options.utrfiles!=[]:
        catalogfilslistlist.append(options.utrfiles+["utr"])
    if options.interGenics!=[]:
        catalogfilslistlist.append(options.interGenics+["intergenic"])
    for catalogfilslist in catalogfilslistlist:
        print("catalogfilslistlist",catalogfilslist)
        tag=catalogfilslist[-1]
        cdsfilenameslist=[];cds_depthfileconfig={};AF_idxlist_cds=[]
        AFintervalMap_SNPcounts=copy.deepcopy(AFintervalMap_SNPcounts_template)
        DAFintervalMap_SNPcounts=copy.deepcopy(DAFintervalMap_SNPcounts_template)    
        for cdsfile,depthfilename in catalogfilslist[:-1]:
            cdsfilenameslist.append(cdsfile)
            if depthfilename.lower()!="none":
                cds_depthfileconfig[cdsfile]=[]
                fp=open(depthfilename,'r')
                for line in fp:
                    depthfile_obj=re.search(r"depthfilename=(.*)",line.strip())
                    if depthfile_obj!=None:
                        cds_depthfileconfig[cdsfile].append(depthfile_obj.group(1).strip())
                    elif line.split():
                        cds_depthfileconfig[cdsfile].append(line.strip())
                fp.close()
            else:
                pass  
        #######prepare chromlist #######
        for cdsfname in cdsfilenameslist:
            t=open(cdsfname,'r')
            AF_idxlist_cds.append(re.split(r"\t",t.readline().strip()).index("AF"))
            t.close()
            os.system("awk 'NR>1{print $1}' "+cdsfname+"|sort|uniq|sort >"+cdsfname+"_chrom")
        for cdsf_idx in range(1,len(cdsfilenameslist)):
            os.system("cat "+cdsfilenameslist[cdsf_idx-1]+"_chrom "+cdsfilenameslist[cdsf_idx]+"_chrom|sort|uniq|sort > temp_chrom")
            os.system("rm "+cdsfilenameslist[cdsf_idx-1]+"_chrom ")
            os.system("mv temp_chrom "+cdsfilenameslist[cdsf_idx]+"_chrom")
        os.system("mv "+cdsfilenameslist[-1]+"_chrom cdschromlist")
        t=open("cdschromlist","r")
        chromlist=t.readlines();t.close()
        os.system("rm cdschromlist")
        ######chrom list readed ##########
        for chrom in chromlist:
            curchrom=chrom.strip()
            allsnprecinAchr_mapbyvcfname={}
            wild_CurRecsLinelist=[];wild_CurPosRecs=[];posOfCurRecwild=[]
            wildcdsfilelist=[];wildcdsfilenamelist=[]#records filename to remove 
            idx=len(cdsfilenameslist)
            depthobjmap={};species_idx_map={};randomstr=Util.random_str()
            for cdsfname in cdsfilenameslist[::-1]:
                depthobjmap[cdsfname]=Util.GATK_depthfile(cds_depthfileconfig[cdsfname][0],cds_depthfileconfig[cdsfname][0]+".index")
                species_idx_map[cdsfname]=[]
                for titlename in cds_depthfileconfig[cdsfname][1:]:
                    species_idx_map[cdsfname].append(depthobjmap[cdsfname].title.index("Depth_for_"+titlename))
                idx-=1
                filename=cdsfname+"_one_chrom"+randomstr
                os.system("awk '$1~/"+curchrom+"/{print $0}' "+cdsfname+">"+filename)
                wildcdsfilelist.append(open(filename,'r'))
                wildcdsfilenamelist.append(filename)
            wildcdsfilelist.reverse()
            ########## collect delta_AF #################################################################
            for wf_idx in range(len(wildcdsfilelist)):
                line=wildcdsfilelist[wf_idx].readline()
                if line.split():
                    wild_CurRecsLinelist.append(re.split(r"\s+",line.strip()))
                    posOfCurRecwild.append(int(wild_CurRecsLinelist[wf_idx][1]))###############
                else:
                    wild_CurRecsLinelist.append("endline")
                    posOfCurRecwild.append(999999999999999999999999999)
            while wild_CurRecsLinelist!=["endline"]*len(wildcdsfilelist):
                wild_CurPosRecs=[]
                curpos=min(posOfCurRecwild)
                for wf_idx in range(len(wildcdsfilelist)):
                    if  wild_CurRecsLinelist[wf_idx]=="endline" or int(wild_CurRecsLinelist[wf_idx][1])>curpos:#None means reach the end of the line
                        if posOfCurRecwild[wf_idx]<=curpos:#endline
                            posOfCurRecwild[wf_idx]=999999999999999999999999999999
                        depth_linelist=depthobjmap[cdsfilenameslist[wf_idx]].getdepthByPos_optimized(curchrom,curpos)
                        sum_depth=0
                        for idx in species_idx_map[cdsfilenameslist[wf_idx]]:
                            print("should never reach here")
                            sum_depth+=int(depth_linelist[idx])
                        wild_CurPosRecs.append(["unknow"]*(AF_idxlist_cds[wf_idx]+1))
                        if sum_depth>mindeptojudgefix_forancestral:
                            wild_CurPosRecs[wf_idx][AF_idxlist_cds[wf_idx]]=0
                    elif int(wild_CurRecsLinelist[wf_idx][1])==curpos:
                        wild_CurPosRecs.append(copy.deepcopy(wild_CurRecsLinelist[wf_idx]))
                        line=wildcdsfilelist[wf_idx].readline()
                        if line.split():
                            wild_CurRecsLinelist[wf_idx]=re.split(r"\s+",line.strip())
                            posOfCurRecwild[wf_idx]=int(wild_CurRecsLinelist[wf_idx][1])
                        else:
                            wild_CurRecsLinelist[wf_idx]="endline"
                if wild_CurPosRecs==[]:
                    continue
                w_af=0;aaf=0
                #determin derived allele
                snp=dbvariantstools.operateDB("select","select * from "+topleveltablename+" where chrID='"+curchrom+"' and snp_pos='"+str(curpos)+"'")
                if not snp:
                    if allsnprecinAchr_mapbyvcfname=={}:
#                         allsnprecinAchr_mapbyvcfname["wigeon"]=VCFobj["wigeon"].getVcfListByChrom(curchrom)
                        allsnprecinAchr_mapbyvcfname["fanya"]=VCFobj["fanya"].getVcfListByChrom(curchrom)
                    low0=0
                    high0=len(allsnprecinAchr_mapbyvcfname["fanya"])-1
                    while low0<=high0:
                        mid0=(low0+high0)>>1
                        if allsnprecinAchr_mapbyvcfname["fanya"][mid0][0]<curpos:
                            low0=mid0+1
                        elif allsnprecinAchr_mapbyvcfname["fanya"][mid0][0]>curpos:
                            high0=mid0-1
                        else:#find the pos
                            pos, REF, archicpop_ALT, INFO,FORMAT,samples = allsnprecinAchr_mapbyvcfname["fanya"][mid0]
                            refdep=0;altalleledep=0
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
                            if refdep!=0  or altalleledep<mindeptojudgefix_forancestral:
                                low0=high0+1
                                A_base_idx=-1#alt_allele is the ancestral allele
                            elif refdep==0 and altalleledep>=mindeptojudgefix_forancestral:
                                A_base_idx=1
                                break
                            
                    else:
                        depth_linelist=vcfnameKEY_depthfilename_titlenameVALUE_tojudgeancestrall["fanya"].getdepthByPos_optimized(curchrom,curpos)
                        if int(depth_linelist[vcfnameKEY_depthobjVALUE_tojudgeancestral["fanya"][1]])>=mindeptojudgefix_forancestral:
                            A_base_idx=0
                        else:
                            A_base_idx=-1
    #                 low1=0
    #                 high1=len(allsnprecinAchr_mapbyvcfname["wigeon"])-1
    #                 while low1<=high1:
    #                     mid1=(low1+high1)>>1
    #                     if allsnprecinAchr_mapbyvcfname["wigeon"][mid1][0]<curpos:
    #                         low1=mid1+1
    #                     elif allsnprecinAchr_mapbyvcfname["wigeon"][mid1][0]>curpos:
    #                         high1=mid1-1
    #                     else:
    #                         pos, REF, archicpop_ALT, INFO,FORMAT,samples = allsnprecinAchr_mapbyvcfname["wigeon"][mid1]
    #                         refdep=0;altalleledep=0
    #                         AD_idx=(re.split(":",FORMAT)).index("AD")#gatk GT:AD:DP:GQ:PL
    #                         for sample in samples:
    #                             if len(re.split(":",sample))==1:# ./.
    #                                 continue
    #                             AD_depth=re.split(",",re.split(":",sample)[AD_idx])
    #                             try :
    #                                 refdep+=int(AD_depth[0])
    #                                 altalleledep+=int(AD_depth[1])
    #                             except ValueError:
    #                                 print("Ancestralallele.fillAncestral except ValueError",sample,end="")
    #                         if  refdep!=0:
    #                             pass
    # #                             A_base_idx=-1
    #                         elif refdep==0 and altalleledep>=mindeptojudgefix and A_base_idx==1:
    #                             A_base_idx=1
    #                         break
    #                 else:
    #                     if A_base_idx==1:
    #                         A_base_idx=-1
    #                     else:
    #                         A_base_idx=0
                    if A_base_idx==-1:
                        print(curchrom,curpos,"snp not find in db and not sufficent in fanya,skip")
                        continue
##########
                else:
                    wigeondepthlist1=re.split(r",",snp[0][13])
                    fanyadepthlist=re.split(r",",snp[0][9])
                    if len(wigeondepthlist1)==2 and len(fanyadepthlist)==2 and (int(wigeondepthlist1[0]) + int(wigeondepthlist1[1])>=mindeptojudgefix_forancestral or int(fanyadepthlist[0]) + int(fanyadepthlist[1])>=mindeptojudgefix_forancestral) and ((wigeondepthlist1[0].strip()=="0" and fanyadepthlist[0].strip()=="0") or (wigeondepthlist1[1].strip()=="0" and fanyadepthlist[1].strip()=="0") ):
                        if wigeondepthlist1[0].strip()=="0" and fanyadepthlist[0].strip()=="0":
                            A_base_idx=1
                        elif wigeondepthlist1[1].strip()=="0" and fanyadepthlist[1].strip()=="0":
                            A_base_idx=0
                        else:
                            print("never get here ,how could it be",wigeondepthlist1,fanyadepthlist)
                    elif (len(fanyadepthlist)==2  and int(fanyadepthlist[0]) + int(fanyadepthlist[1])>=mindeptojudgefix_forancestral+5 and (fanyadepthlist[0].strip()=="0" or fanyadepthlist[1].strip()=="0" )):#   or (snp[0][7]=="no covered" and len(depthlist2)==2 and int(depthlist2[0]) + int(depthlist2[1])>=mindeptojudgefix and (depthlist2[1].strip()=="0" or depthlist2[0].strip()=="0")):
                        if  fanyadepthlist[0].strip()=="0":
                            A_base_idx=1
                        elif fanyadepthlist[1].strip()=="0":
                            A_base_idx=0
                    else:
                        print("skip snp",snp)
                        continue
                #######ancestrall allele judged ,average DAF ################
                idx=0;w_unknowcount=0
                for e in wild_CurPosRecs:
                    if e[AF_idxlist_cds[idx]]=="unknow":
                        idx+=1;w_unknowcount+=1
                        continue
                    aaf+=float(e[AF_idxlist_cds[idx]])
                    if A_base_idx==1:
                        w_af+=(1-float(e[AF_idxlist_cds[idx]]))
                    else:
                        w_af+=float(e[AF_idxlist_cds[idx]])
                    idx+=1
                if len(wild_CurPosRecs)-w_unknowcount==0:
                    continue
                DAF=(w_af/(len(wild_CurPosRecs)-w_unknowcount))
                AAF=(aaf/(len(wild_CurPosRecs)-w_unknowcount))
                for a,b in AFintervalMap_SNPcounts.keys():
                    if AAF>=a and AAF<b and AAF!=0:
                        AFintervalMap_SNPcounts[a,b]+=1
                        break
                for a,b in DAFintervalMap_SNPcounts.keys():
                    if DAF>=a and DAF<b and DAF!=0:
                        DAFintervalMap_SNPcounts[a,b]+=1
                        break
            for f in wildcdsfilelist:
                f.close()
            for filename in wildcdsfilenamelist:
                os.system("rm "+filename)
        AFsMapByFileName[tag]=copy.deepcopy(AFintervalMap_SNPcounts)
        DAFsMapByFileName[tag]=copy.deepcopy(DAFintervalMap_SNPcounts)
    outfile=open(options.outfileprename+"gtfAAF","w")
    print("AFbin",end="\t",file=outfile)
    for tag in sorted(AFsMapByFileName.keys()):
        print(tag,end="\t",file=outfile)
    print("",file=outfile)
    for a,b in AFintervalMap_SNPcounts_template.keys():
        print(str(a),str(b),end="\t",file=outfile)
        for tag in sorted(AFsMapByFileName.keys()):
            print(str(AFsMapByFileName[tag][a,b]),end="\t",file=outfile)
        print("",file=outfile)
    outfile.close()
    outfile=open(options.outfileprename+"gtfDAF","w")
    print("DAFbin",end="\t",file=outfile)
    for tag in sorted(DAFsMapByFileName.keys()):
        print(tag,end="\t",file=outfile)
    print("",file=outfile)
    for a,b in sorted(DAFintervalMap_SNPcounts_template.keys()):
        print(str(a),str(b),end="\t",file=outfile)
        for tag in sorted(DAFsMapByFileName.keys()):
            print(str(DAFsMapByFileName[tag][a,b]),end="\t",file=outfile)
        print("",file=outfile)
    outfile.close()
    #################cds######################################
    if options.cdsfiles!=[]:
        cdsfilenameslist=[];cds_depthfileconfig={};AF_idxlist_cds=[]
        DAFintervalMap_SNPcounts_sys=copy.deepcopy(DAFintervalMap_SNPcounts_template)
        DAFintervalMap_SNPcounts_nonsys=copy.deepcopy(DAFintervalMap_SNPcounts_template)
        DAFintervalMap_SNPcounts_nonsense=copy.deepcopy(DAFintervalMap_SNPcounts_template)
        AFintervalMap_SNPcounts_sys=copy.deepcopy(AFintervalMap_SNPcounts_template)
        AFintervalMap_SNPcounts_nonsys=copy.deepcopy(AFintervalMap_SNPcounts_template)
        AFintervalMap_SNPcounts_nonsense=copy.deepcopy(AFintervalMap_SNPcounts_template)
        for cdsfile,depthfilename in options.cdsfiles:
            cdsfilenameslist.append(cdsfile)
            if depthfilename.lower()!="none":
                cds_depthfileconfig[cdsfile]=[]
                fp=open(depthfilename,'r')
                for line in fp:
                    depthfile_obj=re.search(r"depthfilename=(.*)",line.strip())
                    if depthfile_obj!=None:
                        cds_depthfileconfig[cdsfile].append(depthfile_obj.group(1).strip())
                    elif line.split():
                        cds_depthfileconfig[cdsfile].append(line.strip())
                fp.close()
            else:
                pass  
        #######prepare chromlist #######
        for cdsfname in cdsfilenameslist:
            t=open(cdsfname,'r')
            AF_idxlist_cds.append(re.split(r"\t",t.readline().strip()).index("AF"))
            t.close()
            os.system("awk 'NR>1{print $1}' "+cdsfname+"|sort|uniq|sort >"+cdsfname+"_chrom")
        for cdsf_idx in range(1,len(cdsfilenameslist)):
            os.system("cat "+cdsfilenameslist[cdsf_idx-1]+"_chrom "+cdsfilenameslist[cdsf_idx]+"_chrom|sort|uniq|sort > temp_chrom")
            os.system("rm "+cdsfilenameslist[cdsf_idx-1]+"_chrom ")
            os.system("mv temp_chrom "+cdsfilenameslist[cdsf_idx]+"_chrom")
        os.system("mv "+cdsfilenameslist[-1]+"_chrom cdschromlist")
        t=open("cdschromlist","r")
        chromlist=t.readlines();t.close()
        os.system("rm cdschromlist")
        ######chrom list readed ##########
        for chrom in chromlist:
            allsnprecinAchr_mapbyvcfname={}
            curchrom=chrom.strip()
            wild_CurRecsLinelist=[];wild_CurPosRecs=[];posOfCurRecwild=[]
            wildcdsfilelist=[];wildcdsfilenamelist=[]#records filename to remove 
            idx=len(cdsfilenameslist)
            depthobjmap={};species_idx_map={};randomstr=Util.random_str()
            for cdsfname in cdsfilenameslist[::-1]:
                depthobjmap[cdsfname]=Util.GATK_depthfile(cds_depthfileconfig[cdsfname][0],cds_depthfileconfig[cdsfname][0]+".index")
                species_idx_map[cdsfname]=[]
                for titlename in cds_depthfileconfig[cdsfname][1:]:
                    species_idx_map[cdsfname].append(depthobjmap[cdsfname].title.index("Depth_for_"+titlename))
                idx-=1
                filename=cdsfname+"_one_chrom"+randomstr
                os.system("awk '$1~/"+curchrom+"/{print $0}' "+cdsfname+">"+filename)
                wildcdsfilelist.append(open(filename,'r'))
                wildcdsfilenamelist.append(filename)
            wildcdsfilelist.reverse()
            ########## collect delta_AF #################################################################
            for wf_idx in range(len(wildcdsfilelist)):
                line=wildcdsfilelist[wf_idx].readline()
                if line.split():
                    wild_CurRecsLinelist.append(re.split(r"\s+",line.strip()))
                    posOfCurRecwild.append(int(wild_CurRecsLinelist[wf_idx][1]))###############
                else:
                    wild_CurRecsLinelist.append("endline")
                    posOfCurRecwild.append(999999999999999999999999999)
            while wild_CurRecsLinelist!=["endline"]*len(wildcdsfilelist):
                wild_CurPosRecs=[]
                curpos=min(posOfCurRecwild)
                for wf_idx in range(len(wildcdsfilelist)):
                    if  wild_CurRecsLinelist[wf_idx]=="endline" or int(wild_CurRecsLinelist[wf_idx][1])>curpos:#None means reach the end of the line
                        if posOfCurRecwild[wf_idx]<=curpos:#endline
                            posOfCurRecwild[wf_idx]=999999999999999999999999999999
                        depth_linelist=depthobjmap[cdsfilenameslist[wf_idx]].getdepthByPos_optimized(curchrom,curpos)
                        sum_depth=0
                        for idx in species_idx_map[cdsfilenameslist[wf_idx]]:
                            sum_depth+=int(depth_linelist[idx])
                        wild_CurPosRecs.append(["unknow"]*(AF_idxlist_cds[wf_idx]+1))
                        if sum_depth>mindeptojudgefix_forancestral:
                            wild_CurPosRecs[wf_idx][AF_idxlist_cds[wf_idx]]=0
                    elif int(wild_CurRecsLinelist[wf_idx][1])==curpos:
                        wild_CurPosRecs.append(copy.deepcopy(wild_CurRecsLinelist[wf_idx]))
                        line=wildcdsfilelist[wf_idx].readline()
                        if line.split():
                            wild_CurRecsLinelist[wf_idx]=re.split(r"\s+",line.strip())
                            posOfCurRecwild[wf_idx]=int(wild_CurRecsLinelist[wf_idx][1])
                        else:
                            wild_CurRecsLinelist[wf_idx]="endline"
                if wild_CurPosRecs==[]:
                    continue
                w_af=0;aaf=0
                #determin derived allele
                snp=dbvariantstools.operateDB("select","select * from "+topleveltablename+" where chrID='"+curchrom+"' and snp_pos='"+str(curpos)+"'")
                if not snp:
                    if allsnprecinAchr_mapbyvcfname=={}:
#                         allsnprecinAchr_mapbyvcfname["wigeon"]=VCFobj["wigeon"].getVcfListByChrom(curchrom)
                        allsnprecinAchr_mapbyvcfname["fanya"]=VCFobj["fanya"].getVcfListByChrom(curchrom)
                    low0=0
                    high0=len(allsnprecinAchr_mapbyvcfname["fanya"])-1
                    while low0<=high0:
                        mid0=(low0+high0)>>1
                        if allsnprecinAchr_mapbyvcfname["fanya"][mid0][0]<curpos:
                            low0=mid0+1
                        elif allsnprecinAchr_mapbyvcfname["fanya"][mid0][0]>curpos:
                            high0=mid0-1
                        else:#find the pos
                            pos, REF, archicpop_ALT, INFO,FORMAT,samples = allsnprecinAchr_mapbyvcfname["fanya"][mid0]
                            refdep=0;altalleledep=0
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
                            if refdep!=0  or altalleledep<mindeptojudgefix_forancestral:
                                A_base_idx=-1#alt_allele is the ancestral allele
                                low0=high0+1
                            elif refdep==0 and altalleledep>=mindeptojudgefix_forancestral:
                                A_base_idx=1
                                break
                            
                    else:
                        depth_linelist=vcfnameKEY_depthfilename_titlenameVALUE_tojudgeancestrall["fanya"].getdepthByPos_optimized(curchrom,curpos)
                        if int(depth_linelist[vcfnameKEY_depthobjVALUE_tojudgeancestral["fanya"][1]])>=mindeptojudgefix_forancestral:
                            A_base_idx=0
                        else:
                            A_base_idx=-1
    #                 low1=0
    #                 high1=len(allsnprecinAchr_mapbyvcfname["wigeon"])-1
    #                 while low1<=high1:
    #                     mid1=(low1+high1)>>1
    #                     if allsnprecinAchr_mapbyvcfname["wigeon"][mid1][0]<curpos:
    #                         low1=mid1+1
    #                     elif allsnprecinAchr_mapbyvcfname["wigeon"][mid1][0]>curpos:
    #                         high1=mid1-1
    #                     else:
    #                         pos, REF, archicpop_ALT, INFO,FORMAT,samples = allsnprecinAchr_mapbyvcfname["wigeon"][mid1]
    #                         refdep=0;altalleledep=0
    #                         AD_idx=(re.split(":",FORMAT)).index("AD")#gatk GT:AD:DP:GQ:PL
    #                         for sample in samples:
    #                             if len(re.split(":",sample))==1:# ./.
    #                                 continue
    #                             AD_depth=re.split(",",re.split(":",sample)[AD_idx])
    #                             try :
    #                                 refdep+=int(AD_depth[0])
    #                                 altalleledep+=int(AD_depth[1])
    #                             except ValueError:
    #                                 print("Ancestralallele.fillAncestral except ValueError",sample,end="")
    #                         if  refdep!=0:
    #                             pass
    # #                             A_base_idx=-1
    #                         elif refdep==0 and altalleledep>=mindeptojudgefix and A_base_idx==1:
    #                             A_base_idx=1
    #                         break
    #                 else:
    #                     if A_base_idx==1:
    #                         A_base_idx=-1
    #                     else:
    #                         A_base_idx=0
                    if A_base_idx==-1:
                        print(curchrom,curpos,"snp not find in db and not sufficent in fanya,skip")
                        continue
                else:
                    wigeondepthlist1=re.split(r",",snp[0][13])
                    fanyadepthlist=re.split(r",",snp[0][9])
                    if len(wigeondepthlist1)==2 and len(fanyadepthlist)==2 and (int(wigeondepthlist1[0]) + int(wigeondepthlist1[1])>=mindeptojudgefix_forancestral or int(fanyadepthlist[0]) + int(fanyadepthlist[1])>=mindeptojudgefix_forancestral) and ((wigeondepthlist1[0].strip()=="0" and fanyadepthlist[0].strip()=="0") or (wigeondepthlist1[1].strip()=="0" and fanyadepthlist[1].strip()=="0") ):
                        if wigeondepthlist1[0].strip()=="0" and fanyadepthlist[0].strip()=="0":
                            A_base_idx=1
                        elif wigeondepthlist1[1].strip()=="0" and fanyadepthlist[1].strip()=="0":
                            A_base_idx=0
                        else:
                            print("never get here")
                    elif (len(fanyadepthlist)==2  and int(fanyadepthlist[0]) + int(fanyadepthlist[1])>=mindeptojudgefix_forancestral and (fanyadepthlist[0].strip()=="0" or fanyadepthlist[1].strip()=="0" )):#   or (snp[0][7]=="no covered" and len(depthlist2)==2 and int(depthlist2[0]) + int(depthlist2[1])>=mindeptojudgefix and (depthlist2[1].strip()=="0" or depthlist2[0].strip()=="0")):
                        if  fanyadepthlist[0].strip()=="0":
                            A_base_idx=1
                        elif fanyadepthlist[1].strip()=="0":
                            A_base_idx=0
                    else:
                        print("skip snp",snp)
                        continue
                #######ancestrall allele judged ,average DAF ################
                idx=0;w_unknowcount=0
                for e in wild_CurPosRecs:
                    if e[AF_idxlist_cds[idx]]=="unknow":
                        idx+=1;w_unknowcount+=1
                        continue
                    aaf+=float(e[AF_idxlist_cds[idx]])
                    if A_base_idx==1:
                        w_af+=(1-float(e[AF_idxlist_cds[idx]]))
                    else:
                        w_af+=float(e[AF_idxlist_cds[idx]])
                    idx+=1
                if len(wild_CurPosRecs)-w_unknowcount==0:
                    continue
                AAF=(aaf/(len(wild_CurPosRecs)-w_unknowcount))
                DAF=(w_af/(len(wild_CurPosRecs)-w_unknowcount))
                for a,b in AFintervalMap_SNPcounts_nonsys.keys():
                    if AAF>=a and AAF<b and AAF!=0 and len(wild_CurPosRecs[0])>=14:
                        if wild_CurPosRecs[0][-3]==wild_CurPosRecs[0][-1]:
                            AFintervalMap_SNPcounts_sys[a,b]+=1
                        elif wild_CurPosRecs[0][-3].find("*")!=-1 or wild_CurPosRecs[0][-1].find("*")!=-1:
                            AFintervalMap_SNPcounts_nonsense[a,b]+=1
                        elif wild_CurPosRecs[0][-3]!=wild_CurPosRecs[0][-1]:
                            AFintervalMap_SNPcounts_nonsys[a,b]+=1
                        break
                for a,b in DAFintervalMap_SNPcounts_nonsys.keys():
                    if DAF>=a and DAF<b and DAF!=0 and len(wild_CurPosRecs[0])>=14:
                        if wild_CurPosRecs[0][-3]==wild_CurPosRecs[0][-1]:
                            DAFintervalMap_SNPcounts_sys[a,b]+=1
                        elif wild_CurPosRecs[0][-3].find("*")!=-1 or wild_CurPosRecs[0][-1].find("*")!=-1:
                            DAFintervalMap_SNPcounts_nonsense[a,b]+=1                            
                        elif wild_CurPosRecs[0][-3]!=wild_CurPosRecs[0][-1]:
                            DAFintervalMap_SNPcounts_nonsys[a,b]+=1

                        break
            for f in wildcdsfilelist:
                f.close()
            for filename in wildcdsfilenamelist:
                os.system("rm "+filename)
        DAFsMapByFileName["sys"]=copy.deepcopy(DAFintervalMap_SNPcounts_sys)
        DAFsMapByFileName["nonsys"]=copy.deepcopy(DAFintervalMap_SNPcounts_nonsys)
        DAFsMapByFileName["nonsense"]=copy.deepcopy(DAFintervalMap_SNPcounts_nonsense)
        AFsMapByFileName["sys"]=copy.deepcopy(AFintervalMap_SNPcounts_sys)
        AFsMapByFileName["nonsys"]=copy.deepcopy(AFintervalMap_SNPcounts_nonsys)
        AFsMapByFileName["nonsense"]=copy.deepcopy(AFintervalMap_SNPcounts_nonsense)
    ################## same to INTRON,UTR,INTERGENIC ############

    outfile=open(options.outfileprename+"DAF","w")
    print("DAFbin",end="\t",file=outfile)
    for tag in sorted(DAFsMapByFileName.keys()):
        print(tag,end="\t",file=outfile)
    print("",file=outfile)
    for a,b in DAFintervalMap_SNPcounts_template.keys():
        print(str(a),str(b),end="\t",file=outfile)
        for tag in sorted(DAFsMapByFileName.keys()):
            print(str(DAFsMapByFileName[tag][a,b]),end="\t",file=outfile)
        print("",file=outfile)
    outfile.close()
    outfile=open(options.outfileprename+"AAF","w")
    print("AFbin",end="\t",file=outfile)
    for tag in sorted(AFsMapByFileName.keys()):
        print(tag,end="\t",file=outfile)
    print("",file=outfile)
    for a,b in AFintervalMap_SNPcounts_template.keys():
        print(str(a),str(b),end="\t",file=outfile)
        for tag in sorted(AFsMapByFileName.keys()):
            print(str(AFsMapByFileName[tag][a,b]),end="\t",file=outfile)
        print("",file=outfile)
    outfile.close()
    outfile=open(options.outfileprename+"MAF","w")
    print("MAFbin",end="\t",file=outfile)
    for tag in sorted(MAFsMapByFileName.keys()):
        print(tag,end="\t",file=outfile)
    print("",file=outfile)
    for a,b in MAFintervalMap_SNPcounts_template.keys():
        print(str(a),str(b),end="\t",file=outfile)
        for tag in sorted(MAFsMapByFileName.keys()):
            print(str(MAFsMapByFileName[tag][a,b]),end="\t",file=outfile)
        print("",file=outfile)
    outfile.close()
    print("finished")