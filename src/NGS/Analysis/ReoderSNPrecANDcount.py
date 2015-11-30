# -*- coding: UTF-8 -*-
'''
Created on 2015-5-28

@author: liurui
'''
import copy, re, os
from optparse import OptionParser
from scipy.stats import fisher_exact as fe

from NGS.BasicUtil import Util, VCFutil
import src.NGS.BasicUtil.DBManager as dbm


parser = OptionParser()
parser.add_option("-i", "--interval", dest="interval", nargs=3,
                  help="minvalue maxvalue breaks", metavar="FILE")
parser.add_option("-w", "--wildcdsfilenames", dest="wildcdsfilenames",action="append",nargs=2,help="cds depth")
parser.add_option("-d", "--domesticcdsfilenames", dest="domesticcdsfilenames",action="append",nargs=2,help="cds depth")
parser.add_option("-t", "--topleveltablejudgeancestral", dest="topleveltablejudgeancestral")
parser.add_option("-o", "--outfileprename", dest="outfileprename", help="default infile1_infile2")
parser.add_option("-2","--ancenstral_or_derived",dest="ancenstral_or_derived",default="d",help="ancenstral(a) or derived(d)")
(options, args) = parser.parse_args()
mindeptojudgefix=15
#####################
VCFobj={};vcfnameKEY_depthfilename_titlenameVALUE_tojudgeancestrall={};vcfnameKEY_depthobjVALUE_tojudgeancestral={}
VCFobj["wigeon"]=VCFutil.VCF_Data("/home/bioinfo/liurui/data/vcffiles/uniqmap/taihudomesticgoose/taihudomesticgoose.pool.withindel.vcf")
VCFobj["fanya"]=VCFutil.VCF_Data("/home/bioinfo/liurui/data/vcffiles/uniqmap/fanya/fanya._pool.withindel.vcf")
vcfnameKEY_depthfilename_titlenameVALUE_tojudgeancestrall["wigeon"]=Util.GATK_depthfile("/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth","/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth.index")#here is a temp trick not a error
vcfnameKEY_depthfilename_titlenameVALUE_tojudgeancestrall["fanya"]=Util.GATK_depthfile("/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth","/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth.index")
vcfnameKEY_depthobjVALUE_tojudgeancestral["wigeon"]=["/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth",9]
vcfnameKEY_depthobjVALUE_tojudgeancestral["fanya"]=["/home/bioinfo/liurui/data/depth/g_j_sm_k_l_y_f_w_pool/gjsmklyfw_gatk.depth",3]
####################################
dbvariantstools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.vcfdbname)
minvalue = float(options.interval[0])
maxvalue = float(options.interval[1])
breaks = int(options.interval[2])
d_increase = (maxvalue - minvalue) / breaks
outwildfileName=options.outfileprename+"wild_"+str(len(options.wildcdsfilenames))
outwildfile=open(outwildfileName,'w')
outdomesticfileName=options.outfileprename+"domestic_"+str(len(options.domesticcdsfilenames))
outdomesticfile=open(outdomesticfileName,'w')
statisticsfile=open(options.outfileprename+"_statistics","w")
#make a chrom list contain all chrom of all file,no matter there are one or more file for -w -d each
#extract chroms and uniq
AF_idxlistwild=[]
AF_idxlistdomestic=[]
wildcdsfilenames=[];wild_depthfilenames={}#{ vcftablename1:[depthfilename1,name1,name2] , vcftablename2:[depthfilename2,name1,name2] } or {vcftablename1:None, vcftablename2:None}
for vcffilename,namefile in options.wildcdsfilenames:
    wildcdsfilenames.append(vcffilename)
    
    if namefile.lower()!="none":
        wild_depthfilenames[vcffilename]=[]
        fp=open(namefile,'r')
        for line in fp:
            depthfile_obj=re.search(r"depthfilename=(.*)",line.strip())
            if depthfile_obj!=None:
                wild_depthfilenames[vcffilename].append(depthfile_obj.group(1).strip())
            elif line.split():
                wild_depthfilenames[vcffilename].append(line.strip())
        fp.close()
    else:
        pass
domesticcdsfilenames=[];dom_depthfilenames={}#{ vcftablename1:[depthfilename1,name1,name2] , vcftablename2:[depthfilename2,name1,name2] } or {vcftablename1:None, vcftablename2:None}
for vcffilename,namefile in options.domesticcdsfilenames:
    domesticcdsfilenames.append(vcffilename)
    
    if namefile.lower()!="none":
        dom_depthfilenames[vcffilename]=[]
        fp=open(namefile,'r')
        for line in fp:
            depthfile_obj=re.search(r"depthfilename=(.*)",line.strip())
            if depthfile_obj!=None:
                dom_depthfilenames[vcffilename].append(depthfile_obj.group(1).strip())
            elif line.split():
                dom_depthfilenames[vcffilename].append(line.strip())
        fp.close()
    else:
        pass
for fname in wildcdsfilenames:
    t=open(fname,'r')
    AF_idxlistwild.append(re.split(r"\t",t.readline().strip()).index("AF"))
    t.close()
    os.system("awk 'NR>1{print $1}' "+fname+"|sort|uniq|sort >"+fname+"_chrom")
print("AF_idxlistwild",AF_idxlistwild)
#     wildcdsfilelist.append(open(fname,'r'))
for fname in domesticcdsfilenames:
    t=open(fname,'r')
    AF_idxlistdomestic.append(re.split(r"\t",t.readline().strip()).index("AF"))
    t.close()
    os.system("awk 'NR>1{print $1}' "+fname+"|sort|uniq|sort >"+fname+"_chrom")
print("AF_idxlistdomestic",AF_idxlistdomestic)
#     domesticcdsfilelist.append(open(fname,'r'))
#merge and uniq ,when there is only one file for wild or domestic,the for loop below don't excute
randomstr=Util.random_str()
for f_idx in range(1,len(wildcdsfilenames)):
    os.system("cat "+wildcdsfilenames[f_idx-1]+"_chrom "+wildcdsfilenames[f_idx]+"_chrom|sort|uniq|sort > temp_chrom"+randomstr)
    os.system("rm "+wildcdsfilenames[f_idx-1]+"_chrom ")
    os.system("mv temp_chrom"+randomstr+" "+wildcdsfilenames[f_idx]+"_chrom")
for f_idx in range(1,len(domesticcdsfilenames)):
    os.system("cat "+domesticcdsfilenames[f_idx-1]+"_chrom "+domesticcdsfilenames[f_idx]+"_chrom|sort|uniq|sort > temp_chrom"+randomstr)
    os.system("rm "+domesticcdsfilenames[f_idx-1]+"_chrom ")
    os.system("mv temp_chrom"+randomstr+" "+domesticcdsfilenames[f_idx]+"_chrom")

os.system("comm -12 "+domesticcdsfilenames[-1]+"_chrom "+wildcdsfilenames[-1]+"_chrom > chromtable_containbothwildanddomestic"+randomstr)
os.system("rm "+domesticcdsfilenames[-1]+"_chrom");os.system("rm "+wildcdsfilenames[-1]+"_chrom")

t=open("chromtable_containbothwildanddomestic"+randomstr,"r")
chromlist=t.readlines();t.close()
os.system("rm chromtable_containbothwildanddomestic"+randomstr)
if __name__ == '__main__':
    intervalFileName=options.outfileprename+".interval"
    intervalfile=open(intervalFileName,"w")
    while minvalue + d_increase<= maxvalue-d_increase :
        print(str(minvalue),str(minvalue+d_increase),sep="\t",file=intervalfile)
        minvalue+=d_increase
    else:
        if minvalue<maxvalue-d_increase:
            print("never get here")
        print(str(minvalue),str(1.0),sep="\t",file=intervalfile)
        intervalfile.close()
    ###### data structure ###################################
    intervalMap_dom_SNPrec={}#{(bin_start,bin_end):{sys:[reclist1,reclist2,,,,],nonsys:[],nonsense:[]}
    intervalMap_wild_SNPrec={}
#     wild_CurRecsLinelist=[];dom_CurRecsLinelist=[]#wildcurline_rec_collector record current linelist of all files in the same order as wildcdsfilelist
#     wild_CurPosRecs=[];dom_CurPosRecs=[]#
#     posOfCurRecwild=[];posOfCurRecdom=[]
    statisticMap={}
    ################## init intervalMap_wild_SNPrec and intervalMap_dom_SNPrec##########
    intervalfile=open(intervalFileName,'r')
    for line in intervalfile:
        linelist=re.split(r'\s+',line.strip())
        statisticMap[float(linelist[0]),float(linelist[1])]=[{"sysnonymous":0,"nonsysnonymous":0,"nonsense":0},{"sysnonymous":0,"nonsysnonymous":0,"nonsense":0}]#domestic and wild  
        intervalMap_wild_SNPrec[float(linelist[0]),float(linelist[1])]={"sysnonymous":[],"nonsysnonymous":[],"nonsense":[]}
        intervalMap_dom_SNPrec[float(linelist[0]),float(linelist[1])]={"sysnonymous":[],"nonsysnonymous":[],"nonsense":[]}
        
    ############################# bin the delta AF for all chrom; ###########
    for chrom in chromlist:
        wild_CurRecsLinelist=[];dom_CurRecsLinelist=[]#wildcurline_rec_collector record current linelist of all files in the same order as wildcdsfilelist
        wild_CurPosRecs=[];dom_CurPosRecs=[]#
        posOfCurRecwild=[];posOfCurRecdom=[]
    ###################### prepare chrom specified cdsreds for wild and domestic ###################################
        wildcdsfilelist=[];wildcdsfilenamelist=[]
        domesticcdsfilelist=[];domesticcdsfilenamelist=[]
        curchrom=chrom.strip()
        allsnprecinAchr_mapbyvcfname={}
        
        idx=len(wildcdsfilenames)
        depthobjmap={};species_idx_map={};randomstr=Util.random_str()
        for fname in wildcdsfilenames[::-1]:#bacause popup need from the end to start according by idx
            depthobjmap[fname]=Util.GATK_depthfile(wild_depthfilenames[fname][0],wild_depthfilenames[fname][0]+".index")
            species_idx_map[fname]=[]
            for titlename in wild_depthfilenames[fname][1:]:
                species_idx_map[fname].append(depthobjmap[fname].title.index("Depth_for_"+titlename))
            idx-=1
#             os.system("rm "+fname+"_one_chrom")
            filename=fname+"_one_chrom"+randomstr
            os.system("awk '$1~/"+curchrom+"/{print $0}' "+fname+">"+filename)

            wildcdsfilelist.append(open(filename,'r'))
            wildcdsfilenamelist.append(filename)

        wildcdsfilelist.reverse()
#         print(tempAF_idxlistwild)
        idx=len(domesticcdsfilenames)
        for fname in domesticcdsfilenames[::-1]:
            depthobjmap[fname]=Util.GATK_depthfile(dom_depthfilenames[fname][0],dom_depthfilenames[fname][0]+".index")
            species_idx_map[fname]=[]
            for titlename in dom_depthfilenames[fname][1:]:
                species_idx_map[fname].append(depthobjmap[fname].title.index("Depth_for_"+titlename))
            idx-=1
#             os.system("rm "+fname+"_one_chrom")
            filename=fname+"_one_chrom"+randomstr
            os.system("awk '$1~/"+curchrom+"/{print $0}' "+fname+">"+filename)
#             a=os.popen("less -S "+filename+"|wc -l")
#             if a.readline().strip()=="0":
#                 a.close()
#                 tempAF_idxlistdomestic.pop(idx)
# #                 domesticcdsfilelist.append("norecords")
#                 continue
            domesticcdsfilelist.append(open(filename,'r'))
            domesticcdsfilenamelist.append(filename)
        domesticcdsfilelist.reverse()
#         print(tempAF_idxlistdomestic)
    ########## collect delta_AF #################################################################
        for wf_idx in range(len(wildcdsfilelist)):
#             if wildcdsfilelist[wf_idx]=="norecords":
#                 continue
#             wildcdsfilelist[wf_idx].readline()#title
            line=wildcdsfilelist[wf_idx].readline()
            if line.split():
                wild_CurRecsLinelist.append(re.split(r"\s+",line.strip()))
                posOfCurRecwild.append(int(wild_CurRecsLinelist[wf_idx][1]))###############
            else:
                wild_CurRecsLinelist.append("endline")
                posOfCurRecwild.append(999999999999999999999999999)
        for df_idx in range(len(domesticcdsfilelist)):
#             domesticcdsfilelist[df_idx].readline()#title
#             if domesticcdsfilelist[df_idx]=="norecords":
#                 continue
            line=domesticcdsfilelist[df_idx].readline()
            if line.split():
                dom_CurRecsLinelist.append(re.split(r"\s+",line.strip()))
                posOfCurRecdom.append(int(dom_CurRecsLinelist[df_idx][1]))
            else:
                dom_CurRecsLinelist.append("endline")
                posOfCurRecdom.append(9999999999999999999999999999)
        
        while wild_CurRecsLinelist!=["endline"]*len(wildcdsfilelist) or dom_CurRecsLinelist!=["endline"]*len(domesticcdsfilelist):
            #loop every time clean dom_CurPosRecs ;wild_CurPosRecs it is used to caculate delta AF for every pos
            wild_CurPosRecs=[];dom_CurPosRecs=[]
#             AF_idxlistwild_CurPos=[];AF_idxlistdomestic_CurPos=[]
            curpos=min(posOfCurRecwild+posOfCurRecdom)
            for wf_idx in range(len(wildcdsfilelist)):
                if  wild_CurRecsLinelist[wf_idx]=="endline" or int(wild_CurRecsLinelist[wf_idx][1])>curpos:#None means reach the end of the line
                    if  posOfCurRecwild[wf_idx]<=curpos:#endline 
                        posOfCurRecwild[wf_idx]=999999999999999999999999999999
                    # 
                    depth_linelist=depthobjmap[wildcdsfilenames[wf_idx]].getdepthByPos_optimized(curchrom,curpos)
                    sum_depth=0
                    for idx in species_idx_map[wildcdsfilenames[wf_idx]]:
                        sum_depth+=int(depth_linelist[idx])
                    wild_CurPosRecs.append(["unknow"]*(AF_idxlistwild[wf_idx]+1))
                    if sum_depth>mindeptojudgefix:
                        wild_CurPosRecs[wf_idx][AF_idxlistwild[wf_idx]]=0
                elif int(wild_CurRecsLinelist[wf_idx][1])==curpos:
                    wild_CurPosRecs.append(copy.deepcopy(wild_CurRecsLinelist[wf_idx]))
#                     AF_idxlistwild_CurPos.append(AF_idxlistwild[wf_idx])
                    line=wildcdsfilelist[wf_idx].readline()
                    if line.split():
                        wild_CurRecsLinelist[wf_idx]=re.split(r"\s+",line.strip())
                        posOfCurRecwild[wf_idx]=int(wild_CurRecsLinelist[wf_idx][1])
                    else:
                        wild_CurRecsLinelist[wf_idx]="endline"
            #
            for df_idx in range(len(domesticcdsfilelist)):
                if dom_CurRecsLinelist[df_idx]=="endline" or int(dom_CurRecsLinelist[df_idx][1])>curpos:
                    if  posOfCurRecdom[df_idx]<=curpos:
                        posOfCurRecdom[df_idx]=9999999999999999999999999999999
                    depth_linelist=depthobjmap[domesticcdsfilenames[df_idx]].getdepthByPos_optimized(curchrom,curpos)
                    sum_depth=0
                    for idx in species_idx_map[domesticcdsfilenames[df_idx]]:
                        sum_depth+=int(depth_linelist[idx])
                    dom_CurPosRecs.append(["unknow"]*(AF_idxlistdomestic[df_idx]+1))
                    if sum_depth>mindeptojudgefix:
                        dom_CurPosRecs[df_idx][AF_idxlistdomestic[df_idx]]=0
                elif int(dom_CurRecsLinelist[df_idx][1])==curpos:
                    dom_CurPosRecs.append(copy.deepcopy(dom_CurRecsLinelist[df_idx]))
#                     AF_idxlistdomestic_CurPos.append(AF_idxlistdomestic[df_idx])
                    line=domesticcdsfilelist[df_idx].readline()
                    if line.split():
                        dom_CurRecsLinelist[df_idx]=re.split(r"\s+",line.strip())
                        posOfCurRecdom[df_idx]=int(dom_CurRecsLinelist[df_idx][1])
                    else:
                        dom_CurRecsLinelist[df_idx]="endline"
#             print(curpos)
#             print(posOfCurRecwild,posOfCurRecdom)          
            if wild_CurPosRecs==[] or dom_CurPosRecs==[]:

                continue
            
            w_af=0;d_af=0
            #determin derived allele
            snp=dbvariantstools.operateDB("select","select * from "+options.topleveltablejudgeancestral+" where chrID='"+curchrom+"' and snp_pos='"+str(curpos)+"'")
            if not snp:
                if allsnprecinAchr_mapbyvcfname=={}:
                    allsnprecinAchr_mapbyvcfname["wigeon"]=VCFobj["wigeon"].getVcfListByChrom(curchrom)
                    allsnprecinAchr_mapbyvcfname["fanya"]=VCFobj["wigeon"].getVcfListByChrom(curchrom)
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
                        if refdep!=0  or altalleledep<mindeptojudgefix:
                            A_base_idx=-1#alt_allele is the ancestral allele
                            low0=high0+1
                        elif refdep==0 and altalleledep>=mindeptojudgefix:
                            A_base_idx=1
                            break
                        
                else:
                    depth_linelist=vcfnameKEY_depthfilename_titlenameVALUE_tojudgeancestrall["fanya"].getdepthByPos_optimized(curchrom,curpos)
                    if int(depth_linelist[vcfnameKEY_depthobjVALUE_tojudgeancestral["fanya"][1]])>=mindeptojudgefix:
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
#                 if snp[0][11]==None or snp[0][7]==None:
#                     continue
                wigeondepthlist1=re.split(r",",snp[0][13])
                fanyadepthlist=re.split(r",",snp[0][9])
#                 barheaddepthlist=re.split(r",",snp[0][11])
#                 if len(fanyadepthlist)==2 and int(fanyadepthlist[1]) >=mindeptojudgefix and fanyadepthlist[0].strip()=="0":
#                     A_base_idx=1
#                 elif len(fanyadepthlist)==2 and int(fanyadepthlist[0])>=mindeptojudgefix and fanyadepthlist[1].strip()=="0":
#                     A_base_idx=0
#                 else:
#                     print("skip snp",snp[0][1],snp[0][7:])
#                     continue
                #############
                if len(wigeondepthlist1)==2 and len(fanyadepthlist)==2 and (int(wigeondepthlist1[0]) + int(wigeondepthlist1[1])>=mindeptojudgefix or int(fanyadepthlist[0]) + int(fanyadepthlist[1])>=mindeptojudgefix) and ((wigeondepthlist1[0].strip()=="0" and fanyadepthlist[0].strip()=="0") or (wigeondepthlist1[1].strip()=="0" and fanyadepthlist[1].strip()=="0") ):
                    if wigeondepthlist1[0].strip()=="0" and fanyadepthlist[0].strip()=="0":
                        A_base_idx=1
                    elif wigeondepthlist1[1].strip()=="0" and fanyadepthlist[1].strip()=="0":
                        A_base_idx=0
                    else:
                        print(snp,"never get here!")
                elif (len(fanyadepthlist)==2  and int(fanyadepthlist[0]) + int(fanyadepthlist[1])>=mindeptojudgefix+5 and (fanyadepthlist[0].strip()=="0" or fanyadepthlist[1].strip()=="0" )):#   or (snp[0][7]=="no covered" and len(depthlist2)==2 and int(depthlist2[0]) + int(depthlist2[1])>=mindeptojudgefix and (depthlist2[1].strip()=="0" or depthlist2[0].strip()=="0")):
                    if  fanyadepthlist[0].strip()=="0":
                        A_base_idx=1
                    elif fanyadepthlist[1].strip()=="0":
                        A_base_idx=0
                    else:
                        print(snp,"never get here!")
                else:
                    print(snp,"skip snp")
                    continue
            idx=0;w_unknowcount=0
            if options.ancenstral_or_derived=="a":
                A_base_idx=1-A_base_idx
            for e in wild_CurPosRecs:
                if e[AF_idxlistwild[idx]]=="unknow":
                    idx+=1;w_unknowcount+=1
                    continue
                if A_base_idx==1:
                    w_af+=(1-float(e[AF_idxlistwild[idx]]))
                else:
                    w_af+=float(e[AF_idxlistwild[idx]])
                idx+=1
            idx=0;d_unknowcount=0
            for e in dom_CurPosRecs:
                if e[AF_idxlistdomestic[idx]]=="unknow":
                    idx+=1;d_unknowcount+=1
                    continue
                if A_base_idx==1:
                    d_af+=(1-float(e[AF_idxlistdomestic[idx]]))
                else:
                    d_af+=float(e[AF_idxlistdomestic[idx]])
                
                idx+=1
            if len(wild_CurPosRecs)-w_unknowcount==0 or len(dom_CurPosRecs)-d_unknowcount==0:
                continue
            delta_af=(w_af/(len(wild_CurPosRecs)-w_unknowcount))-(d_af/(len(dom_CurPosRecs)-d_unknowcount))
#             print("wild_CurPosRecs",wild_CurPosRecs,"\n","dom_CurPosRecs",dom_CurPosRecs)
#            count by bin
            for a,b in intervalMap_wild_SNPrec.keys():
                if delta_af>=a and delta_af<b and len(wild_CurPosRecs[0])>=14:#so we can judge it is sysnonymous or nonsysnonymous or nonsense
                    if wild_CurPosRecs[0][-3]==wild_CurPosRecs[0][-1]:
#                         print("sysnonymous")
                        intervalMap_wild_SNPrec[a,b]["sysnonymous"].append(wild_CurPosRecs)
                    elif wild_CurPosRecs[0][-3].find("*")!=-1 or wild_CurPosRecs[0][-1].find("*")!=-1:
                        
                        intervalMap_wild_SNPrec[a,b]["nonsense"].append(wild_CurPosRecs)
                    elif wild_CurPosRecs[0][-3]!=wild_CurPosRecs[0][-1]:
#                         print(w_af,"=",wild_CurPosRecs,d_af,"=",dom_CurPosRecs,delta_af,">",a,delta_af,"<",b)
                        intervalMap_wild_SNPrec[a,b]["nonsysnonymous"].append(wild_CurPosRecs)
            #reverse
            delta_af=(d_af/(len(dom_CurPosRecs)-d_unknowcount))-(w_af/(len(wild_CurPosRecs)-w_unknowcount))
            for a,b in intervalMap_dom_SNPrec.keys():
                if delta_af>=a and delta_af<b and len(dom_CurPosRecs[0])>=14:
                    if dom_CurPosRecs[0][-3]==dom_CurPosRecs[0][-1]:
                        intervalMap_dom_SNPrec[a,b]["sysnonymous"].append(dom_CurPosRecs)
                    elif dom_CurPosRecs[0][-3].find("*")!=-1 or dom_CurPosRecs[0][-1].find("*")!=-1:
                        intervalMap_dom_SNPrec[a,b]["nonsense"].append(dom_CurPosRecs)
                    elif dom_CurPosRecs[0][-3]!=dom_CurPosRecs[0][-1]:
                        intervalMap_dom_SNPrec[a,b]["nonsysnonymous"].append(dom_CurPosRecs)

        for f in domesticcdsfilelist:
            f.close()
        for f in wildcdsfilelist:
            f.close()
        for filename in domesticcdsfilenamelist:
            os.system("rm "+filename)
        for filename in wildcdsfilenamelist:
            os.system("rm "+filename)
#output    
    for a,b in intervalMap_wild_SNPrec.keys():
        for rec in intervalMap_wild_SNPrec[a,b]["sysnonymous"]:
            statisticMap[a,b][1]["sysnonymous"]+=1
            print(round(a,3),round(b,3),"sysnonymous",*rec,sep="\t",file=outwildfile)
        for rec in intervalMap_wild_SNPrec[a,b]["nonsysnonymous"]:
            statisticMap[a,b][1]["nonsysnonymous"]+=1
            print(round(a,3),round(b,3),"nonsysnonymous",*rec,sep="\t",file=outwildfile)
        for rec in intervalMap_wild_SNPrec[a,b]["nonsense"]:
            statisticMap[a,b][1]["nonsense"]+=1
            print(round(a,3),round(b,3),"nonsense",*rec,sep="\t",file=outwildfile)
    for a,b in intervalMap_dom_SNPrec.keys():
        for rec in intervalMap_dom_SNPrec[a,b]["sysnonymous"]:
            statisticMap[a,b][0]["sysnonymous"]+=1
            print(round(a,3),round(b,3),"sysnonymous",*rec,sep="\t",file=outdomesticfile)
        for rec in intervalMap_dom_SNPrec[a,b]["nonsysnonymous"]:
            statisticMap[a,b][0]["nonsysnonymous"]+=1
            print(round(a,3),round(b,3),"nonsysnonymous",*rec,sep="\t",file=outdomesticfile)
        for rec in intervalMap_dom_SNPrec[a,b]["nonsense"]:
            statisticMap[a,b][0]["nonsense"]+=1
            print(round(a,3),round(b,3),"nonsense",*rec,sep="\t",file=outdomesticfile)
    print("delta_AFbins\t\tdomesticpops\t\twildpops\t\t",file=statisticsfile)
    print("            \tsysnonymous\tnonsysnonymous\tnonsense\tsysnonymous\tnonsysnonymous\tnonsense\tfishertestpvalue",file=statisticsfile)
    for a,b in sorted(statisticMap.keys()):
        print(round(a,3),round(b,3),statisticMap[a,b][0]["sysnonymous"],statisticMap[a,b][0]["nonsysnonymous"],statisticMap[a,b][0]["nonsense"],statisticMap[a,b][1]["sysnonymous"],statisticMap[a,b][1]["nonsysnonymous"],statisticMap[a,b][1]["nonsense"],str(fe([[statisticMap[a,b][0]["nonsysnonymous"],statisticMap[a,b][1]["nonsysnonymous"]],[statisticMap[a,b][0]["sysnonymous"],statisticMap[a,b][1]["sysnonymous"]]],"two-sided")[1]),sep="\t",file=statisticsfile)
    statisticsfile.close()
    outdomesticfile.close()
    outwildfile.close()
    print("finish")