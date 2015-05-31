# -*- coding: UTF-8 -*-
'''
Created on 2015-5-28

@author: liurui
'''
import copy, re, os
from optparse import OptionParser
import random

from NGS.BasicUtil import Util
import src.NGS.BasicUtil.DBManager as dbm


parser = OptionParser()
parser.add_option("-i", "--interval", dest="interval", nargs=3,
                  help="minvalue maxvalue breaks", metavar="FILE")
parser.add_option("-w", "--wildcdsfilenames", dest="wildcdsfilenames",action="append",default=[],help="")
parser.add_option("-d", "--domesticcdsfilenames", dest="domesticcdsfilenames",action="append",default=[],help="")
parser.add_option("-t", "--topleveltablejudgeancestral", dest="topleveltablejudgeancestral")
parser.add_option("-o", "--outfileprename", dest="outfileprename", help="default infile1_infile2")
(options, args) = parser.parse_args()
mindeptojudgefix=15

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
for fname in options.wildcdsfilenames:
    os.system("awk 'NR>1{print $1}' "+fname+"|sort|uniq|sort >"+fname+"_chrom")
#     wildcdsfilelist.append(open(fname,'r'))
for fname in options.domesticcdsfilenames:
    os.system("awk 'NR>1{print $1}' "+fname+"|sort|uniq|sort >"+fname+"_chrom")
#     domesticcdsfilelist.append(open(fname,'r'))
#merge and uniq ,when there is only one file for wild or domestic,the for loop below don't excute
for f_idx in range(1,len(options.wildcdsfilenames)):
    os.system("cat "+options.wildcdsfilenames[f_idx-1]+"_chrom "+options.wildcdsfilenames[f_idx]+"_chrom|sort|uniq|sort > temp_chrom")
    os.system("rm "+options.wildcdsfilenames[f_idx-1]+"_chrom ")
    os.system("mv temp_chrom "+options.wildcdsfilenames[f_idx]+"_chrom")
for f_idx in range(1,len(options.domesticcdsfilenames)):
    os.system("cat "+options.domesticcdsfilenames[f_idx-1]+"_chrom "+options.domesticcdsfilenames[f_idx]+"_chrom|sort|uniq|sort > temp_chrom")
    os.system("rm "+options.domesticcdsfilenames[f_idx-1]+"_chrom ")
    os.system("mv temp_chrom "+options.domesticcdsfilenames[f_idx]+"_chrom")
os.system("comm -12 "+options.domesticcdsfilenames[-1]+"_chrom "+options.wildcdsfilenames[-1]+"_chrom > chromtable_containbothwildanddomestic")
os.system("rm "+options.domesticcdsfilenames[-1]+"_chrom");os.system("rm "+options.wildcdsfilenames[-1]+"_chrom")
t=open("chromtable_containbothwildanddomestic","r")
chromlist=t.readlines();t.close()
if __name__ == '__main__':
    intervalFileName=options.outfileprename+".interval"
    intervalfile=open(intervalFileName,"w")
    while minvalue + d_increase<= maxvalue :
        print(str(minvalue),str(minvalue+d_increase),sep="\t",file=intervalfile)
        minvalue+=d_increase
    else:
        if minvalue<maxvalue:
            print(str(minvalue),str(maxvalue),sep="\t",file=intervalfile)
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
        wildcdsfilelist=[]
        domesticcdsfilelist=[]
        curchrom=chrom.strip()
        for fname in options.wildcdsfilenames:
            os.system("rm "+fname+"_one_chrom")
            os.system("awk '$1~/"+curchrom+"/{print $0}' "+fname+">"+fname+"_one_chrom")
            a=os.popen("less -S "+fname+"_one_chrom|wc -l")
            if a.readline().strip()=="0":
                a.close()
                continue
            wildcdsfilelist.append(open(fname+"_one_chrom",'r'))
        for fname in options.domesticcdsfilenames:
            os.system("rm "+fname+"_one_chrom")
            os.system("awk '$1~/"+curchrom+"/{print $0}' "+fname+">"+fname+"_one_chrom")
            a=os.popen("less -S "+fname+"_one_chrom|wc -l")
            if a.readline().strip()=="0":
                a.close()
                continue
            domesticcdsfilelist.append(open(fname+"_one_chrom",'r'))
    ########## collect delta_AF #################################################################
        for wf_idx in range(len(wildcdsfilelist)):
#             wildcdsfilelist[wf_idx].readline()#title
            line=wildcdsfilelist[wf_idx].readline()
            if line.split():
                wild_CurRecsLinelist.append(re.split(r"\s+",line.strip()))
                posOfCurRecwild.append(int(wild_CurRecsLinelist[wf_idx][1]))###############
            else:
                wild_CurRecsLinelist.append(None)
        for df_idx in range(len(domesticcdsfilelist)):
#             domesticcdsfilelist[df_idx].readline()#title
            line=domesticcdsfilelist[df_idx].readline()
            if line.split():
                dom_CurRecsLinelist.append(re.split(r"\s+",line.strip()))
                posOfCurRecdom.append(int(dom_CurRecsLinelist[df_idx][1]))
            else:
                dom_CurRecsLinelist.append(None)
        
        while wild_CurRecsLinelist!=[None]*len(wildcdsfilelist) or dom_CurRecsLinelist!=[None]*len(domesticcdsfilelist):
            #loop every time clean dom_CurPosRecs ;wild_CurPosRecs it is used to caculate delta AF for every pos
            wild_CurPosRecs=[];dom_CurPosRecs=[]
            curpos=min(posOfCurRecwild+posOfCurRecdom)
            for wf_idx in range(len(wildcdsfilelist)):
                if  wild_CurRecsLinelist[wf_idx]==None :
                    if  posOfCurRecwild[wf_idx]<=curpos:
                        posOfCurRecwild[wf_idx]=999999999999999999999999999999
                elif int(wild_CurRecsLinelist[wf_idx][1])==curpos:
                    if wild_CurRecsLinelist[wf_idx][0]!=curchrom:
                        print(wild_CurRecsLinelist[wf_idx][0],curchrom)
                        exit(-1)
                    wild_CurPosRecs.append(copy.deepcopy(wild_CurRecsLinelist[wf_idx]))
                    line=wildcdsfilelist[wf_idx].readline()
                    if line.split():
                        wild_CurRecsLinelist[wf_idx]=re.split(r"\s+",line.strip())
                        posOfCurRecwild[wf_idx]=int(wild_CurRecsLinelist[wf_idx][1])
                    else:
                        wild_CurRecsLinelist[wf_idx]=None
            for df_idx in range(len(domesticcdsfilelist)):
                if dom_CurRecsLinelist[df_idx]==None :
                    if  posOfCurRecdom[df_idx]<=curpos:
                        posOfCurRecdom[df_idx]=9999999999999999999999999999999
                elif int(dom_CurRecsLinelist[df_idx][1])==curpos:
                    if dom_CurRecsLinelist[df_idx][0]!=curchrom:
                        print(curchrom,dom_CurRecsLinelist[df_idx][0])
                        exit(-1)
                    dom_CurPosRecs.append(copy.deepcopy(dom_CurRecsLinelist[df_idx]))
                    line=domesticcdsfilelist[df_idx].readline()
                    if line.split():
                        dom_CurRecsLinelist[df_idx]=re.split(r"\s+",line.strip())
                        posOfCurRecdom[df_idx]=int(dom_CurRecsLinelist[df_idx][1])
                    else:
                        dom_CurRecsLinelist[df_idx]=None
#             print(curpos)
#             print(posOfCurRecwild,posOfCurRecdom)          
            if wild_CurPosRecs==[] or dom_CurPosRecs==[]:
#                 print("wild_CurPosRecs",wild_CurPosRecs)
#                 print("dom_CurPosRecs",dom_CurPosRecs)
#                 print("dom_CurRecsLinelist",dom_CurRecsLinelist)
#                 print("wild_CurRecsLinelist",wild_CurRecsLinelist)
                continue
            
            w_af=0;d_af=0
            #determin derived allele
            snp=dbvariantstools.operateDB("select","select * from "+options.topleveltablejudgeancestral+" where chrID='"+curchrom+"' and snp_pos='"+str(curpos)+"'")

            if snp and snp[0][6]!=None and snp[0][6]!="no covered" and re.search(r'[\w\W]+[,][\w\W]+:\d+,\d+',snp[0][6])==None:
                archicpop=re.search(r'([ATCGNatcg]+):(\d+),(\d+)',snp[0][6])
                archic_base=archicpop.group(1).strip().upper()
            if not snp or re.search(r'[\w\W]+[,][\w\W]+:\d+,\d+',snp[0][6])!=None or (archicpop.group(2).strip()!='0' and archicpop.group(3).strip()!='0') or snp[0][6]=="no covered":
                A_base_idx=random.randint(0,1)
#                 print("can't judge which is derived allele,random select ref or alt as ancenstral")
            if not snp or re.search(r'[\w\W]+[,][\w\W]+:\d+,\d+',snp[0][6])!=None or (int(archicpop.group(2).strip())+int(archicpop.group(3).strip())<=mindeptojudgefix):
                A_base_idx=random.randint(0,1)
#                 print("can't judge which is derived allele,random select ref or alt as ancenstral")
            if snp and snp[0][6]!="no covered" and archicpop.group(2).strip()=='0':
#                 print("determin derived allele")
                A_base_idx=1#alt_allele is the ancestral allele
            elif snp and snp[0][6]!="no covered" and  archicpop.group(3).strip()=='0':
#                 print("determin derived allele")
                A_base_idx=0#ref_allele is the ancestral allele
            #caculate allele freq for each pop
            for e in wild_CurPosRecs:
                if A_base_idx==1:
                    w_af+=(1-float(e[5]))
                else:
                    w_af+=float(e[5])
            for e in dom_CurPosRecs:
                if A_base_idx==1:
                    d_af+=(1-float(e[5]))
                else:
                    d_af+=float(e[5])
            
            delta_af=(w_af/len(wild_CurPosRecs))-(d_af/len(dom_CurPosRecs))
#             print("wild_CurPosRecs",wild_CurPosRecs,"\n","dom_CurPosRecs",dom_CurPosRecs)
            for a,b in intervalMap_wild_SNPrec.keys():
                if delta_af>=a and delta_af<b and len(wild_CurPosRecs[0])>=14:
                    if wild_CurPosRecs[0][-3]==wild_CurPosRecs[0][-1]:
#                         print("sysnonymous")
                        intervalMap_wild_SNPrec[a,b]["sysnonymous"].append(wild_CurPosRecs)
                    elif wild_CurPosRecs[0][-3].find("*")!=-1 or wild_CurPosRecs[0][-1].find("*")!=-1:
#                         print("nonsense")
                        intervalMap_wild_SNPrec[a,b]["nonsense"].append(wild_CurPosRecs)
                    elif wild_CurPosRecs[0][-3]!=wild_CurPosRecs[0][-1]:
#                         print("nonsysnonymous")
                        intervalMap_wild_SNPrec[a,b]["nonsysnonymous"].append(wild_CurPosRecs)
            #reverse
            delta_af=(d_af/len(dom_CurPosRecs))-(w_af/len(wild_CurPosRecs))
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
#output    
    for a,b in intervalMap_wild_SNPrec.keys():
        for rec in intervalMap_wild_SNPrec[a,b]["sysnonymous"]:
            statisticMap[a,b][1]["sysnonymous"]+=1
            print(a,b,"sysnonymous",*rec,sep="\t",file=outwildfile)
        for rec in intervalMap_wild_SNPrec[a,b]["nonsysnonymous"]:
            statisticMap[a,b][1]["nonsysnonymous"]+=1
            print(a,b,"nonsysnonymous",*rec,sep="\t",file=outwildfile)
        for rec in intervalMap_wild_SNPrec[a,b]["nonsense"]:
            statisticMap[a,b][1]["nonsense"]+=1
            print(a,b,"nonsense",*rec,sep="\t",file=outwildfile)
    for a,b in intervalMap_dom_SNPrec.keys():
        for rec in intervalMap_dom_SNPrec[a,b]["sysnonymous"]:
            statisticMap[a,b][0]["sysnonymous"]+=1
            print(a,b,"sysnonymous",*rec,sep="\t",file=outdomesticfile)
        for rec in intervalMap_dom_SNPrec[a,b]["nonsysnonymous"]:
            statisticMap[a,b][0]["nonsysnonymous"]+=1
            print(a,b,"nonsysnonymous",*rec,sep="\t",file=outdomesticfile)
        for rec in intervalMap_dom_SNPrec[a,b]["nonsense"]:
            statisticMap[a,b][0]["nonsense"]+=1
            print(a,b,"nonsense",*rec,sep="\t",file=outdomesticfile)
    print("delta_AFbins\t\tdomesticpops\t\twildpops\t\t",file=statisticsfile)
    print("            \tsysnonymous\tnonsysnonymous\tnonsense\tsysnonymous\tnonsysnonymous\tnonsense",file=statisticsfile)
    for a,b in statisticMap.keys():
        print(a,b,statisticMap[a,b][0]["sysnonymous"],statisticMap[a,b][0]["nonsysnonymous"],statisticMap[a,b][0]["nonsense"],statisticMap[a,b][1]["sysnonymous"],statisticMap[a,b][1]["nonsysnonymous"],statisticMap[a,b][1]["nonsense"],sep="\t",file=statisticsfile)
    statisticsfile.close()
    outdomesticfile.close()
    outwildfile.close()
    print("finish")