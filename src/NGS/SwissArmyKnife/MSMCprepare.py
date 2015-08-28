# -*- coding: UTF-8 -*-
'''
Created on 2015-8-10

@author: liurui
'''

from itertools import combinations
from optparse import OptionParser
import re, sys, copy, os

# from NGS.BasicUtil import *
# import src.NGS.BasicUtil.DBManager as dbm


parser = OptionParser()
parser.add_option("-b","--bamlist_eachpop",dest="bamlist_eachpop",action="append",help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-c","--configfile",dest="configfile",help="R(r)/G(g)")
parser.add_option("-m","--minlength",dest="minlength")
parser.add_option("-C","--chromlistfilename",dest="chromlistfilename")
(options, args) = parser.parse_args()
# chromtable = Util.pekingduckchromtable
primaryID = "chrID"
configfile=open(options.configfile,'r')
for line in configfile:
    if re.search(r"snpablegenome\s*=\s*(.*)",line.strip())!=None:
        snpablegenome=re.search(r"snpablegenome\s*=\s*(.*)",line.strip()).group(1)
    if re.search(r"genome\s*=\s*(.*)",line.strip())!=None:
        refgenome=re.search(r"genome\s*=\s*(.*)",line.strip()).group(1)
    if re.search(r"outputpath\s*=\s*(.*)",line.strip())!=None:
        outputpath=re.search(r"outputpath\s*=\s*(.*)",line.strip()).group(1)
    if re.search(r"apply_mask_l\s*=\s*(.*)",line.strip())!=None:
        apply_mask_l=re.search(r"apply_mask_l\s*=\s*(.*)",line.strip()).group(1)
    if re.search(r"samtoolspath\s*=\s*(.*)",line.strip())!=None:
        samtoolspath=re.search(r"samtoolspath\s*=\s*(.*)",line.strip()).group(1)
bamfilehandlers=[]
samplebamlistVALUEpopfilenameKEY={}
result_chrlistMAPbypopfilename={}
mappabilitychrlist=[]#order by chromID
allkindofpaire=list(combinations(options.bamlist_eachpop,2))
print(allkindofpaire)
msmsinfile_chrVALUE_popKEY={}
for popfilename in options.bamlist_eachpop:
    msmsinfile_chrVALUE_popKEY[popfilename]=[]
    bamfilehandlers.append(open(popfilename.strip(),'r'))
    samplebamlistVALUEpopfilenameKEY[popfilename.strip()]=[]
    result_chrlistMAPbypopfilename[popfilename.strip()]={}

for idx in range(len(options.bamlist_eachpop)):
    for bamfile in bamfilehandlers[idx]:
        if bamfile.startswith('#') or not bamfile.split():
            continue
        else:
            result_chrlistMAPbypopfilename[options.bamlist_eachpop[idx]][bamfile.strip()]=[]
            samplebamlistVALUEpopfilenameKEY[options.bamlist_eachpop[idx]].append(bamfile.strip())
            if not os.path.isfile(bamfile.strip()+".bai"):
                a=os.system(samtoolspath+"/samtools index "+bamfile.strip())
                if a!=0:
                    print(samtoolspath+"/samtools index "+bamfile.strip(),"index wrong")
                    exit(-1)
msmsinfilesplit_chrVALUE_popsKEY={}
for pop1,pop2 in allkindofpaire:
    msmsinfilesplit_chrVALUE_popsKEY[(pop1,pop2)]=[]
if __name__ == '__main__':
#     genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
#     currentsql="select * from " + chromtable+" where chrlength>="+options.minlength+" order by "+primaryID
#     result=genomedbtools.operateDB("select",currentsql)
    chromlistfile=open(options.chromlistfilename,"r")
    chromlist=[]
    for chrrow in chromlistfile:
        chrrowlist=re.split(r'\s+',chrrow.strip())
        chromlist.append((chrrowlist[0].strip(),int(chrrowlist[1].strip())))
    chromlistfile.close()

#     for i in range(0,totalChroms,20):
    for currentID,currentchrLen in chromlist:
#     for row in result:
#         print(row)
#         currentID=row[0].strip()
#         currentchrLen=row[1]
        a=os.system("""awk 'BEGIN{firstfind="false"}{if(firstfind=="false"){c=$0}if(c~/^>"""+currentID+"""/){firstfind="true";if($0!~/^>"""+currentID+"""/ && $0~/^>/){c=$0}else{print $0}}}' """+snpablegenome+""" |awk 'BEGIN{OFS=" ";FS=" "}{if(NR==1){sub(/[.]1/,"",$1)};print $0}' >"""+outputpath+"""/"""+currentID[:-2]+"""currentchrsnpablegenome""")
        if a!=0:
            print("""awk 'BEGIN{firstfind="false"}{if(firstfind=="false"){c=$0}if(c~/^>"""+currentID+"""/){firstfind="true";if($0!~/^>"""+currentID+"""/ && $0~/^>/){c=$0}else{print $0}}}' """+snpablegenome+""">"""+outputpath+"""/currentchrsnpablegenome""")
            exit(-1)
        print("finish awk produce currentchrsnpablegenome")
        a=os.system("""for((i=0,j=1;i<"""+str(currentchrLen)+""";i++,j++)); do printf """+currentID[:-2]+""""\t"$i"\t"$j"\n"; done >"""+ outputpath+"""/"""+currentID[:-2]+"""duck_currentchrom_position_raw_converse.bed """)
        if a!=0:
            print("""for((i=0,j=1;i<"""+str(currentchrLen)+""";i++,j++)); do printf """+currentID[:-2]+""""\t"$i"\t"$j"\n"; done >"""+ outputpath+"""/"""+currentID[:-2]+"""duck_currentchrom_position_raw_converse.bed """)
            exit(-1)
        print("finish for loop produce duck_currentchrom_position_raw_converse.bed")
        a=os.system(apply_mask_l+" "+outputpath+"/"+currentID[:-2]+"currentchrsnpablegenome  "+ outputpath+"/"+currentID[:-2]+"duck_currentchrom_position_raw_converse.bed |gzip -c >"+outputpath+"/mappability_mask"+currentID[:-2]+".bed.gz")
        mappabilitychrlist.append("mappability_mask"+currentID[:-2]+".bed.gz")
        if a!=0:
            print(apply_mask_l+" "+outputpath+"/currentchrsnpablegenome  "+ outputpath+"/duck_currentchrom_position_raw_converse.bed ")
            exit(-1)
        print(currentID[:-2])
        for pop in samplebamlistVALUEpopfilenameKEY.keys():
            print("popfilename:",pop)
            for pathtosample in samplebamlistVALUEpopfilenameKEY[pop]:
                samplename=re.search(r"[^/]*$",pathtosample).group(0)
                print("samplename",samplename)
                f=os.popen(samtoolspath+"/samtools depth -r "+currentID+" "+pathtosample+" | awk '{sum += $3} END {print sum / NR}'")
                meandepth=f.readline().strip()
                print("meandepth",meandepth)
                f.close()
#                 print(samtoolspath+"/samtools mpileup -q 20 -Q 20 -C 50 -u -r "+currentID+" -f "+refgenome+" "+pathtosample+" |"+samtoolspath+"/bcftools view -cgI - |sed 's/^"+currentID+"/"+currentID[:-2]+"/g'|/home/bioinfo/liurui/software/Python-3.4.3/python /pub/tool/msmc/tools/bamCaller.py "+meandepth+" "+outputpath+"/"+samplename+"_covered_sites_"+currentID[:-2]+".bed.gz |gzip -c > "+outputpath+"/"+samplename+"_"+currentID[:-2]+".vcf.gz")
                a=os.system(samtoolspath+"/samtools mpileup -q 20 -Q 20 -C 50 -u -r "+currentID+" -f "+refgenome+" "+pathtosample+" |"+samtoolspath+"/bcftools view -cgI - |sed 's/^"+currentID+"/"+currentID[:-2]+"/g'|/home/bioinfo/liurui/software/Python-3.4.3/python /pub/tool/msmc/tools/bamCaller.py "+meandepth+" "+outputpath+"/"+samplename+"_covered_sites_"+currentID[:-2]+".bed.gz |gzip -c > "+outputpath+"/"+samplename+"_"+currentID[:-2]+".vcf.gz" )
                if a!=0:
                    print("skip this chrom",currentID)
                    exit(-1)
                print("samtools finished")
                result_chrlistMAPbypopfilename[pop][pathtosample].append((samplename+"_covered_sites_"+currentID[:-2]+".bed.gz  ",samplename+"_"+currentID[:-2]+".vcf.gz  "))
        generate_multihetsep_statepart1="/home/bioinfo/liurui/software/Python-3.4.3/python /pub/tool/msmc/tools/generate_multihetsep.py "
        generate_multihetsep_statepart2=" "
        for pop in samplebamlistVALUEpopfilenameKEY.keys():
            generate_multihetsep_statepart1="/home/bioinfo/liurui/software/Python-3.4.3/python /pub/tool/msmc/tools/generate_multihetsep.py "
            generate_multihetsep_statepart2=" "
            popname=re.search(r"[^/]*$",pop).group(0).strip()
            if not os.path.exists(outputpath+"/"+popname):
                os.makedirs(outputpath+"/"+popname)
            for pathtosample in samplebamlistVALUEpopfilenameKEY[pop]:
                generate_multihetsep_statepart1+=" --mask="+outputpath+"/"+result_chrlistMAPbypopfilename[pop][pathtosample][-1][0]
                generate_multihetsep_statepart2+=outputpath+"/"+result_chrlistMAPbypopfilename[pop][pathtosample][-1][1]
            a=os.system(generate_multihetsep_statepart1+" --mask="+outputpath+"/mappability_mask"+currentID[:-2]+".bed.gz "+generate_multihetsep_statepart2+" > "+outputpath+"/"+popname+"/"+popname+"_"+currentID[:-2]+".msmc.infile")
            msmsinfile_chrVALUE_popKEY[pop].append(outputpath+"/"+popname+"/"+popname+"_"+currentID[:-2]+".msmc.infile  ")
            if a!=0:
                print(generate_multihetsep_statepart1+" --mask="+outputpath+"/mappability_mask"+currentID[:-2]+".bed.gz "+generate_multihetsep_statepart2+" > "+outputpath+"/"+popname+"/"+popname+"_"+currentID[:-2]+".msmc.infile")
                print(samtoolspath+"/samtools mpileup -q 20 -Q 20 -C 50 -u -r "+currentID+" -f "+refgenome+" "+pathtosample+" |"+samtoolspath+"/bcftools view -cgI - |sed 's/^"+currentID+"/"+currentID[:-2]+"/g'|/home/bioinfo/liurui/software/Python-3.4.3/python /pub/tool/msmc/tools/bamCaller.py "+meandepth+" "+outputpath+"/"+samplename+"_covered_sites_"+currentID[:-2]+".bed.gz |gzip -c > "+outputpath+"/"+samplename+"_"+currentID[:-2]+".vcf.gz")
                print("error generate_multihetsep.py")
                os.system("rm "+outputpath+"/"+popname+"/"+popname+"_"+currentID[:-2]+".msmc.infile")
                continue
                exit(-1)
#         for pop1,pop2 in allkindofpaire:
#             pop1name=re.search(r"[^/]*$",pop1).group(0).strip()
#             pop1sample1=samplebamlistVALUEpopfilenameKEY[pop1][0]
#             pop2name=re.search(r"[^/]*$",pop2).group(0).strip()
#             pop2sample1=samplebamlistVALUEpopfilenameKEY[pop2][0]
#             if not os.path.exists(outputpath+"/"+pop1name+"_"+pop2name):
#                 os.makedirs(outputpath+"/"+pop1name+"_"+pop2name)
#             print("/home/bioinfo/liurui/software/Python-3.4.3/python /pub/tool/msmc/tools/generate_multihetsep.py --mask="+outputpath+"/"+result_chrlistMAPbypopfilename[pop1][pop1sample1][-1][0]+" --mask="+outputpath+"/"+result_chrlistMAPbypopfilename[pop2][pop2sample1][-1][0]+" --mask="+outputpath+"/mappability_mask"+currentID[:-2]+".bed.gz "+outputpath+"/"+result_chrlistMAPbypopfilename[pop1][pop1sample1][-1][1]+outputpath+"/"+result_chrlistMAPbypopfilename[pop2][pop2sample1][-1][1]+" > "+outputpath+"/"+pop1name+"_"+pop2name+"/"+pop1name+"_"+pop2name+"_"+currentID[:-2]+".msmc.infile")
#             a=os.system("/home/bioinfo/liurui/software/Python-3.4.3/python /pub/tool/msmc/tools/generate_multihetsep.py --mask="+outputpath+"/"+result_chrlistMAPbypopfilename[pop1][pop1sample1][-1][0]+" --mask="+outputpath+"/"+result_chrlistMAPbypopfilename[pop2][pop2sample1][-1][0]+" --mask="+outputpath+"/mappability_mask"+currentID[:-2]+".bed.gz "+outputpath+"/"+result_chrlistMAPbypopfilename[pop1][pop1sample1][-1][1]+outputpath+"/"+result_chrlistMAPbypopfilename[pop2][pop2sample1][-1][1]+" > "+outputpath+"/"+pop1name+"_"+pop2name+"/"+pop1name+"_"+pop2name+"_"+currentID[:-2]+".msmc.infile")
#             msmsinfilesplit_chrVALUE_popsKEY[(pop1,pop2)].append(outputpath+"/"+pop1name+"_"+pop2name+"/"+pop1name+"_"+pop2name+"_"+currentID[:-2]+".msmc.infile  ")
    print(result_chrlistMAPbypopfilename)
    runmsmc="/pub/tool/msmc/msmc_linux_64bit  --fixedRecombination  -t 16 -o "
    print("all msmc prepare done,go into runing msmc")
#     for pop in msmsinfile_chrVALUE_popKEY.keys():
#         popname=re.search(r"[^/]*$",pop).group(0).strip()
#         runmsmc="/pub/tool/msmc/msmc_linux_64bit  --fixedRecombination  -t 16 -o "+outputpath+"/"+popname+"/"+popname+"result "
#         for msmcinfile_a_chr in msmsinfile_chrVALUE_popKEY[pop]:
#             runmsmc+=msmcinfile_a_chr
#         print(runmsmc)
#         a=os.system(runmsmc)
#         if a!=0:
#             print("error")
#             exit(-1)
#     #the command blow are modyfied after programm running,so break execution when programm reach here
#     for pop1,pop2 in msmsinfilesplit_chrVALUE_popsKEY.keys():
#         pop1name=re.search(r"[^/]*$",pop1).group(0).strip()
#         pop2name=re.search(r"[^/]*$",pop2).group(0).strip()
#         runmsmc="/pub/tool/msmc/msmc_linux_64bit  --fixedRecombination --skipAmbiguous -P 0,0,1,1 -t 16 -o "+outputpath+"/"+pop1name+"_"+pop2name+"/"+pop1name+"_"+pop2name+"result "
#         for msmcinfile_a_chr in msmsinfilesplit_chrVALUE_popsKEY[(pop1,pop2)]:
#             runmsmc+=msmcinfile_a_chr
#         print(runmsmc)
#         a=os.system(runmsmc)
#         if a!=0:
#             print("error")
#             exit(-1)