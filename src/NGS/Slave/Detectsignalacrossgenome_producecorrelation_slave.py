'''
Created on 2015-8-21

@author: liurui
'''
import fractions, re, os, copy,pysam
from optparse import OptionParser
from os.path import sys

from NGS.BasicUtil import Util, VCFutil
import NGS.BasicUtil.DBManager as dbm


parser = OptionParser()
parser.add_option("-T","--targetpopvcfconfig",dest="targetpopvcfconfig",action="append",help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-R","--refpopvcffileconfig",dest="refpopvcffileconfig",action="append",help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-t","--topleveltablejudgeancestral",dest="topleveltablejudgeancestral",help="R(r)/G(g)")
parser.add_option("-c","--chromlistfilename",dest="chromlistfilename")
parser.add_option("-n","--numberofindvdoftargetpop_todividintobin",dest="numberofindvdoftargetpop_todividintobin",default="o",help="conflit with correlationfile")
parser.add_option("-o","--outfileprewithpath",dest="outfileprewithpath")
(options, args) = parser.parse_args()
def make_freq_xaxisKEY_yaxisseqVALUERelation(a):
    chromlistfilename=a[0];topleveltablename=a[1];targetpopvcffile_withdepthconfig=a[2];refpopvcffile_withdepthconfig=a[3];numberofindvdoftargetpop_todividintobin=int(a[4])
    mindepthtojudefixed=20
    d_increase=fractions.Fraction(1, (2*int(numberofindvdoftargetpop_todividintobin)))
    d_increase=round(d_increase,11)
    minvalue=0.000000000000
    freq_xaxisKEY_yaxisVALUE_seq_list={}
    for i in range(numberofindvdoftargetpop_todividintobin*2-1):
        freq_xaxisKEY_yaxisVALUE_seq_list[(minvalue,minvalue+d_increase+0.00000000004)]=[]
        minvalue+=d_increase
    else:
        freq_xaxisKEY_yaxisVALUE_seq_list[(minvalue,1)]=[]
    for a,b in sorted(freq_xaxisKEY_yaxisVALUE_seq_list.keys()):
        print(str(a),str(b))
#     while minvalue+d_increase<=1:
#         freq_xaxisKEY_yaxisVALUE_seq_list[(minvalue,minvalue+d_increase+0.00000000004)]=[]
#         print('%.12f'%minvalue,'%.12f'%(minvalue+d_increase+0.00000000004))
#         minvalue+=d_increase
#     else:
#         freq_xaxisKEY_yaxisVALUE_seq_list[]
    print("process ID:",os.getpid(),"start",chromlistfilename) 
    dbvariantstools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.vcfdbname)
    chromlistfile=open(chromlistfilename,"r")
    chromlistfilelines=chromlistfile.readlines()
    chromlistfile.close()
    chromlist=[]
    for chrrow in chromlistfilelines:
        chrrowlist=re.split(r'\s+',chrrow.strip())
        chromlist.append((chrrowlist[0].strip(),int(chrrowlist[1].strip())))

    vcfnamelist=[];listofpopvcfmapOfAChr=[];methodlist=[]
    vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
    N_of_targetpop=len(targetpopvcffile_withdepthconfig)
    N_of_refpop=len(refpopvcffile_withdepthconfig)
    #{ vcftablename1:[depthfilename1,name1,name2] , vcftablename2:[depthfilename2,name1,name2] } or {vcftablename1:None, vcftablename2:None}
    for vcfconfigfilename in targetpopvcffile_withdepthconfig[:]+refpopvcffile_withdepthconfig[:]:
            listofpopvcfmapOfAChr.append({})
            vcfconfig=open(vcfconfigfilename,"r")
            for line in vcfconfig:
                vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                if vcffilename_obj!=None:
                    vcfname=vcffilename_obj.group(1).strip()
                    vcfnamelist.append(vcfname)
                    vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                    vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
                elif line.split():
                    vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
            vcfconfig.close()
            if re.search(r"indvd[^/]+",vcfname)!=None:
                methodlist.append("indvd")
    
            elif re.search(r"pool[^/]+",vcfname)!=None:
                methodlist.append("pool")
    
            else:
                print("vcfname must with 'pool' or 'indvd'")
                exit(-1)   
    for currentchrID,currentchrLen in chromlist:
        for vcfname in vcfnamelist:
            if currentchrID in vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][0].VcfIndexMap:
                break
        else:
            print("this chr doesn't exist in anypop")
            continue
        for vcfobj_idx in range(len(vcfnamelist)):
            listofpopvcfmapOfAChr[vcfobj_idx]={}
            listofpopvcfmapOfAChr[vcfobj_idx][currentchrID]=vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfnamelist[vcfobj_idx]][0].getVcfListByChrom(currentchrID)
        target_ref_SNPs=Util.alinmultPopSnpPos(listofpopvcfmapOfAChr, "o")
        for snp_aligned in target_ref_SNPs[currentchrID]:
            if len(snp_aligned[1])!=1 or len(snp_aligned[2])!=1:
                print("multple allele",snp_aligned)
                continue
            curpos=int(snp_aligned[0])
            snp=dbvariantstools.operateDB("select","select * from "+topleveltablename+" where chrID='"+currentchrID+"' and snp_pos="+str(curpos)+"")
            if not snp or snp==0:
                print(currentchrID,curpos,"snp not find in db,skip")
                continue
            else:#judge the ancenstrall allele
                fanyadepthlist=re.split(r",",snp[0][9])
                if len(fanyadepthlist)==2 and int(fanyadepthlist[1]) >=mindepthtojudefixed and fanyadepthlist[0].strip()=="0":
                    A_base_idx=1
                elif len(fanyadepthlist)==2 and int(fanyadepthlist[0])>=mindepthtojudefixed and fanyadepthlist[1].strip()=="0":
                    A_base_idx=0
                else:
                    print("skip snp",snp[0][1],snp[0][7:])
                    continue
            ancestrallcontext=snp[0][5].strip()[0].upper()+snp[0][3+A_base_idx].strip().upper()+snp[0][5].strip()[2].upper()
            if "CG" in ancestrallcontext or "GC" in ancestrallcontext:
                print("skip CG site",ancestrallcontext)
                continue
            ##########x-axis
            countedAF=0;target_DAF_sum=0#;noofnocoveredsample=0
            for i in range(3,N_of_targetpop+3):
                if snp_aligned[i]==None:
                    if len(vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfnamelist[i-3]])==1:
                        print("no depth file")
                        continue
                    else:
                        sum_depth=0
                        for samfile in vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfnamelist[i-3]][1:]:
                            ACGTdep=samfile.count_coverage(currentchrID,curpos-1,curpos)
                            for dep in ACGTdep:
                                sum_depth+=dep[0]
                        if sum_depth>=mindepthtojudefixed:
                            AF=0
                        else:
                            continue
                else:
                    if methodlist[i-3]=="indvd":
                        AF = float(re.search(r"AF=([\d\.]+);", snp_aligned[i][0]).group(1))
                    elif methodlist[i-3]=="pool":
                        refdep = 0;altalleledep = 0
                        AD_idx = (re.split(":", snp_aligned[i][1])).index("AD")  # gatk GT:AD:DP:GQ:PL
                        for sample in snp_aligned[i][2]:
                            if len(re.split(":", sample)) == 1:  # ./.
                                continue
                            AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                            try :
                                refdep += int(AD_depth[0])
                                altalleledep += int(AD_depth[1])
                            except ValueError:
                                print(sample, end="|")
                        if refdep==altalleledep and altalleledep==0:
                            print("no sample available in this pop")
#                                 noofnocoveredsample+=1
                            continue
                        AF=altalleledep/(altalleledep+refdep)
                if A_base_idx==0:
                    DAF=1-AF
                elif A_base_idx==1:
                    DAF=AF
                target_DAF_sum+=DAF;countedAF+=1
            if countedAF==0 :#or target_DAF_sum==0:
                print("skip this snp,because it fiexd as ancestral or no covered in this pos in target pops",snp_aligned,snp)
                continue
            target_DAF=target_DAF_sum/countedAF
            ###############y-axis
            countedAF=0;rer_DAF_sum=0
            for i in range(3+N_of_targetpop,N_of_refpop+N_of_targetpop+3):
                if snp_aligned[i]==None:
                    if len(vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfnamelist[i-3]])==1:
                        continue
                    else:
#                         depth_linelist=vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfnamelist[i-3-N_of_targetpop]].getdepthByPos_optimized(currentchrID,curpos)
                        sum_depth=0
                        for samfile in vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfnamelist[i-3]][1:]:
                            ACGTdep=samfile.count_coverage(currentchrID,curpos-1,curpos)
                            for dep in ACGTdep:
                                sum_depth+=dep[0]
                        if sum_depth>=mindepthtojudefixed:
                            AF=0
                        else:
                            continue
                else:
                    if methodlist[i-3]=="indvd":
                        AF = float(re.search(r"AF=([\d\.]+);", snp_aligned[i][0]).group(1))
                        AN = float(re.search(r"AN=([\d\.]+);", snp_aligned[i][0]).group(1))
                        if AN<5:
                            continue
                    elif methodlist[i-3]=="pool":
                        refdep = 0;altalleledep = 0
                        AD_idx = (re.split(":", snp_aligned[i][1])).index("AD")  # gatk GT:AD:DP:GQ:PL
                        for sample in snp_aligned[i][2]:
                            if len(re.split(":", sample)) == 1:  # ./.
                                continue
                            AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                            try :
                                refdep += int(AD_depth[0])
                                altalleledep += int(AD_depth[1])
                            except ValueError:
                                print(sample, end="|")
                        if (refdep==altalleledep and altalleledep==0) or altalleledep+ refdep<10:
                            continue
                        AF=altalleledep/(altalleledep+refdep)
                if A_base_idx==0:
                    DAF=1-AF
                elif A_base_idx==1:
                    DAF=AF
                rer_DAF_sum+=DAF;countedAF+=1
            if countedAF==0 or rer_DAF_sum==0:
                print("skip this snp,because it  no covered in this pos in ref pops",snp_aligned,snp)
                continue
            ######collect according bins
            for a,b in sorted(freq_xaxisKEY_yaxisVALUE_seq_list.keys()):
                if target_DAF >a and target_DAF<=b:
                    freq_xaxisKEY_yaxisVALUE_seq_list[(a,b)].append(rer_DAF_sum/countedAF)
                    break
#     freq_xaxisKEY_yaxisVALUERelation={}
#     for a,b in sorted(freq_xaxisKEY_yaxisVALUE_seq_list.keys()):
#         freq_xaxisKEY_yaxisVALUERelation[(a,b)]=numpy.mean(freq_xaxisKEY_yaxisVALUE_seq_list[(a,b)])
#         print('%.12f'%a,'%.12f'%(b),'%.12f'%(freq_xaxisKEY_yaxisVALUERelation[(a,b)]),"process ID:",os.getpid(),"done",sep="\t")
    print("process ID:",os.getpid(),"done")
    return copy.deepcopy(freq_xaxisKEY_yaxisVALUE_seq_list)
if __name__ == '__main__':
    filenamelistfilename=options.outfileprewithpath+".freqcorrelationfilenamelist"
    parameterstuples=(options.chromlistfilename,options.topleveltablejudgeancestral,options.targetpopvcfconfig,options.refpopvcffileconfig,options.numberofindvdoftargetpop_todividintobin)
    print(parameterstuples,options.outfileprewithpath)
    freq_xaxisKEY_yaxisVALUE_seq_list=make_freq_xaxisKEY_yaxisseqVALUERelation(parameterstuples)
    outfilename=options.outfileprewithpath+"_part_"+str(os.getpid())+Util.random_str()
    outfile=open(outfilename,'w')
    filenamelistfile=open(filenamelistfilename,'a')
    for a,b in sorted(freq_xaxisKEY_yaxisVALUE_seq_list.keys()):
        print(str(a),str(b),*freq_xaxisKEY_yaxisVALUE_seq_list[(a,b)],sep="\t",file=outfile)
    outfile.close()
    print(outfilename,file=filenamelistfile)
    print(sys.argv,outfilename)
    filenamelistfile.close()
    print("process ID:",os.getpid(),"finished")
    exit(0)