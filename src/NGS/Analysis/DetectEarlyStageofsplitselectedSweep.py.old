'''
Created on 2015-8-1

@author: liurui
'''
from optparse import OptionParser
import re,numpy,fractions,copy

from NGS.BasicUtil import *
import NGS.BasicUtil.Util
import src.NGS.BasicUtil.DBManager as dbm


primaryID = "chrID"
mindeptojudgefix=20
parser = OptionParser()
# parser.add_option("-c", "--chromtable", dest="chromtable",# action="callback",type="string",callback=useoptionvalue_previous2,
#                   help="write report to FILE")
parser.add_option("-T","--targetpopvcffile_withdepth",dest="targetpopvcffile_withdepth",action="append",nargs=2,help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-R","--refpopvcffile_withdepth",dest="refpopvcffile_withdepth",action="append",nargs=2,help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-t","--topleveltablejudgeancestral",dest="topleveltablejudgeancestral",help="R(r)/G(g)")
parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
parser.add_option("-s","--slideSize",dest="slideSize",help="default infile2_infile1")#
parser.add_option("-m","--minlength",dest="minlength")
parser.add_option("-n","--numberofindvdoftargetpop_todividintobin",dest="numberofindvdoftargetpop_todividintobin",default="o",help="conflit with correlationfile")
parser.add_option("-o","--outfileprewithpath",dest="outfileprewithpath")
parser.add_option("-c","--correlationfile",dest="correlationfile",default=None,help="conflit with numberofindvdoftargetpop_todividintobin")

parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

windowWidth=int(options.winwidth)
slideSize=int(options.slideSize)
mindepthtojudefixed=20
if __name__ == '__main__':
    if options.correlationfile==None:
        d_increase=fractions.Fraction(1, (2*int(options.numberofindvdoftargetpop_todividintobin)))
        d_increase=round(d_increase,11)
        intervalFileName=options.outfileprewithpath+".freq_correlation"
        freq_correlation_config=open(intervalFileName,"w")
        minvalue=0.000000000000
        freq_xaxisKEY_yaxisseqVALUERelation={}
        
        while minvalue+d_increase<=1:
            freq_xaxisKEY_yaxisseqVALUERelation[(minvalue,minvalue+d_increase+0.00000000004)]=[]
            
            minvalue+=d_increase


#     for a,b in sorted(freq_xaxisKEY_yaxisseqVALUERelation.keys()):
#         print(a,b)
    outputname=options.outfileprewithpath
    poplist=[];methodlist=[];listofpopvcfmapOfAChr=[];vcfnameKEY_depthobjVALUE={}
    vcfnameKEY_depthfilename_titlenameVALUE={}#{ vcftablename1:[depthfilename1,name1,name2] , vcftablename2:[depthfilename2,name1,name2] } or {vcftablename1:None, vcftablename2:None}
    vcfnamelist=[]
    N_of_targetpop=len(options.targetpopvcffile_withdepth)
    N_of_refpop=len(options.refpopvcffile_withdepth)
    dbvariantstools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.vcfdbname)
    genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    obsexpcaculator=Caculators.Caculate_S_ObsExp_difference(mindeptojudgefix,N_of_targetpop,N_of_refpop,dbvariantstools,options.topleveltablejudgeancestral)
    for vcf,depthconfig in options.targetpopvcffile_withdepth[:]+options.refpopvcffile_withdepth[:]:
        listofpopvcfmapOfAChr.append({})
        outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcf).group(0))[0])
        vcfnamelist.append(vcf)
        poplist.append(VCFutil.VCF_Data(vcf))  # new a class
        if depthconfig.lower()!="none":
            vcfnameKEY_depthfilename_titlenameVALUE[vcf]=[]
            fp=open(depthconfig,'r')
            for line in fp:
                depthfile_obj=re.search(r"depthfilename=(.*)",line.strip())
                if depthfile_obj!=None:
                    #the two for make freq_correlation config file
                    vcfnameKEY_depthfilename_titlenameVALUE[vcf].append(depthfile_obj.group(1).strip())#append depthfilename1
                    vcfnameKEY_depthobjVALUE[vcf]=Util.GATK_depthfile(vcfnameKEY_depthfilename_titlenameVALUE[vcf][0],vcfnameKEY_depthfilename_titlenameVALUE[vcf][0]+".index")
                    #the two for slide window that code exist blow
                    obsexpcaculator.species_idx_list.append([])
                    obsexpcaculator.depthobjlist.append(Util.GATK_depthfile(vcfnameKEY_depthfilename_titlenameVALUE[vcf][0],vcfnameKEY_depthfilename_titlenameVALUE[vcf][0]+".index"))                
                elif line.split():
                    titlename=line.strip()
                    idx=vcfnameKEY_depthobjVALUE[vcf].title.index("Depth_for_"+titlename)
                    vcfnameKEY_depthfilename_titlenameVALUE[vcf].append(idx)#append idxname
                    obsexpcaculator.species_idx_list[-1].append(idx)
            fp.close()
        if re.search(r"indvd[^/]+",vcf)!=None:
            obsexpcaculator.MethodToSeqpoplist.append("indvd")
            methodlist.append("indvd")
        elif re.search(r"pool[^/]+",vcf)!=None:
            obsexpcaculator.MethodToSeqpoplist.append("pool")
            methodlist.append("pool")
        else:
            print("vcfname must with 'pool' or 'indvd'")
            exit(-1)
    print("vcfnameKEY_depthfilename_titlenameVALUE",vcfnameKEY_depthfilename_titlenameVALUE,"obsexpcaculator.MethodToSeqpoplist",obsexpcaculator.MethodToSeqpoplist,"obsexpcaculator.species_idx_list",obsexpcaculator.species_idx_list,sep="\n")   
    outfile = open(outputname + ".earlypostiveselected"+str(windowWidth)+"_"+str(slideSize), 'w')

    
    
    print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile)
    print("select","select count(*) from "+Util.pekingduckchromtable + " where chrlength>="+options.minlength)
    totalChroms = genomedbtools.operateDB("select","select count(*) from "+Util.pekingduckchromtable + " where chrlength>="+options.minlength)[0][0]
    if options.correlationfile==None:
#make correlation file     
        for i in range(0,totalChroms,20):
            currentsql="select * from " + Util.pekingduckchromtable+" where chrlength>="+options.minlength+" order by "+primaryID+" limit "+str(i)+",20"
            result=genomedbtools.operateDB("select",currentsql)
             
            for row in result:
                currentchrID=row[0].strip()
                currentchrLen=int(row[1])
                for vcfobj in poplist:
                    if currentchrID in vcfobj.VcfIndexMap:
                        break
                else:
                    print("this chr doesn't exist in anypop")
                    continue
                #this chr exist in one of the vcffile
                for vcfobj_idx in range(len(poplist)):
                    listofpopvcfmapOfAChr[vcfobj_idx]={}
                    listofpopvcfmapOfAChr[vcfobj_idx][currentchrID]=poplist[vcfobj_idx].getVcfListByChrom(currentchrID)
                target_ref_SNPs=Util.alinmultPopSnpPos(listofpopvcfmapOfAChr, "o")
                for snp_aligned in target_ref_SNPs[currentchrID]:
                    if len(snp_aligned[1])!=1 or len(snp_aligned[2])!=1:
                        print("multple allele",snp_aligned)
                        continue
                    curpos=int(snp_aligned[0])
                    snp=dbvariantstools.operateDB("select","select * from "+options.topleveltablejudgeancestral+" where chrID='"+currentchrID+"' and snp_pos='"+str(curpos)+"'")
                    if not snp:
                        print(currentchrID,curpos,"snp not find in db,skip")
                        continue
                    else:#judge the ancenstrall allele
                        depthlist1=re.split(r",",snp[0][7])
                        depthlist2=re.split(r",",snp[0][9])
                        if len(depthlist1)==2 and len(depthlist2)==2 and (int(depthlist1[0]) + int(depthlist1[1])>=mindeptojudgefix or int(depthlist2[0]) + int(depthlist2[1])>=mindeptojudgefix) and ((depthlist1[0].strip()=="0" and depthlist2[0].strip()=="0") or (depthlist1[1].strip()=="0" and depthlist2[1].strip()=="0") ):
                            if depthlist1[0].strip()=="0" and depthlist2[0].strip()=="0":
                                A_base_idx=1
                            elif depthlist1[1].strip()=="0" and depthlist2[1].strip()=="0":
                                A_base_idx=0
                            else:
                                print(snp,"never get here!")
                        elif (len(depthlist1)==2 and  snp[0][9] == "no covered" and int(depthlist1[0]) + int(depthlist1[1])>=mindeptojudgefix and (depthlist1[0].strip()=="0" or depthlist1[1].strip()=="0" ))   or (snp[0][7]=="no covered" and len(depthlist2)==2 and int(depthlist2[0]) + int(depthlist2[1])>=mindeptojudgefix and (depthlist2[1].strip()=="0" or depthlist2[0].strip()=="0")):
                            if (snp[0][9] == "no covered" and depthlist1[0].strip()=="0") or (snp[0][7]=="no covered" and depthlist2[0].strip()=="0"):
                                A_base_idx=1
                            elif (snp[0][9] == "no covered" and depthlist1[1].strip()=="0") or (snp[0][7]=="no covered" and depthlist2[1].strip()=="0"):
                                A_base_idx=0
                            else:
                                print(snp,"never get here!")
                        else:
                            print(snp,"skip snp")
                            continue
                    #juge the ancenstral context
                    ancestrallcontext=snp[0][5].strip()[0].upper()+snp[0][3+A_base_idx].strip().upper()+snp[0][5].strip()[2].upper()
                    if "CG" in ancestrallcontext or "GC" in ancestrallcontext:
                        print("skip CG site",ancestrallcontext)
                        continue
                    ##########x-axis
                    countedAF=0;target_DAF_sum=0#;noofnocoveredsample=0
                    for i in range(3,N_of_targetpop+3):
                        if snp_aligned[i]==None:
                            if vcfnameKEY_depthobjVALUE==None:
                                print("no depth file")
                                continue
                            else:
                                depth_linelist=vcfnameKEY_depthobjVALUE[vcfnamelist[i-3]].getdepthByPos_optimized(currentchrID,curpos)
                                sum_depth=0
                                for idx in vcfnameKEY_depthfilename_titlenameVALUE[vcfnamelist[i-3]][1:]:
                                    sum_depth+=int(depth_linelist[idx])
                                if sum_depth>=mindepthtojudefixed:
                                    AF=0
                                else:
    #                                 noofnocoveredsample+=1
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
                    if target_DAF_sum==0 or countedAF==0:
                        print("skip this snp,because it fiexd as ancestral or no covered in this pos in target pops",snp_aligned,snp)
                        continue
                    target_DAF=target_DAF_sum/countedAF
                    ###############y-axis
                    countedAF=0;rer_DAF_sum=0
                    for i in range(3+N_of_targetpop,N_of_refpop+N_of_targetpop+3):
                        if snp_aligned[i]==None:
                            if vcfnameKEY_depthobjVALUE==None:
                                continue
                            else:
                                depth_linelist=vcfnameKEY_depthobjVALUE[vcfnamelist[i-3-N_of_targetpop]].getdepthByPos_optimized(currentchrID,curpos)
                                sum_depth=0
                                for idx in vcfnameKEY_depthfilename_titlenameVALUE[vcfnamelist[i-3-N_of_targetpop]][1:]:
                                    sum_depth+=int(depth_linelist[idx])
                                if sum_depth>=mindepthtojudefixed:
                                    AF=0
                                else:
                                    continue
                        else:
                            if methodlist[i-3-N_of_targetpop]=="indvd":
                                AF = float(re.search(r"AF=([\d\.]+);", snp_aligned[i][0]).group(1))
                            elif methodlist[i-3-N_of_targetpop]=="pool":
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
                                    continue
                                AF=altalleledep/(altalleledep+refdep)
                        if A_base_idx==0:
                            DAF=1-AF
                        elif A_base_idx==1:
                            DAF=AF
                        rer_DAF_sum+=DAF;countedAF+=1
                    if countedAF==0:
                        print("skip this snp,because it  no covered in this pos in target pops",snp_aligned,snp)
                        continue
                    ######collect according bins
                    for a,b in sorted(freq_xaxisKEY_yaxisseqVALUERelation.keys()):
                        if target_DAF >a and target_DAF<=b:
                            freq_xaxisKEY_yaxisseqVALUERelation[(a,b)].append(rer_DAF_sum/countedAF)
                            break
        freq_xaxisKEY_yaxisVALUERelation={}
        for a,b in sorted(freq_xaxisKEY_yaxisseqVALUERelation.keys()):
            freq_xaxisKEY_yaxisVALUERelation[(a,b)]=numpy.mean(freq_xaxisKEY_yaxisseqVALUERelation[(a,b)])
            print('%.12f'%a,'%.12f'%(b),'%.12f'%(freq_xaxisKEY_yaxisVALUERelation[(a,b)]),sep="\t",file=freq_correlation_config)
        freq_correlation_config.close()
    else:
        correlationfile=open(options.correlationfile,'r')
        freq_xaxisKEY_yaxisVALUERelation={}
        for line in correlationfile:
            linelist=re.split(r"\s+",line.strip())
            freq_xaxisKEY_yaxisVALUERelation[float(linelist[0]),float(linelist[1])]=float(linelist[2])
        correlationfile.close()
#     for a,b in sorted(freq_xaxisKEY_yaxisVALUERelation.keys()):
#         print('%.12f'%a,'%.12f'%(b),'%.12f'%(freq_xaxisKEY_yaxisVALUERelation[(a,b)]),sep="\t")
    #slide window to caculate S
#     exit()
    obsexpcaculator.freq_xaxisKEY_yaxisVALUERelation=freq_xaxisKEY_yaxisVALUERelation
    win = Util.Window()
    obsexpsignalmapbychrom={}

    for i in range(0,totalChroms,20):
        currentsql="select * from " + Util.pekingduckchromtable+" where chrlength>="+options.minlength+" order by "+primaryID+" limit "+str(i)+",20"
        result=genomedbtools.operateDB("select",currentsql)
        for row in result:
            currentchrID=row[0].strip()
            currentchrLen=int(row[1])
            for vcfobj in poplist:
                if currentchrID in vcfobj.VcfIndexMap:
                    break
            else:
                print("this chr doesn't exist in anypop")
                fillNA=[(0,0,0,'NA')]
                for i in range(int(currentchrLen/slideSize)):
                    fillNA.append((0,0,0,'NA'))
                obsexpsignalmapbychrom[currentchrID]=fillNA
                continue
            #this chr exist in one of the vcffile,then alinmultPopSnpPos
            for vcfobj_idx in range(len(poplist)):
                listofpopvcfmapOfAChr[vcfobj_idx]={}
                listofpopvcfmapOfAChr[vcfobj_idx][currentchrID]=poplist[vcfobj_idx].getVcfListByChrom(currentchrID)
            target_ref_SNPs=Util.alinmultPopSnpPos(listofpopvcfmapOfAChr, "o")
            obsexpcaculator.currentchrID=currentchrID
            ##########
            win.slidWindowOverlap(target_ref_SNPs[currentchrID], currentchrLen, windowWidth, slideSize, obsexpcaculator)
            obsexpsignalmapbychrom[currentchrID]=copy.deepcopy(win.winValueL)
    for i in range(0,totalChroms,20):
        currentsql="select * from " + Util.pekingduckchromtable+" where chrlength>="+options.minlength+" order by "+primaryID+" limit "+str(i)+",20"
        result=genomedbtools.operateDB("select",currentsql)
        for row in result:
            currentchrID=row[0].strip()
            if currentchrID in obsexpsignalmapbychrom:
                for i in range(len(obsexpsignalmapbychrom[currentchrID])):
                    if obsexpsignalmapbychrom[currentchrID][i][3]=="NA":
                        print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][0]) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][1]) + "\t"+str(obsexpsignalmapbychrom[currentchrID][i][2])+"\t" + "NA" + "\t" + "NA", file=outfile)
                    else:
                        print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][0]) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][1]) + "\t" +str(obsexpsignalmapbychrom[currentchrID][i][2])+"\t"+ '%.15f'%(obsexpsignalmapbychrom[currentchrID][i][3][0]) + "\t" + '%.12f'%(obsexpsignalmapbychrom[currentchrID][i][3][1]), file=outfile)
    outfile.close()
    dbvariantstools.disconnect()
    genomedbtools.disconnect()
#             for vcfobj in poplist:
#                 if currentchrID in vcfobj.VcfIndexMap:
#                     break
#             else:
#                 continue
#             #this chr exist in one of the vcffile
#             for vcfobj_idx in range(len(poplist)):
#                 listofpopvcfmapOfAChr[vcfobj_idx]={}
#                 listofpopvcfmapOfAChr[vcfobj_idx][currentchrID]=poplist[vcfobj_idx].getVcfListByChrom(currentchrID)
#                 if vcfnameKEY_depthobjVALUE!={}:
#                     vcfnameKEY_depthobjVALUE[(options.targetpopvcffile_withdepth[:]+options.refpopvcffile_withdepth[:])[vcfobj_idx][0]].depthfilefp.seek(vcfnameKEY_depthobjVALUE[(options.targetpopvcffile_withdepth[:]+options.refpopvcffile_withdepth[:])[vcfobj_idx][0]].covfileidx[currentchrID])
#             target_ref_SNPs=Util      

    
    
    