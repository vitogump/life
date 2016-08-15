'''
Created on 2015-8-21

@author: liurui
'''
from optparse import OptionParser
import re, numpy, fractions, copy, os, pysam
from src.NGS.Service import Ancestralallele
from NGS.BasicUtil import *
import src.NGS.BasicUtil.DBManager as dbm


parser = OptionParser()
parser.add_option("-c","--chromlistfilename",dest="chromlistfilename",help="early,pairfst,pbs,lsbl,is")
parser.add_option("-p","--typeOfcalculate",dest="typeOfcalculate")
parser.add_option("-t","--topleveltablejudgeancestral",dest="topleveltablejudgeancestral",help="assigned only if -p early")
parser.add_option("-T","--targetpopvcfconfig",dest="targetpopvcfconfig",action="append",help="firstline is vcffilename=,the rest lines can be none or bamfilename per line")
parser.add_option("-R","--refpopvcffileconfig",dest="refpopvcffileconfig",action="append",help="firstline is vcffilename=,the rest lines can be none or bamfilename per line")
parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
parser.add_option("-s","--slideSize",dest="slideSize",help="default infile2_infile1")#
parser.add_option("-C","--correlationfile",dest="correlationfile",default=None,help="conflit with numberofindvdoftargetpop_todividintobin")
parser.add_option("-o","--outfileprewithpath",dest="outfileprewithpath")
parser.add_option("-m","--masterpid",dest="masterpid")
# parser.add_option("-f","--reffa",dest="reffa",help="used for fill toplevel context")
(options, args) = parser.parse_args()
mindeptojudgefix=20
windowWidth=int(options.winwidth)
slideSize=int(options.slideSize)
if __name__ == '__main__':
    print("runSlave_slidewin process ID",os.getpid(),"start")
    flankseqfafilename=options.chromlistfilename+str(os.getpid())+"snpflankseq.fa"
    outputname=options.outfileprewithpath
    listofpopvcfmapOfAChr=[];vcfnamelist=[]
    N_of_targetpop=len(options.targetpopvcfconfig)
    N_of_refpop=len(options.refpopvcffileconfig)
    
#     genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    if options.typeOfcalculate=="early":
        dbvariantstools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.vcfdbname)
        dynamicIU_toptable_obj=Ancestralallele.dynamicInsertUpdateAncestralContext(dbvariantstools,Util.beijingreffa,options.topleveltablejudgeancestral)
        
        obsexpcaculator=Caculators.Caculate_S_ObsExp_difference(mindeptojudgefix,N_of_targetpop,N_of_refpop,dbvariantstools,options.topleveltablejudgeancestral)
        obsexpcaculator.dynamicIU_toptable_obj=dynamicIU_toptable_obj
        obsexpcaculator.flankseqfafile=open(flankseqfafilename,"w")
    elif options.typeOfcalculate=="pairfst":
        obsexpcaculator=Caculators.Caculate_pairFst(mindeptojudgefix,N_of_targetpop,N_of_refpop)
    elif options.typeOfcalculate=="is":
        obsexpcaculator
    for vcfconfigfilename in options.targetpopvcfconfig[:]+options.refpopvcffileconfig[:]:
        listofpopvcfmapOfAChr.append({})
        vcfconfig=open(vcfconfigfilename,"r")
        for line in vcfconfig:
            vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
            if vcffilename_obj!=None:
                vcfname=vcffilename_obj.group(1).strip()
                vcfnamelist.append(vcfname)
                outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcfname).group(0))[0])[:3]
                obsexpcaculator.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                obsexpcaculator.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
            elif line.split():
                obsexpcaculator.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
        vcfconfig.close()
        if re.search(r"indvd[^/]+",vcfname)!=None:
            obsexpcaculator.MethodToSeqpoplist.append("indvd")

        elif re.search(r"pool[^/]+",vcfname)!=None:
            obsexpcaculator.MethodToSeqpoplist.append("pool")

        else:
            print("vcfname must with 'pool' or 'indvd'")
            exit(-1)   
    chrlistfilewithoutpath=re.search(r"[^/]*$",options.chromlistfilename).group(0)
    plainname=re.search(r"[^/]*$",outputname).group(0)
    if len(plainname)>=250:
        outputname=outputname[:-(len(plainname)-250)]
    outfile = open(outputname + ".earlypostiveselected"+str(windowWidth)+"_"+str(slideSize)+chrlistfilewithoutpath, 'w')

    aaaa=open(options.outfileprewithpath+".slidwin_filelist"+options.masterpid,'a')
    print(outputname + ".earlypostiveselected"+str(windowWidth)+"_"+str(slideSize)+chrlistfilewithoutpath,file=aaaa)
    aaaa.close()
    print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile)
    freq_correlation_configFileName=options.outfileprewithpath+".freq_correlation_merged"
    if options.correlationfile!=freq_correlation_configFileName:
        print("warning !",freq_correlation_configFileName," is not equal to ",options.correlationfile)
    freq_correlation_config=open(options.correlationfile,"r")
    final_freq_xaxisKEY_yaxisVALUERelation={}
    for line in freq_correlation_config:
        if line.split():
            linelist=re.split(r"\t",line.strip())
            a=float(linelist[0]);b=float(linelist[1]);yaxisfreq=float(linelist[2])
            final_freq_xaxisKEY_yaxisVALUERelation[(a,b)]=yaxisfreq
    obsexpcaculator.freq_xaxisKEY_yaxisVALUERelation=final_freq_xaxisKEY_yaxisVALUERelation
    freq_correlation_config.close()
    win = Util.Window()
    obsexpsignalmapbychrom={}
    chromlistfile=open(options.chromlistfilename,"r")
    chromlist=[]
    for chrrow in chromlistfile:
        chrrowlist=re.split(r'\s+',chrrow.strip())
        chromlist.append((chrrowlist[0].strip(),int(chrrowlist[1].strip())))


    for currentchrID,currentchrLen in chromlist:
        for vcfname in vcfnamelist:
            vcfobj=obsexpcaculator.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][0]
#         for vcfobj in poplist:
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
#         for vcfobj_idx in range(len(poplist)):
        for vcfobj_idx in range(len(vcfnamelist)):
            listofpopvcfmapOfAChr[vcfobj_idx]={}
            vcfobj=obsexpcaculator.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][0]
            print(vcfnamelist[vcfobj_idx],"getvcf")
            listofpopvcfmapOfAChr[vcfobj_idx][currentchrID]=vcfobj.getVcfListByChrom(currentchrID)
        target_ref_SNPs=Util.alinmultPopSnpPos(listofpopvcfmapOfAChr, "o")
        obsexpcaculator.currentchrID=currentchrID
        obsexpcaculator.alignedSNP_absentinfo={}
        obsexpcaculator.alignedSNP_absentinfo[currentchrID]=[]
        ##########
        win.slidWindowOverlap(target_ref_SNPs[currentchrID], currentchrLen, windowWidth, slideSize, obsexpcaculator)
        obsexpsignalmapbychrom[currentchrID]=copy.deepcopy(win.winValueL)

    for currentchrID,currentchrLen in chromlist:

        if currentchrID in obsexpsignalmapbychrom:
            for i in range(len(obsexpsignalmapbychrom[currentchrID])):
                if obsexpsignalmapbychrom[currentchrID][i][3]=="NA":
                    print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][0]) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][1]) + "\t"+str(obsexpsignalmapbychrom[currentchrID][i][2])+"\t" + "NA" + "\t" + "NA", file=outfile)
                else:
#                     zS=(obsexpsignalmapbychrom[currentchrID][i][3][0]-exception)/std1
                    print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][0]) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][1]) + "\t" +str(obsexpsignalmapbychrom[currentchrID][i][2])+"\t"+ '%.15f'%(obsexpsignalmapbychrom[currentchrID][i][3][0]) + "\t" + '%.12f'%(obsexpsignalmapbychrom[currentchrID][i][3][1]), file=outfile)
    outfile.close()
    chromlistfile.close()
    dbvariantstools.disconnect()
#     genomedbtools.disconnect()
    print("runSlave_slidewin process ID",os.getpid(),"done")