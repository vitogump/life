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
mindeptojudgefix=20#for pool only
windowWidth=int(options.winwidth)
slideSize=int(options.slideSize)
if __name__ == '__main__':
    print("runSlave_slidewin process ID",os.getpid(),"start")
    flankseqfafilename=options.chromlistfilename+str(os.getpid())+"snpflankseq.fa"
    
#     genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    chrlistfilewithoutpath=re.search(r"[^/]*$",options.chromlistfilename).group(0)
    if options.typeOfcalculate=="early":
        dbvariantstools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.vcfdbname)
        dynamicIU_toptable_obj=Ancestralallele.dynamicInsertUpdateAncestralContext(dbvariantstools,Util.beijingreffa,options.topleveltablejudgeancestral)
        
        obsexpcaculator=Caculators.Caculate_S_ObsExp_difference(mindeptojudgefix,options.targetpopvcfconfig,options.refpopvcffileconfig,dbvariantstools,options.topleveltablejudgeancestral,options.outfileprewithpath)
        obsexpcaculator.dynamicIU_toptable_obj=dynamicIU_toptable_obj
        obsexpcaculator.flankseqfafile=open(flankseqfafilename,"w")
        plainname=re.search(r"[^/]*$",obsexpcaculator.outputname).group(0)
        if len(plainname)>=250:
            obsexpcaculator.outputname=obsexpcaculator.outputname[:-(len(plainname)-250)]
        outputname=obsexpcaculator.outputname
        outfile = open( outputname+ ".earlypostiveselected"+str(windowWidth)+"_"+str(slideSize)+chrlistfilewithoutpath, 'w')
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
    elif options.typeOfcalculate=="pairfst":
        obsexpcaculator=Caculators.Caculate_pairFst(mindeptojudgefix,options.targetpopvcfconfig,options.refpopvcffileconfig)
        plainname=re.search(r"[^/]*$",obsexpcaculator.outputname).group(0)
        if len(plainname)>=250:
            obsexpcaculator.outputname=obsexpcaculator.outputname[:-(len(plainname)-250)]
        outfile = open(obsexpcaculator.outputname + ".earlypostiveselected"+str(windowWidth)+"_"+str(slideSize)+chrlistfilewithoutpath, 'w')
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile)

    elif options.typeOfcalculate=="is":
        obsexpcaculator=Caculators.Caculate_IS(mindeptojudgefix/2,options.targetpopvcfconfig,options.refpopvcffileconfig)
        plainname=re.search(r"[^/]*$",options.outfileprewithpath).group(0)
        if len(plainname)>=250:
            outputname=options.outfileprewithpath[:-(len(plainname)-250)]
        outputname=options.outfileprewithpath
        outfile = open(outputname + "."+options.typeOfcalculate+str(windowWidth)+"_"+str(slideSize)+chrlistfilewithoutpath, 'w')
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos",*obsexpcaculator.vcfname_combination,sep="\t",file=outfile)

    aaaa=open(options.outfileprewithpath+".slidwin_filelist"+options.masterpid,'a')
    print(outputname + "."+options.typeOfcalculate+str(windowWidth)+"_"+str(slideSize)+chrlistfilewithoutpath,file=aaaa)
    aaaa.close()


    win = Util.Window()
    obsexpsignalmapbychrom={}
    chromlistfile=open(options.chromlistfilename,"r")
    chromlist=[]
    for chrrow in chromlistfile:
        chrrowlist=re.split(r'\s+',chrrow.strip())
        chromlist.append((chrrowlist[0].strip(),int(chrrowlist[1].strip())))


    for currentchrID,currentchrLen in chromlist:
        for vcfname in obsexpcaculator.vcfnamelist:
            vcfobj=obsexpcaculator.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][0]
#         for vcfobj in poplist:
            if currentchrID in vcfobj.VcfIndexMap:
                break
        else:
            print("this chr doesn't exist in anypop")
            fillNA=obsexpcaculator.getResult()#[(0,0,0,'NA')]
            for i in range(int(currentchrLen/slideSize)):
                fillNA.append(obsexpcaculator.getResult())#(0,0,0,'NA')
            obsexpsignalmapbychrom[currentchrID]=fillNA
            continue
        #this chr exist in one of the vcffile,then alinmultPopSnpPos
#         for vcfobj_idx in range(len(poplist)):
        for vcfobj_idx in range(len(obsexpcaculator.vcfnamelist)):
            obsexpcaculator.listOfpopvcfRecsmapByAChr[vcfobj_idx]={}
            vcfobj=obsexpcaculator.vcfnameKEY_vcfobj_pyBAMfilesVALUE[obsexpcaculator.vcfnamelist[vcfobj_idx]][0]
            print(obsexpcaculator.vcfnamelist[vcfobj_idx],"getvcf")
            obsexpcaculator.listOfpopvcfRecsmapByAChr[vcfobj_idx][currentchrID]=vcfobj.getVcfListByChrom(currentchrID)
        target_ref_SNPs=Util.alinmultPopSnpPos(obsexpcaculator.listOfpopvcfRecsmapByAChr, "o")
        if options.typeOfcalculate=="early":
            obsexpcaculator.currentchrID=currentchrID
            obsexpcaculator.dynamicIU_toptable_obj.currentchrLen=currentchrLen
            obsexpcaculator.alignedSNP_absentinfo={}
            obsexpcaculator.alignedSNP_absentinfo[currentchrID]=[]
        ##########
        win.slidWindowOverlap(target_ref_SNPs[currentchrID], currentchrLen, windowWidth, slideSize, obsexpcaculator)
        obsexpsignalmapbychrom[currentchrID]=copy.deepcopy(win.winValueL)
    
    for currentchrID,currentchrLen in chromlist:

        if currentchrID in obsexpsignalmapbychrom:
            for i in range(len(obsexpsignalmapbychrom[currentchrID])):
                if options.typeOfcalculate=="early":
                    if obsexpsignalmapbychrom[currentchrID][i][3]=="NA":
                        print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][0]) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][1]) + "\t"+str(obsexpsignalmapbychrom[currentchrID][i][2])+"\t" + "NA" + "\t" + "NA", file=outfile)
                    else:
    #                     zS=(obsexpsignalmapbychrom[currentchrID][i][3][0]-exception)/std1
                        print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][0]) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][1]) + "\t" +str(obsexpsignalmapbychrom[currentchrID][i][2])+"\t"+ '%.15f'%(obsexpsignalmapbychrom[currentchrID][i][3][0]) + "\t" + '%.12f'%(obsexpsignalmapbychrom[currentchrID][i][3][1]), file=outfile)
                elif options.typeOfcalculate=="pairfst":
                    pass
                elif options.typeOfcalculate=="is":
                    print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][0])+ "\t" + str(obsexpsignalmapbychrom[currentchrID][i][1]),*obsexpsignalmapbychrom[currentchrID][i][3],sep="\t", file=outfile)# + "\t"+str()
    outfile.close()
    chromlistfile.close()
    if options.topleveltablejudgeancestral!=None:
        dbvariantstools.disconnect()
#     genomedbtools.disconnect()
    print("runSlave_slidewin process ID",os.getpid(),"done")