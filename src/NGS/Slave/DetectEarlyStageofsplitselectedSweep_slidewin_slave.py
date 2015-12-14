'''
Created on 2015-8-21

@author: liurui
'''
from optparse import OptionParser
from NGS.BasicUtil import *
import src.NGS.BasicUtil.DBManager as dbm
import re, numpy, fractions, copy,os


parser = OptionParser()
parser.add_option("-c","--chromlistfilename",dest="chromlistfilename")
parser.add_option("-t","--topleveltablejudgeancestral",dest="topleveltablejudgeancestral",help="R(r)/G(g)")
parser.add_option("-T","--targetpopvcffile_withdepth",dest="targetpopvcffile_withdepth",action="append",nargs=2,help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-R","--refpopvcffile_withdepth",dest="refpopvcffile_withdepth",action="append",nargs=2,help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
parser.add_option("-s","--slideSize",dest="slideSize",help="default infile2_infile1")#
parser.add_option("-C","--correlationfile",dest="correlationfile",default=None,help="conflit with numberofindvdoftargetpop_todividintobin")
parser.add_option("-o","--outfileprewithpath",dest="outfileprewithpath")
(options, args) = parser.parse_args()
mindeptojudgefix=20
windowWidth=int(options.winwidth)
slideSize=int(options.slideSize)
if __name__ == '__main__':
    print("runSlave_slidewin process ID",os.getpid(),"start")
    outputname=options.outfileprewithpath
    poplist=[];listofpopvcfmapOfAChr=[]
    N_of_targetpop=len(options.targetpopvcffile_withdepth)
    N_of_refpop=len(options.refpopvcffile_withdepth)
    dbvariantstools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.vcfdbname)
    genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    obsexpcaculator=Caculators.Caculate_S_ObsExp_difference(mindeptojudgefix,N_of_targetpop,N_of_refpop,dbvariantstools,options.topleveltablejudgeancestral)
    for vcf,depthconfig in options.targetpopvcffile_withdepth[:]+options.refpopvcffile_withdepth[:]:
        listofpopvcfmapOfAChr.append({})
        outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcf).group(0))[0])[:2]
        poplist.append(VCFutil.VCF_Data(vcf))  # new a class
        if depthconfig.lower()!="none":

            fp=open(depthconfig,'r')
            for line in fp:
                depthfilename_obj=re.search(r"depthfilename=(.*)",line.strip())
                if depthfilename_obj!=None:
                    #the two for make freq_correlation config file

                    #the two for slide window that code exist blow
                    obsexpcaculator.species_idx_list.append([])
                    obsexpcaculator.depthobjlist.append(Util.GATK_depthfile(depthfilename_obj.group(1).strip(),depthfilename_obj.group(1).strip()+".index"))                                
                elif line.split():
                    titlename=line.strip()
                    idx=obsexpcaculator.depthobjlist[-1].title.index("Depth_for_"+titlename)

                    obsexpcaculator.species_idx_list[-1].append(idx)

            fp.close()
        if re.search(r"indvd[^/]+",vcf)!=None:
            obsexpcaculator.MethodToSeqpoplist.append("indvd")

        elif re.search(r"pool[^/]+",vcf)!=None:
            obsexpcaculator.MethodToSeqpoplist.append("pool")

        else:
            print("vcfname must with 'pool' or 'indvd'")
            exit(-1)   
    chrlistfilewithoutpath=re.search(r"[^/]*$",options.chromlistfilename).group(0)
    plainname=re.search(r"[^/]*$",outputname).group(0)
    if len(plainname)>=250:
        outputname=outputname[:-(len(plainname)-250)]
    outfile = open(outputname + ".earlypostiveselected"+str(windowWidth)+"_"+str(slideSize)+chrlistfilewithoutpath, 'w')

    
    
    print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile)
    freq_correlation_configFileName=options.outfileprewithpath+".freq_correlation_merged"
    if options.correlationfile!=freq_correlation_configFileName:
        print("what's wrong??",freq_correlation_configFileName,options.correlationfile)
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
#     for i in range(0,totalChroms,20):
    winCrossGenome=[]
    for chrom in obsexpsignalmapbychrom.keys():
        for i in range(len(obsexpsignalmapbychrom[chrom])):
            if obsexpsignalmapbychrom[chrom][i][3]!="NA":
                winCrossGenome.append(obsexpsignalmapbychrom[chrom][i][3][0])
    exception =numpy.mean(winCrossGenome)
    std0=numpy.std(winCrossGenome,ddof=0)
    std1=numpy.std(winCrossGenome,ddof=1)
    for currentchrID,currentchrLen in chromlist:
#         currentsql="select * from " + Util.pekingduckchromtable+" where chrlength>="+options.minlength+" order by "+primaryID+" limit "+str(i)+",20"
#         result=genomedbtools.operateDB("select",currentsql)
#         for row in result:
#             currentchrID=row[0].strip()
        if currentchrID in obsexpsignalmapbychrom:
            for i in range(len(obsexpsignalmapbychrom[currentchrID])):
                if obsexpsignalmapbychrom[currentchrID][i][3]=="NA":
                    print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][0]) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][1]) + "\t"+str(obsexpsignalmapbychrom[currentchrID][i][2])+"\t" + "NA" + "\t" + "NA", file=outfile)
                else:
                    zS=(obsexpsignalmapbychrom[currentchrID][i][3][0]-exception)/std1
                    print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][0]) + "\t" + str(obsexpsignalmapbychrom[currentchrID][i][1]) + "\t" +str(obsexpsignalmapbychrom[currentchrID][i][2])+"\t"+ '%.15f'%(obsexpsignalmapbychrom[currentchrID][i][3][0]) + "\t" + '%.12f'%(zS), file=outfile)
    outfile.close()
    chromlistfile.close()
    dbvariantstools.disconnect()
    genomedbtools.disconnect()
    print("runSlave_slidewin process ID",os.getpid(),"done")