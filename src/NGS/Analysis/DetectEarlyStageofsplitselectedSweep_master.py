'''
Created on 2015-8-1

@author: liurui
'''
from multiprocessing.dummy import Pool
from optparse import OptionParser
import os
import re, numpy, fractions, copy

from NGS.BasicUtil import *
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
parser.add_option("-c","--chromlistfilename",dest="chromlistfilename",action="append")
parser.add_option("-n","--numberofindvdoftargetpop_todividintobin",dest="numberofindvdoftargetpop_todividintobin",default="o",help="conflit with correlationfile")
parser.add_option("-o","--outfileprewithpath",dest="outfileprewithpath")
parser.add_option("-C","--correlationfile",dest="correlationfile",default=None,help="conflit with numberofindvdoftargetpop_todividintobin")
parser.add_option("-m","--numberofthreads",dest="numberofthreads")
parser.add_option("-1","--pathtoslave_config",dest="pathtoslave_config",default=None)
parser.add_option("-2","--pathtoslave_slidewin",dest="pathtoslave_slidewin",default=None)
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

windowWidth=int(options.winwidth)
slideSize=int(options.slideSize)

pathtoPython="/home/bioinfo/liurui/software/Python-3.3.3/python "
def runSlave_makecorrelationfile(a):
    chromlistfilename=a[0];topleveltablename=a[1];targetpopvcffile_withdepthconfig=a[2];refpopvcffile_withdepthconfig=a[3];numberofindvdoftargetpop_todividintobin=a[4];outfileprewithpath=a[5]
    command=pathtoPython+options.pathtoslave_config+" -c "+chromlistfilename+" -t "+topleveltablename
    for vcf,depthconfig in targetpopvcffile_withdepthconfig[:]:
        command+=(" -T "+vcf+" "+depthconfig)
    for vcf,depthconfig in refpopvcffile_withdepthconfig[:]:
        command+=(" -R "+vcf+" "+depthconfig)
    a=os.system(command+" -n "+numberofindvdoftargetpop_todividintobin+" -o "+outfileprewithpath)
def runSlave_slidewin(a):
    chromlistfilename=a[0];topleveltablename=a[1];targetpopvcffile_withdepthconfig=a[2];refpopvcffile_withdepthconfig=a[3];winwidth=a[4];slideSize=a[5];correlationfile=a[6];outfileprewithpath=a[7]
    command=pathtoPython+options.pathtoslave_slidewin+" -c "+chromlistfilename+" -t "+topleveltablename
    for vcf,depthconfig in targetpopvcffile_withdepthconfig[:]:
        command+=(" -T "+vcf+" "+depthconfig)
    for vcf,depthconfig in refpopvcffile_withdepthconfig[:]:
        command+=(" -R "+vcf+" "+depthconfig)
    a=os.system(command+" -w "+winwidth+" -s "+slideSize+" -o "+outfileprewithpath+" -C "+correlationfile)
if __name__ == '__main__':
    if options.correlationfile==None:
        d_increase=fractions.Fraction(1, (2*int(options.numberofindvdoftargetpop_todividintobin)))
        d_increase=round(d_increase,11)
        minvalue=0.000000000000
        final_freq_xaxisKEY_yaxisVALUE_seq_list={}
        
        for i in range(int(options.numberofindvdoftargetpop_todividintobin)*2-1):
            final_freq_xaxisKEY_yaxisVALUE_seq_list[(minvalue,minvalue+d_increase+0.00000000004)]=[]
            minvalue+=d_increase
        else:
            final_freq_xaxisKEY_yaxisVALUE_seq_list[(minvalue,1)]=[]

        print(final_freq_xaxisKEY_yaxisVALUE_seq_list)
#     for a,b in sorted(freq_xaxisKEY_yaxisseqVALUERelation.keys()):
#         print(a,b)


    
    if options.correlationfile==None:
        if int(options.numberofthreads)!=len(options.chromlistfilename):
            print("int(options.numberofthreads)!=len(options.chromlistfilename)")
            exit(-1)
        freq_xaxisKEY_yaxisVALUERelation_maplist=[]
        pool=Pool(int(options.numberofthreads))
        parameterstuples_list=[]
        for chromlistfile in options.chromlistfilename:
            parameterstuples_list.append((chromlistfile,options.topleveltablejudgeancestral,options.targetpopvcffile_withdepth,options.refpopvcffile_withdepth,options.numberofindvdoftargetpop_todividintobin,options.outfileprewithpath))
            print(len(parameterstuples_list[-1]),parameterstuples_list[-1])
        print(len(parameterstuples_list),parameterstuples_list)
#         exit()
        pool.map(runSlave_makecorrelationfile,parameterstuples_list)
        pool.close()
        pool.join()
        f=open(options.outfileprewithpath+".freqcorrelationfilenamelist",'r')
        for freqseq_cor_filename in f:# for every file
            freqseqmap={}
            if freqseq_cor_filename.split():
                freqseq_cor_file=open(freqseq_cor_filename,'r')
                for line in freqseq_cor_file:#for every freq seq bin
                    if line.split():
                        linelist=re.split(r"\t",line.strip)
                        a=float(linelist[0]);b=float(linelist[1])
                        freqseqmap[(a,b)]=[]
                        for freq in linelist[2:]:
                            freqseqmap[(a,b)].append(float(freq))
                freq_xaxisKEY_yaxisVALUERelation_maplist.append(copy.deepcopy(freqseqmap))
        for freq_xaxisKEY_yaxisseqVALUERelation_part in freq_xaxisKEY_yaxisVALUERelation_maplist:
            for xaxis in sorted(final_freq_xaxisKEY_yaxisVALUE_seq_list.keys()):
                final_freq_xaxisKEY_yaxisVALUE_seq_list[xaxis]+=freq_xaxisKEY_yaxisseqVALUERelation_part[xaxis]
        final_freq_xaxisKEY_yaxisVALUERelation={}
        freq_correlation_configFileName=options.outfileprewithpath+".freq_correlation_merged"
        freq_correlation_config=open(freq_correlation_configFileName,"w")
        for a,b in sorted(final_freq_xaxisKEY_yaxisVALUE_seq_list.keys()):
            final_freq_xaxisKEY_yaxisVALUERelation[(a,b)]=numpy.mean(final_freq_xaxisKEY_yaxisVALUE_seq_list[(a,b)])
            print('%.12f'%a,'%.12f'%b,'%.12f'%(final_freq_xaxisKEY_yaxisVALUERelation[(a,b)]),sep="\t",file=freq_correlation_config)
        freq_correlation_config.close()
        print("freq_correlation_config is produced")
        exit()
        
    else:
        if len(options.chromlistfilename)!=1:
            print("need only one chromlistfilename")
            exit(-1)
        correlationfile=open(options.correlationfile,'r')
        final_freq_xaxisKEY_yaxisVALUERelation={}
        for line in correlationfile:
            linelist=re.split(r"\s+",line.strip())
            final_freq_xaxisKEY_yaxisVALUERelation[float(linelist[0]),float(linelist[1])]=float(linelist[2])
        correlationfile.close()
#     for a,b in sorted(final_freq_xaxisKEY_yaxisVALUE_seq_list.keys()):
#         print('%.12f'%a,'%.12f'%(b),'%.12f'%(final_freq_xaxisKEY_yaxisVALUE_seq_list[(a,b)]),sep="\t")
    #slide window to caculate S
    print("all final_freq_xaxisKEY_yaxisVALUERelation done ,slide window now")
    pool=Pool(int(options.numberofthreads))
    parameterstuples_list=[]
    for chromlistfile in options.chromlistfilename:
        parameterstuples_list.append((chromlistfile,options.topleveltablejudgeancestral,options.targetpopvcffile_withdepth,options.refpopvcffile_withdepth,options.winwidth,options.slideSize,freq_correlation_configFileName,options.outfileprewithpath))
        print(len(parameterstuples_list[-1]),parameterstuples_list[-1])
    print(len(parameterstuples_list),parameterstuples_list)
#         exit()
    pool.map(runSlave_slidewin,parameterstuples_list)
    pool.close()
    pool.join()
    print("finished")

    

    
    
    