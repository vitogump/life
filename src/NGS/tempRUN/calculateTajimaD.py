'''
Created on 2016-1-10

@author: liurui
'''
from optparse import OptionParser
import re,os


parser = OptionParser()
parser.add_option("-d","--tajimaDdir",dest="tajimaDdir",help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
# parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
parser.add_option("-p","--prestr",dest="prestr",default="tajimaDandwattersons_theta",help="the files name were consist of prestr + chrNO")
parser.add_option("-o","--outfilename",dest="outfilename",help="default infile2_infile1")#

(options, args) = parser.parse_args()
# windowWidth=int(options.winwidth)
# slideSize=int(options.slideSize)
pathtovcftools="/pub/tool/vcftools_0.1.12b/bin/vcftools "
if __name__ == '__main__':
    
    datafiles = os.listdir(path=options.tajimaDdir)
    outfile=open(options.outfilename,"w")
    for datafile in datafiles:
        if datafile.startswith(options.prestr):
            i=0
            chrID=re.search(r""+options.prestr.strip()+"([\w\W]+)",datafile.strip()).group(1)
            f=open(options.tajimaDdir+"/"+datafile,'r')
            f.readline()
            for line in f:
                linelist=re.split(r",",line)
                print(chrID,str(i),str(i*20000),str(i*20000+40000),linelist[2],linelist[1],linelist[11],sep="\t",file=outfile)#1:"tajima's D" 11:"theta_Watterson"
                i+=1
            f.close()
    outfile.close()
                
#     vcffile=open(options.vcffile,"r")
#     for titleline in vcffile:
#         if re.search(r'^#',titleline)!=None:
#             break
#         if re.search(r"^##contig=<ID=([\w\W]+),length=[\d]+>",titleline)!=None:
#             curchr=re.search(r"^##contig=<ID=([\w\W]+),length=[\d]+>",titleline).group(1)