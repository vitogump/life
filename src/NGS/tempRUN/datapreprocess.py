'''
Created on 2019年1月8日
pick up each data unit into one folder
@author: RuiLiu
'''
from optparse import OptionParser
import os,re
# from NGS.Analysis.IntervalValueToBoxplot import path

parser = OptionParser()
parser.add_option("-d", "--rootdir", dest="rootdir",help="data files directly under this dir ")
parser.add_option("-r", "--rmoriginalfile", dest="rmoriginalfile",action="store_true",default=False)
parser.add_option("-s", "--suffix", dest="suffix",action="append",default=[],help='suffix for a folder to package,["1.fq.gz","2.fq.gz"]')
parser.add_option("-o", "--outdir", dest="outdir",default=None,help="outpath,default is the same fold as -d assigned ")
(options, args) = parser.parse_args()

# cf=open(options.rootdir,'r')

filepath=options.rootdir.strip()

files=os.listdir(filepath)
print(filepath,files,options.suffix)
folderMap={}#{"AS10_":[AS10_1.fq.gz,AS10_1.fq.gz]}
# os.system("cp "+os.path.join(filepath,"md5*.txt")+" "+os.path.join(filepath,"md5*.txt.copy"))
if options.outdir==None:
    outputdir=filepath
else:
    outputdir=options.outdir.strip()
if __name__ == '__main__':
    for fi in files:
        
        fi_d=os.path.join(filepath,fi)
        if os.path.isdir(fi_d):
            pass
        else:
            for s in options.suffix:
                if re.search(r"."+s,fi)!=None:
                    l_pre=len(fi)-len(s)
                    data_d=os.path.join(outputdir,fi[:l_pre]).strip("_").strip("-")
                    if not os.path.exists(data_d):
                        os.makedirs(data_d)
                    os.system("cp "+fi_d+" "+data_d)
                    if options.rmoriginalfile:
                        os.system("rm "+fi_d)