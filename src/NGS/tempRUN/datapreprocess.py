'''
Created on 2019年1月8日
pick up each data unit into one folder
@author: RuiLiu
'''
from optparse import OptionParser
import os,re
# from NGS.Analysis.IntervalValueToBoxplot import path

parser = OptionParser()
parser.add_option("-f", "--configfile", dest="configfile",help="configfile ")
parser.add_option("-s", "--suffix", dest="suffix",action="append",default=["1.fq.gz","2.fq.gz"],help="suffix for a folder to package")
(options, args) = parser.parse_args()

cf=open(options.configfile,'r')

filepath=cf.readline().strip()
files=os.listdir(filepath)
folderMap={}#{"AS10_":[AS10_1.fq.gz,AS10_1.fq.gz]}
os.system("cp "+os.path.join(filepath,"md5.txt")+" "+os.path.join(filepath,"md5.txt.copy"))

if __name__ == '__main__':
    for fi in files:
        l_pre=len(fi)-len(options.suffix[0])
        fi_d=os.path.join(filepath,fi)
        if os.path.isdir(fi_d):
            pass
        else:
            for s in options.suffix:
                if re.search(r"."+s,fi)!=None:
                    data_d=os.path.join(filepath,fi[:l_pre])
                    if not os.path.exists(data_d):
                        os.makedirs(data_d)
                    os.system("cp "+fi_d+" "+data_d)