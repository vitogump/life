# -*- coding: UTF-8 -*-
'''
Created on 2014-11-27

@author: liurui
'''
from optparse import OptionParser

from NGS.BasicUtil import VCFutil , Util
from NGS.Service.Ancestralallele import AncestralAlleletabletools
from src.pipelinecontrol.Util import upTodownTravelDir, \
    OperatorWithData_loadintodatabase


parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."
parser.add_option("-i", "--vcffiledir", dest="vcffiledir",help="oneline scriptexamplefile")
parser.add_option("-s","--suffix",dest="vcfsuffix",help="suffix")
# parser.add_option("-o", "--outputpath", dest="outputpath", help="outputpath")
parser.add_option("-d", "--datadepth", dest="datadepth", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")
parser.add_option("-p", "--droptable", dest="droptable", action="store_true",default=False,help="")
parser.add_option("-I","--Interceptor_depth",dest="Interceptor_depth",default="0",help="depth of the folder to output")
parser.add_option("-l", "--interceptdirs", dest="interceptdirs",action="append", default=[], help="winvalue or zvalue")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
Interceptor_depth=int(options.Interceptor_depth)
datadepth=int(options.datadepth)
inputdatafilesrootpath=options.vcffiledir
interceptdirs=options.interceptdirs
vcfsuffix=options.vcfsuffix
if __name__ == '__main__':
    ancestralalleletabletools=AncestralAlleletabletools(database="ninglabvariantdata", ip=Util.ip, usrname=Util.username, pw=Util.password)
    operatorWithData_loadintodatabase=OperatorWithData_loadintodatabase(inputdatafilesrootpath,ancestralalleletabletools, interceptdirs,vcfsuffix,options.droptable)
    upTodownTravelDir(inputdatafilesrootpath,operatorWithData_loadintodatabase,datadepth,Interceptor_depth,0,datadepth)
    print("finish")