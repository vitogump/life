'''
Created on 2014-11-11

@author: liurui
'''
from optparse import OptionParser

from src.pipelinecontrol.Util import *


#just fo
parser = OptionParser()



#scriptDir,mode="series",logfile
parser.add_option("-d", "--scriptDir", dest="scriptDir",help="scriptDir")
parser.add_option("-s", "--scripts", dest="scripts",action="append",help="oneline scriptexamplefile not in use so far")
parser.add_option("-p","--logicalpurpose",dest="logicalpurpose",help="a little note message")
# parser.add_option("-l", "--logfile", dest="logfile", help="bam bai sam sorted.bam vcf blast and so on. note this is just used in the cmdline output parameter")
parser.add_option("-t", "--ThreadsNum", dest="ThreadsNum",help="p:parallel s:series")

                                                                                                                                                          
(options, args) = parser.parse_args()

if __name__ == '__main__':
#     jk=myJobTracker(scriptDir=options.scriptDir,NumOfThread=int(options.ThreadsNum))
    callsh_updateDB(scriptDir=options.scriptDir,NumOfThread=int(options.ThreadsNum),logicalpurpose=options.logicalpurpose)
    print("finish")