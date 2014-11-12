'''
Created on 2014-11-11

@author: liurui
'''
from optparse import OptionParser

from src.pipelinecontrol.Util import JobTracker


#just fo
parser = OptionParser()



#scriptDir,mode="series",logfile
parser.add_option("-d", "--scriptDir", dest="scriptDir",help="scriptDir")
# parser.add_option("-o", "--outputpath", dest="outputpath", help="outputpath")

parser.add_option("-l", "--logfile", dest="logfile", help="bam bai sam sorted.bam vcf blast and so on. note this is just used in the cmdline output parameter")
parser.add_option("-m", "--mode", dest="mode",help="parallel series")

                                                                                                                                                          
(options, args) = parser.parse_args()
if options.mode.lower()=="s":
    print("series")
    mode="series"
elif options.mode.lower()=="p":
    print("parallel")
    mode="parallel"
if __name__ == '__main__':
    jk=JobTracker(scriptDir=options.scriptDir,mode=mode,logfile=options.logfile)
    jk.run()
    print("finish")