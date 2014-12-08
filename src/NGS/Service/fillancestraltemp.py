'''
Created on 2014-12-1

@author: liurui
'''
from optparse import OptionParser


from NGS.Service.Ancestralallele import AncestralAlleletabletools

parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."
parser.add_option("-a", "--archicpopVcfFile", dest="archicpopVcfFile",help="oneline scriptexamplefile")
# parser.add_option("-o", "--outputpath", dest="outputpath", help="outputpath")
parser.add_option("-c", "--chromtable", dest="chromtable", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")
parser.add_option("-d", "--depthfile", dest="depthfile", help="winvalue or zvalue")
parser.add_option("-l", "--interceptdirs", dest="interceptdirs", help="winvalue or zvalue")
parser.add_option("-t","--toplevelsnptable",dest="toplevelsnptable",default="0",help="depth of the folder to output")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
archicpopVcfFile=options.archicpopVcfFile
chromtable=options.chromtable
toplevelsnptablename=options.toplevelsnptable
archicpopNameindepthFile=options.interceptdirs
depthFile=options.depthFile
if __name__ == '__main__':
    ancestralalleletabletools=AncestralAlleletabletools(database="ninglabvariantdata_tmp", ip="10.2.48.140", usrname="root", pw="1234567")
    ancestralalleletabletools.fillAncestral(archicpopVcfFile=archicpopVcfFile.strip(),depthFile,archicpopNameindepthFile=archicpopNameindepthFile,chromtable=chromtable.strip(),toplevelsnptablename=toplevelsnptablename)