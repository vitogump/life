'''
Created on 2014-12-1

@author: liurui
'''
parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."
parser.add_option("-a", "--archicpopVcfFile", dest="archicpopVcfFile",help="oneline scriptexamplefile")
# parser.add_option("-o", "--outputpath", dest="outputpath", help="outputpath")
parser.add_option("-c", "--chromtable", dest="chromtable", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")

parser.add_option("-t","--toplevelsnptable",dest="toplevelsnptable",default="0",help="depth of the folder to output")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
if __name__ == '__main__':
    pass