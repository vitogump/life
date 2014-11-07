
#import Make_Picture
from NGS.RUtil import *
from optparse import OptionParser

'''
Created on 2013-8-11

@author: rui
'''
parser = OptionParser()
parser.add_option("-p", "--chromPrefix", dest="chromPrefix",# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="chromPrefix + number")
parser.add_option("-s", "--positive_negtive_a", dest="positive_negtive_a",# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="write report to FILE")
# (options, args) = parser.parse_args()
parser.add_option("-T","--dataType",dest="dataType",help="default infile1_infile2")#
parser.add_option("-c","--column",dest="column",default=6,help="default infile2_infile1")#
parser.add_option("-f","--fillvalue",dest="fillvalue",default=0,help="default infile2_infile1")#
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
#if len(sys.argv) < 5:
#    print("python test.py [inputfile1] [inputfile2] [inputfile3]....-p [chromPrefix] -s [positive_negtive(or a)] -T [dataType] -c column_num")
#    exit(-1)
dataType=options.dataType
if options.positive_negtive_a!= "a":
    
    positive_negtive=options.positive_negtive_a
    print(positive_negtive)
else:
    positive_negtive=None
if options.fillvalue=="0":
    options.fillvalue=0
chromPrefix=options.chromPrefix
if __name__ == '__main__':
    for inputfileName in args[:]:
        makeMhtGraph = Make_Picture.MakeMhtGraph()
        if options.positive_negtive_a!= "a":
            makeMhtGraph.makeMhtPicture_HistonPicture(inputfileName,dataType,chromPrefix,positive_negtive,fillvalue=options.fillvalue)
        else:
            makeMhtGraph.makeMhtPicture_HistonPicture(inputfileName,dataType,chromPrefix,positive_negtive=None,fillvalue=options.fillvalue)