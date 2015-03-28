
#import Make_Picture
from NGS.RUtil import *
from optparse import OptionParser

'''
Created on 2013-8-11

@author: rui
'''
parser = OptionParser()

parser.add_option("-o","--pathoutputfilename",dest="pathoutputfilename",help="default infile2_infile1")#

parser.add_option("-p","--positive",dest="multiple_positive_winfiles",action="append",default=[],help="on top")#
parser.add_option("-n","--negtive",dest="multiple_negtive_winfiles",action="append",default=[],help="at bottom")#
(options, args) = parser.parse_args()
#if len(sys.argv) < 5:
#    print("python test.py [inputfile1] [inputfile2] [inputfile3]....-p [chromPrefix] -s [positive_negtive(or a)] -T [dataType] -c column_num")
#    exit(-1)


if __name__ == '__main__':
    makeMhtGraph = Make_Picture.MakeMhtGraph()
    if options.multiple_positive_winfiles!=[]:
        for p_inputfileName in options.multiple_positive_winfiles[:]:
            makeMhtGraph.makeHistonPicture(p_inputfileName, "Fst","c(0,2000)","c(0,45)")
    if options.multiple_negtive_winfiles!=[]:
        for n_inputfileName in options.multiple_negtive_winfiles[:]:
            makeMhtGraph.makeHistonPicture(n_inputfileName, "Hp","c(0,2000)","c(0,45)")
    makeMhtGraph.makeMhtplots_compareInOnePicture(options.pathoutputfilename, options.multiple_positive_winfiles, options.multiple_negtive_winfiles, 0)
#         if options.positive_negtive_a!= "a":
#             makeMhtGraph.makeMhtPicture_HistonPicture(inputfileName,dataType,positive_negtive,fillvalue=options.fillvalue)
#         else:
#             makeMhtGraph.makeMhtPicture_HistonPicture(inputfileName,dataType,positive_negtive=None,fillvalue=options.fillvalue)