
#import Make_Picture
from NGS.RUtil import *
from optparse import OptionParser

'''
Created on 2013-8-11

@author: rui
'''
parser = OptionParser()

parser.add_option("-o","--pathoutputfilename",dest="pathoutputfilename",help="default infile2_infile1")#
parser.add_option("-P","--positive_withgenename",dest="multiple_positive_winfiles_withgenename",action="append",nargs=2,default=[],help="on top,file name and threshold")#
parser.add_option("-N","--negtive_withgenename",dest="multiple_negtive_winfiles_withgenename",action="append",nargs=2,default=[],help="at bottom,file name and threshold")#

parser.add_option("-p","--positive",dest="multiple_positive_winfiles",action="append",default=[],help="on top")#
parser.add_option("-n","--negtive",dest="multiple_negtive_winfiles",action="append",default=[],help="at bottom")#
parser.add_option("-c","--columnname",dest="columnname",default="zvalue")

(options, args) = parser.parse_args()

columnname=options.columnname.strip()
if (options.multiple_negtive_winfiles!=[] or options.multiple_positive_winfiles!=[]) and (options.multiple_negtive_winfiles_withgenename!=[] or options.multiple_positive_winfiles_withgenename!=[]):
    print("require all winfile with gene or all winfile without gene")
    exit(-1)
if __name__ == '__main__':
    makeMhtGraph = Make_Picture.MakeMhtGraph()
    if options.multiple_positive_winfiles_withgenename!=[]:
        for p_inputfileName,t in options.multiple_positive_winfiles_withgenename[:]:
            makeMhtGraph.makeHistonPicture(p_inputfileName, "Fst")#,"c(0,2000)","c(0,45)"
    if options.multiple_negtive_winfiles_withgenename!=[]:
        for n_inputfileName,t in options.multiple_negtive_winfiles_withgenename[:]:
            makeMhtGraph.makeHistonPicture(n_inputfileName, "Hp")#,"c(0,2000)","c(0,45)"
    
    if options.multiple_positive_winfiles!=[]:
        for p_inputfileName in options.multiple_positive_winfiles[:]:
            makeMhtGraph.makeHistonPicture(p_inputfileName, "Fst")#,"c(0,2000)","c(0,45)"
    if options.multiple_negtive_winfiles!=[]:
        for n_inputfileName in options.multiple_negtive_winfiles[:]:
            makeMhtGraph.makeHistonPicture(n_inputfileName, "Hp")#,"c(0,2000)","c(0,45)"
    if (options.multiple_negtive_winfiles_withgenename!=[] or options.multiple_positive_winfiles_withgenename!=[]):
        makeMhtGraph.makeMhtplots_compareInOnePicture_withgeneName(options.pathoutputfilename, options.multiple_positive_winfiles_withgenename, options.multiple_negtive_winfiles_withgenename, 0,columnname)
    else:
        makeMhtGraph.makeMhtplots_compareInOnePicture(options.pathoutputfilename, options.multiple_positive_winfiles, options.multiple_negtive_winfiles, 0,columnname)
