import sys
#import Make_Picture
from NGS.RUtil import *
'''
Created on 2013-8-11

@author: rui
'''
if len(sys.argv) < 5:
    print("python test.py [inputfile1] [inputfile2] [inputfile3]....[chromPrefix] [postive_negtive(or a)][dataType]")
    exit(-1)
dataType=sys.argv[-1]
if sys.argv[-2]!= "a":
    
    postive_negtive=sys.argv[-2]
    print(postive_negtive)
else:
    postive_negtive=None
chromPrefix=sys.argv[-3]
if __name__ == '__main__':
    for inputfileName in sys.argv[1:-3]:
        makeMhtGraph = Make_Picture.MakeMhtGraph()
        if sys.argv[-2]!= "a":
            makeMhtGraph.makeMhtPicture_HistonPicture(inputfileName,dataType,chromPrefix,postive_negtive,fillvalue='NA')
        else:
            makeMhtGraph.makeMhtPicture_HistonPicture(inputfileName,dataType,chromPrefix,postive_negtive=None,fillvalue='NA')