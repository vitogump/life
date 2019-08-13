'''
Created on 2019年8月12日

@author: RuiLiu
'''
import pysam
import numpy

parser = OptionParser()
parser.add_option("-c", "--bamfile", dest="bamfile",help="write report to FILE")
parser.add_option("-r", "--reffa", dest="reffa",help="far sure but with few locs")
parser.add_option("-o","--outfileprename",dest="outfilepreName",help="outfilepreName with path")
parser.add_option("-a","--ancestralspeciescolname",dest="ancestralspeciescolname",help="ancestralspecisname")
parser.add_option("-m","--mindepth",dest="mindepth",help="mindepth for both archicpop and ancestralallel")#
# parser.add_option("-e", "--bedfile", action="append",dest="bedfile",default=[],help="measure region of bedfile")
parser.add_option("-s","--speciesesName",action="append",dest="speciesesName",default=[],help="speciesName in table")#
                                                                                                                                                          
(options, args) = parser.parse_args()
if __name__ == '__main__':
    bamfile=pysam.Samfile(options.bamfile.strip(),'rb')
    bamfile=pysam.AlignmentFile(options.bamfile.strip(),'rb')