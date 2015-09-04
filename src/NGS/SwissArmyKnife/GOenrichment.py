from optparse import OptionParser
import os
import re, sys, time
from scipy import stats
from scipy.stats import hypergeom
from NGS.BasicUtil import *


'''
Created on 2015-4-10

@author: liurui
'''
parser = OptionParser()
parser.add_option("-g", "--gotablefile", dest="gotablefile", help="gotable title with :Ensembl Gene ID    Ensembl Transcript ID    GO Term Accession    GO Term Evidence Code    GO domain    GO Term Name    GO Term Definition,order and upper/lower case is arbitrarily")
parser.add_option("-G","--genelist",dest="genelist",default=None)
parser.add_option("-T","--trscptlist",dest="trscptlist",default=None)
parser.add_option("-o","--outpre",dest="outpre")

(options, args) = parser.parse_args()
if __name__ == '__main__':
    geneUtil.GOenrichment(options.gotablefile,options.outpre,options.genelist=None,options.trscptlist=None):

