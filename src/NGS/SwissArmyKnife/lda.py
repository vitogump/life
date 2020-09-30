'''
Created on 2020年9月24日

@author: RuiLiu
'''
from optparse import OptionParser
import re
from scipy import stats
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis



        
        
if __name__ == '__main__':
    

    print("############perform lda analysis#####################")
    X=pd.read_table('KOlefsefracPvalue.rmtitle.KO.txt',header=None)
    print('select significant rec & trans to a new df')
    Xselect=X.loc[X[11]<=0.05]
#     Xfeaturetile=Xselect.T.loc[1:10]
#     Xfeaturetile.columns=Xselect.T.loc[0].tolist()
    #new group column
    y=pd.DataFrame({'group':['G','G','G','G','G','S','S','S','S','S']})
    y=['G','G','G','G','G','S','S','S','S','S']
    lda = LinearDiscriminantAnalysis()
    X_lda = lda.fit_transform(Xfeaturetile, y)
    print(X_lda)
    print('X_lda[:,0]')
    print(X_lda[:,0])
    print("X_lda[:,1]")
    print(X_lda[:,1])