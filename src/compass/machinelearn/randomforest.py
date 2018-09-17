'''
Created on 2018年9月5日

@author: Dr.liu
'''

from sklearn import cross_validation, metrics
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import accuracy_score, make_scorer
from sklearn.model_selection import train_test_split
import numpy as np
import pandas as pd


# from sklearn.model_selection import GridSearchCV# for second part
def encode_features(df_train,df_test):
    features=["色泽","根蒂","敲声","纹理","脐部","触感"]
    df_combined=pd.concat([df_train[features], df_test[features]])
    for feature in features:
        le=preprocessing.LabelEncoder()
        le=le.fit(df_combined[feature])
        df_train[feature] = le.transform(df_train[feature])
        df_test[feature] = le.transform(df_test[feature])
    return df_train,df_test
def simplify_interval_info(df):
    bins_density = (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
    bins_sugar = (0, 0.1, 0.2, 0.3, 0.4, 0.5)
    
    group_name_density = [0, 1, 2, 3, 4, 5, 6, 7]
    group_name_sugar = [0, 1, 2, 3, 4]
    
    category_density = pd.cut(df['密度'], bins_density, labels=group_name_density)
    categroy_sugar = pd.cut(df['含糖率'], bins_sugar, labels=group_name_sugar)
    
    df['密度'] = category_density
    df['含糖率'] = categroy_sugar
    
    return df
if __name__ == '__main__':
#     train=pd.read_csv("D:\\mechinelearn\\train_modified.csv")
#     target='Disbursed'#Disbursed的值就是二元分类的输出
#     IDcol='ID'
# #     print(train['Disbursed'].value_counts())
#     x_columns=[x for x in train.columns if x not in [target,IDcol]]
#     X=train[x_columns]
#     y=train['Disbursed']
#     print(y)
#     rf0=RandomForestClassifier(oob_score=True,random_state=10)
#     rf0.fit(X,y)
#     print(rf0.oob_score_)
#     
#     exit()
#second part 
    train_data=pd.read_csv("D:\\mechinelearn\\watermelon.csv")
    test_data=pd.read_csv("D:\\mechinelearn\\watermelon.csv")
    features=["色泽","根蒂","敲声","纹理","脐部","触感"]
    
    train_data, test_data = encode_features(train_data, test_data)
    train_data = simplify_interval_info(train_data)
    test_data = simplify_interval_info(test_data)

    
    
    
      
    X_all = train_data.drop(['好瓜'], axis=1)
    y_all = test_data['好瓜']
    y_result = [1,0,0]
    print("====X_all=====")
    print(X_all)
    print("====y_all=====")
    print(y_all)
    num_test=0.8
    X_train,X_test,y_train,y_test=train_test_split(X_all, y_all, test_size=num_test, random_state=3)
    #Choose some parameter combinations to try
    parameters={'n_estimators':[5,6,7],'criterion':['entropy','gini']}
    #type of scoring used to compare parameter combinations
    acc_scorer=make_scorer(accuracy_score)
    clf = RandomForestClassifier()
    #Run the grid search
    print("y_train")
    print(y_train)
    grid_obj=GridSearchCV(clf,parameters,scoring=acc_scorer)

    
    grid_obj=grid_obj.fit(X_train,y_train)
    #set the clf to the best combination of parameters
    clf=grid_obj.best_estimator_
    clf=clf.fit(X_train,y_train)
    print("X_test\n",X_test)
    test_predictions=clf.predict(X_test)
    print("test_predictions\n",test_predictions)
    print("测试集准确率:  %s " % accuracy_score(y_test, test_predictions))
    
    print(test_data)
    print(y_test)
    
    predictions = clf.predict(test_data.drop(['好瓜'], axis=1))
    print("最终准确率:  %s " % accuracy_score(y_result, predictions))