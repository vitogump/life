from mysql.connector import errorcode
import mysql.connector
import re
import time
'''
Created on 2013-8-22

@author: liurui
'''
SLEEP_FOR_NEXT_TRY = 3
class DBTools():
    '''
    classdocs
    '''
    

    def __init__(self, host, user, passwd, db):

        self.host = host
        self.user = user
        self.passwd = passwd
        self.db = db
        self.conn = None
        '''
        Constructor
        '''
    def connect(self):
        if self.conn:
            return
        while True:
            try:
                self.conn=mysql.connector.connect(host=self.host,user=self.user,passwd=self.passwd,database=self.db)
                break
            except mysql.connector.Error as e:
                print('connect fails!{}'.format(e))
                print("sleep %d seconds for next try"% SLEEP_FOR_NEXT_TRY)
                time.sleep(SLEEP_FOR_NEXT_TRY)
    
    def disconnect(self):
        if not self.conn:
            return
        try:
            self.conn.close()
            self.conn=None
        except:
            print("conn can't close ")
    def operateDB(self,sqltype,*sqls,data=None):
        if not self.conn:
            self.connect()
        try:
            cursor=self.conn.cursor()
            result=[]
            if sqltype=='select':
                cursor.execute(sqls[0])
                result=cursor.fetchall()
                return result
            elif sqltype=='updata' or sqltype=='insert':
                for sql in sqls:
                    if data==None:
                        cursor.execute(sql)
                    else:
                        cursor.execute(sql,data)
                    self.conn.commit()
            elif sqltype=="callproc":
                cursor.callproc(sqls[0],data)
        except mysql.connector.Error as e:
            print('query error!{}'.format(e))
        finally:
            cursor.close()
    def create_table(self,table):
        """  创建表
             TABLES = {}
                TABLES['employees'] = (
                    "CREATE TABLE `employees` ("
                    " `emp_no` int(11) NOT NULL AUTO_INCREMENT,"
                    " `birth_date` date NOT NULL,"
                    " `first_name` varchar(14) NOT NULL,"
                    " `last_name` varchar(16) NOT NULL,"
                    " `gender` enum('M','F') NOT NULL,"
                    " `hire_date` date NOT NULL,"
                    " PRIMARY KEY (`emp_no`)"
                    ")"
                    )"""
        if not self.conn:
            self.connect()
        cursor =self.conn.cursor()

        for key in table:
            try:
                print ("Creating table {}:".format(key))
                cursor.execute(table[key])
            except mysql.connector.Error as err:
                if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
                    print('already exist.')
                    return 0
                else:
                    print(err.errmsg)
            else:
                print('OK')
        cursor.close()
#    def load_file(self,sql):
#        print(sql)
#        if not self.conn:
#            self.connect()
#        try:
#            cursor = self.conn.cursor()
#            cursor.execute(sql)目前 mysql connector python 不支持这个load data local infile功能
#            self.conn.commit()
#        except mysql.connector.Error as e:
#            print("loadFile sql fails {}".format(e))
#        finally:
#            cursor.close()
#    def load_file(self,tableName,*tabletitles,fileName):
##        print(tabletitles)
##        exit()
#        if not self.conn:
#            self.connect()
#        cursor=self.conn.cursor()
#        filetoload=open(fileName,'r')
#        filedata=[]
#        for line in filetoload:
#            linelist=re.split(r'\s+',line.strip())
#            filedata.append((linelist[0],linelist[1],linelist[2],linelist[3],linelist[4],linelist[5]))
#        print(filedata)
#        try:
#            sql="insert into "+tableName+" "+str(tabletitles)+" values (%s,%s,%s,%s,%s,%s)"
#            print(sql)
#            cursor.executemany("insert into "+tableName+" "+str(tabletitles)+" values (%s,%s,%s,%s,%s,%s)",filedata)
#        except mysql.connector.Error as e:
#            print("loadFile sql fails {}".format(e))
#        filetoload.close()
#        cursor.close()
    
    def drop_table(self,tableName):
        if not self.conn:
            self.connect()
        cursor=self.conn.cursor()
        cursor.execute("drop table if exists "+tableName)
        cursor.close()
        
        
        
        
        
        
        #self.disconnect()