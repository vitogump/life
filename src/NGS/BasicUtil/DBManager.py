import mysql.connector
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
    def operateDB(self,sqltype,*sqls):
        if not self.conn:
            self.connect()
        try:
            cursor=self.conn.cursor()
            result=[]
            if sqltype=='select':
                cursor.execute(sqls[0])
                result=cursor.fetchall()
                cursor.close()
                return result
            else:
                for sql in sqls:
                    cursor.execute(sql)
        except mysql.connector.Error as e:
            print('query error!{}'.format(e))
        #self.disconnect()