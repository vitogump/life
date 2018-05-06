'''
Created on 2018-3-19

@author: Dr.liu
'''
from web import web
from web import views
import socket

def get_host_ip():
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(('8.8.8.8', 80))
        ip = s.getsockname()[0]
    finally:
        s.close()
    return ip
if __name__ == '__main__':
    web.run(host=get_host_ip(),port=8314,debug=True)