'''
Created on 2017年10月19日

@author: liurui
'''
from web import web
from flask import request,jsonify,send_from_directory,abort,render_template
import os
from werkzeug.routing import BaseConverter
class RegexConverter(BaseConverter):
    def __init__(self, map, *args):
        self.map = map
        self.regex = args[0]
web.url_map.converters['regex'] = RegexConverter

@web.route('/')
def index():
    if "/"=="/":
        return render_template('downloadlist.html')
        print("1",os.getcwd())
    else:
        print(os.getcwd())
        return send_from_directory(os.getcwd(),'../../toDownload/pythonstudy/MITjisuanjikexuebianchengdaolunlec01.mp4',as_attachment=True)