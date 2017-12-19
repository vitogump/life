__all__=["Action","DBA","Entity","Service","config"]
from flask import Flask
from web import config
web = Flask(__name__)
from web import views
web.config.from_pyfile('config.py')