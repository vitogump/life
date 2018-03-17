__all__=["Action","DBA","Entity","Service","config"]
from flask import Flask
from web import config
from flask_bootstrap import Bootstrap
web = Flask(__name__)
Bootstrap(web)
from web import views
web.config.from_pyfile('config.py')