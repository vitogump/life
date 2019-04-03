__all__=["Action","DBA","Entity","Service"]
from flask import Flask
import config
from flask_bootstrap import Bootstrap
web = Flask(__name__)
Bootstrap(web)

web.config.from_pyfile('../config.py')