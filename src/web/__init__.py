__all__=["Action","DBA","Entity","Service"]
from flask import Flask
import config
from flask_bootstrap import Bootstrap
web = Flask(__name__,
            template_folder="templates",#by default
            static_folder=".",
            static_url_path="")
Bootstrap(web)

web.config.from_pyfile('../config.py')