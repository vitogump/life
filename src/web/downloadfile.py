'''
Created on 2017年10月19日

@author: liurui
'''
from web import web
from web import views




#     if request.method=='GET':
#         if os.path.isfile(filename):
#             return send_from_directory(filename, as_attachment=True)
# def server_static(urlpath):
#     print("sever_static")
#     herfs=""""""
#     catajudge=re.search(r'[^/]*$',urlpath).group(0)
#     if urlpath.endswith("GBS"):
#         path='../../toDownload/classical_paper/GBS/'
#     elif urlpath.endswith("Python_study"):
#         path='../../toDownload/pythonstudy/'
#     elif urlpath.endswith(".html"):
#         print("sssss")
#         return static_file(urlpath,root='../../index')
#     else:
#         print(unquote(urlpath))
#         path = "../../"+re.search(r'catalog/(.*)',unquote(urlpath)).group(1)+"/"
#         print(path)
#    else:#download file
#        path="../../"+re.search(r'^.*/',unquote(urlpath)).group(0)
#        print(urlpath,re.search(r'^.*/',urlpath).group(0),re.search(r'[^/]*$',urlpath).group(0))
#        return static_file(re.search(r'[^/]*$',urlpath).group(0), root='../../classical_paper/'+re.search(r'^.*/',urlpath).group(0),download=re.search(r'[^/]*$',urlpath).group(0))
#     l=os.listdir(path=path)
#     print(l,path)
#     for a in l:
#         if os.path.isdir(path+a) and not urlpath.endswith("rar"):
#             url=quote('/download/catalog/'+re.search(r'\.\./\.\./(.*)',path+a).group(1))
#             herfs+="""<a href="""+url+""" target="_self">"""+a+"""</a></p>
#             """
#         elif os.path.isfile(path+a):
#             url=quote('/downloadfile/'+re.search(r'\.\./\.\./(.*)',path+a).group(1))
#             herfs+="""<a href="""+url+""" target="_self">"""+a+"""</a></p>
#             """
#     return herfs     
if __name__ == '__main__':
    web.run(host="localhost",port=8080,debug = True)