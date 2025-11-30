from pydrive2.auth import GoogleAuth
from pydrive2.drive import GoogleDrive

gauth = GoogleAuth()
gauth.LocalWebserverAuth()

drive = GoogleDrive(gauth)



file1 = drive.CreateFile({'title': 'Hello1.org'})
file1.SetContentString('* Hello1')
file1.Upload() # Files.insert()
