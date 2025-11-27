from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build
import os

# Scopes required for Drive access
SCOPES = ['https://www.googleapis.com/auth/drive']


def get_drive_service():
    creds = Credentials.from_service_account_file(
        os.environ.get("GOOGLE_APPLICATION_CREDENTIALS"), scopes=SCOPES)
    service = build('drive', 'v3', credentials=creds)
    return service

# Example: List files in Google Drive
def list_files():
    service = get_drive_service()
    results = service.files().list(pageSize=10, fields="nextPageToken, files(id, name)").execute()
    items = results.get('files', [])
    if not items:
        print('No files found.')
    else:
        print('Files:')
        for item in items:
            print(u'{0} ({1})'.format(item['name'], item['id']))

if __name__ == '__main__':
    list_files()
