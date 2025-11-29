import os
import json
from pydrive2.auth import GoogleAuth
from pydrive2.drive import GoogleDrive

def authenticate_gdrive():
    """Authenticate with Google Drive using service account credentials."""
    gauth = GoogleAuth()
    
    # Check if credentials are provided as environment variable (for GitHub Actions)
    creds_json = os.environ.get('GDRIVE_CREDENTIALS')
    
    if creds_json:
        # Parse JSON credentials from environment variable
        creds_dict = json.loads(creds_json)
        
        # Save credentials temporarily
        with open('credentials.json', 'w') as f:
            json.dump(creds_dict, f)
        
        # Configure settings for service account
        gauth.settings['client_config_backend'] = 'file'
        gauth.settings['client_config_file'] = 'credentials.json'
        gauth.settings['service_config'] = {
            'client_json_file_path': 'credentials.json',
        }
        
        # Authenticate with service account
        gauth.ServiceAuth()
        
        # Clean up temporary file
        os.remove('credentials.json')
    else:
        # Local authentication (interactive)
        gauth.LocalWebserverAuth()
    
    return GoogleDrive(gauth)

def list_files(drive, max_results=100, folder_id=None):
    """
    List files from Google Drive.
    
    Args:
        drive: Authenticated GoogleDrive instance
        max_results: Maximum number of files to retrieve
        folder_id: Optional folder ID to list files from specific folder
    """
    # Build query
    if folder_id:
        query = f"'{folder_id}' in parents and trashed=false"
    else:
        query = "trashed=false"
    
    # List files
    file_list = drive.ListFile({
        'q': query,
        'maxResults': max_results
    }).GetList()
    
    print(f"Found {len(file_list)} files:\n")
    
    for file in file_list:
        print(f"Title: {file['title']}")
        print(f"ID: {file['id']}")
        print(f"Type: {file['mimeType']}")
        print(f"Modified: {file['modifiedDate']}")
        print("-" * 50)
    
    return file_list

def main():
    """Main function to authenticate and list files."""
    try:
        # Authenticate
        print("Authenticating with Google Drive...")
        drive = authenticate_gdrive()
        print("Authentication successful!\n")
        
        # List files
        files = list_files(drive)
        print(f"\nTotal files listed: {len(files)}")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        raise

if __name__ == "__main__":
    main()




# from pydrive2.auth import GoogleAuth
# from pydrive2.drive import GoogleDrive

# gauth = GoogleAuth()
# gauth.LocalWebserverAuth()

# drive = GoogleDrive(gauth)



# file1 = drive.CreateFile({'title': 'Hello1.txt'})
# file1.SetContentString('Hello1')
# file1.Upload() # Files.insert()
