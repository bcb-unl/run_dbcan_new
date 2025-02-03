import requests
import os
import tarfile
from tqdm import tqdm

class DBDownloader:
    def __init__(self, config):
        self.config = config
        self.db_dir = config.get('db_dir')
        self.filename = os.path.join(self.db_dir, "dbCAN_db.tar.gz")
        if not os.path.exists(self.db_dir):
            os.makedirs(self.db_dir)

    def download_file(self):
        with requests.get("https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/dbCAN_db.tar.gz", stream=True) as response:
            response.raise_for_status()
            total_size = int(response.headers.get('content-length', 0))
            with open(self.filename, 'wb') as f, tqdm(
                desc=self.filename,
                total=total_size,
                unit='iB',
                unit_scale=True,
                leave=True
            ) as bar:
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        f.write(chunk)
                        bar.update(len(chunk))

    def extract_tar_file(self):
        with tarfile.open(self.filename, 'r:gz') as tar:
            def member_filter(member):
                if member.isdir():
                    return None
                if '/' in member.name:
                    member.name = '/'.join(member.name.split('/')[1:])
                return member if member.name else None

            filtered_members = [m for m in tar.getmembers() if member_filter(m)]
            tar.extractall(path=self.db_dir, members=filtered_members)
        print(f"Extracted {self.filename} to {self.db_dir}")
        os.remove(self.filename)
        print(f"Deleted the archive: {self.filename}")

        self.remove_empty_dirs(self.db_dir)

    def remove_empty_dirs(self, path):
        for root, dirs, files in os.walk(path, topdown=False):
            for dir in dirs:
                full_dir_path = os.path.join(root, dir)
                if not os.listdir(full_dir_path):
                    os.rmdir(full_dir_path)
                    print(f"Removed empty directory: {full_dir_path}")



# if __name__ == '__main__':
#     config = {
#         'db_dir': './dbCAN_databases'

#     }
#     db_downloader = DBDownloader(config)
#     db_downloader.download_file()
#     db_downloader.extract_tar_file()
