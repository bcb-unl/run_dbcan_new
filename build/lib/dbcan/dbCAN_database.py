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
        response = requests.get("https://bcb.unl.edu/dbCAN2/download/test/dbCAN_db.tar.gz", stream=True)
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
            tar.extractall(path=self.db_dir)
        print(f"Extracted {self.filename} to {self.db_dir}")



# if __name__ == '__main__':
#     config = {
#         'db_dir': './dbCAN_databases'

#     }
#     db_downloader = DBDownloader(config)
#     db_downloader.download_file()
#     db_downloader.extract_tar_file()
