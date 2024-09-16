import os
from urllib.request import urlretrieve
from queue import Queue
from threading import Thread


def setup_download_dir(directory=None):
    """Setup the directory to download files to. If no path/directory specified,
    defaults to current working directory.
    """
    if directory is not None:
        if not os.path.exists(directory):
            os.makedirs(directory)
    else:
        directory = os.getcwd()


def download_file(directory, link, verbose=True):
    """Download a given link/file to the given download directory"""
    # Get the file name of the link to download
    file = os.path.basename(link)
    
    # Check that the download directory exists, and if it doesn't create it
    setup_download_dir(directory)
    
    # Download the file to the download directory
    download_path = "/".join((directory, file))
    if verbose:
        print(f"Downloading {link} to {download_path} \n")
    urlretrieve(link, download_path)


def ftp_download(ftp, file, directory, verbose=True):
    """Download a given FTP link to the given directory."""
    filename = os.path.basename(file)
    
    # Check that the download directory exists, and if it doesn't create it
    setup_download_dir(directory)
    
    # Download the file to the download directory
    download_path = os.path.join(directory, filename)
    if verbose:
        print(f"Downloading {file} to {download_path} \n")
        
    with open(download_path, 'wb') as f:
        ftp.retrbinary('RETR %s' % file, f.write)


class DownloadWorker(Thread):
    """Class utilizing multithreading in python to speed up downloads."""
    def __init__(self, queue):
        Thread.__init__(self)
        self.queue = queue
        
    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            directory, link = self.queue.get()
            try:
                download_file(directory, link)
            finally:
                self.queue.task_done()


class FTPWorker(Thread):
    """Class utilizing multithreading in python to speed up downloads."""
    def __init__(self, queue):
        Thread.__init__(self)
        self.queue = queue
        
    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            ftp, directory, link = self.queue.get()
            try:
                ftp_download(ftp, directory, link)
            finally:
                self.queue.task_done()