from os.path import dirname,join,abspath
import sys
import hashlib
from functools import partial

def get_platform():
    platforms = {
        'linux': 'Linux',
        'linux1': 'Linux',
        'linux2': 'Linux',
        'darwin': 'Mac',
        'win32': 'Windows'
    }
    if sys.platform not in platforms:
        return sys.platform

    return platforms[sys.platform]


def md5sum(filename):
    with open(filename, mode='rb') as f:
        d = hashlib.md5()
        for buf in iter(partial(f.read, 128), b''):
            d.update(buf)
    return d.hexdigest()

##Returnt the path of the Rscript directory in the project
def find_rscript():
    return join(dirname(dirname(abspath(__file__))),"Rscript")

def find_data():
    return join(dirname(dirname(abspath(__file__))),"data")
