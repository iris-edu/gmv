import sys
import os
from urllib.request import urlopen

"""     
Name:
    gmv_param.py
    
Description:
    GMV Parameter file
    
History:
    2020-09-03 Manoch: V.2020.241 R1.1 Public release.
    2020-05-11 Manoch: Initial release.
"""

version = 'V.2020.241'


class ObjDict(dict):
    """Accessing dictionary items as object attributes:
            https://goodcode.io/articles/python-dict-object/
    """

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: {}".format(name))

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: {}".format(name))


def print_message(flag, text, file):
    """Print out a message. Force the flush and write to the file handle"""

    if flag == 'ERR':
        print(f'\n\n{60 * "="}\n', file=sys.stderr, flush=True)
        print(f'[{flag}] {text}', sep='\n', file=sys.stderr, flush=True)
        print(f'\n\n{60 * "="}\n', file=sys.stderr, flush=True)
    elif file is not None:
        print(f'[{flag}] {text}', file=file, flush=True)
    else:
        print(f'[{flag}] {text}', flush=True)


def read_url(target_url, file, verbose=False):
    """Read content of a URL."""
    if verbose:
        print_message('INFO', f'Opening URL: {target_url}', file)

    with urlopen(target_url) as url:
        content = url.read().decode()
    return content


def get_log_file(log_file_name, region, comp, log_to_screen):
    """Open the log file."""
    if log_to_screen:
        log_file = sys.stdout
    else:
        log_file = open(f'{log_file_name}_{region}_{comp}.log', 'a')
    return log_file


def mkdir(target_directory, file=None):
    """ Make a directory if it does not exist."""
    directory = None
    try:
        directory = target_directory
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory
    except Exception as _er:
        print_message('ERR', f'failed to create directory {directory}\n{_er}', file)
        return None
