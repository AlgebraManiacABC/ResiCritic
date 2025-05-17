from pathlib import Path
import tempfile
import pickle
import os

def save_log(data, filename: str):
    dir_name = Path(filename).absolute().parent.resolve()
    with tempfile.NamedTemporaryFile('wb', delete=False, dir=dir_name, prefix='._', suffix='.tmp') as tmp_file:
        pickle.dump(data,tmp_file)
        tmp_name = tmp_file.name
    os.replace(tmp_name, filename)

def init_log(filename: str) -> dict:
    log_dict = {
        'step': 'init',
    }
    try:
        Path(filename).touch()
        save_log(log_dict, filename)
    except Exception as e:
        print(f"Unexpected error while initializing log! {e}")
        return None
    return log_dict

def load_log(filename: str) -> dict:
    filepath = Path(filename)
    if filepath.exists() and filepath.is_file():
        with open(filename,'rb') as f:
            return pickle.load(f) or init_log(filename)
    return init_log(filename)

def filename_to_log(filename: str) -> str:
    return Path(filename).with_suffix('.pkl')
