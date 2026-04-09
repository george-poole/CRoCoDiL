import shutil
import os

DEST_DIR = './_build/html/'

SOURCE_DIRS = (
    './system_a/',
    './system_aMu/',
    './system_aTheta/',
    './system_b/',
    './system_c/',
    './system_d/',
)

IGNORE = (
    '*.py', 
    '*.ipynb', 
    '*.pdf', 
    '*.h5', 
    '*.xdmf', 
    '*.npy', 
    '*.npz',
    '*.txt',
    '*full*',
)

if __name__ == "__main__":
    for s in SOURCE_DIRS:
        shutil.copytree(
            s,
            os.path.join(DEST_DIR, s),
            ignore=shutil.ignore_patterns(*IGNORE),
            symlinks=True,
            dirs_exist_ok=True,   
        )