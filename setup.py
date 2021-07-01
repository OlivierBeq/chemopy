# -*- coding: utf-8 -*-

"""Setup module."""

import os
import shutil

import setuptools

if __name__ == "__main__":

    # The following ensures external tools and config files
    # are being copied to the source code folder

    # List of folders under root of project directory
    ROOT_FOLDERS = ['MOPAC', 'external']
    # List of files under root of project directory
    ROOT_FILES = ['MOPAC.cfg', 'external.cfg']

    # Absolute path to the root folder of project
    HERE = os.path.dirname(__file__)

    # Relative path to src directory
    SRC_DIR = "./src/chemopy"
    SRC_PATH = os.path.realpath(os.path.join(HERE, SRC_DIR))

    # Store paths for reverse operation
    matchings = []
    here = os.path.realpath(HERE)
    for folder in ROOT_FOLDERS:
        if len(folder) == 0:
            raise RuntimeError(f'Folder does not exist: {folder}')
        match = {'from': here,
                 'to': SRC_PATH,
                 'name': folder}
        # Skip if folder already in the temporary folder being compiled
        # Needed for proper functioning of tox environments
        if os.path.exists(os.path.join(SRC_PATH, folder)):
            continue
        matchings.append(match)
    for file_ in ROOT_FILES:
        match = {'from': here,
                 'to': SRC_PATH,
                 'name': file_}
        # Skip if file already in the temporary folder being compiled
        # Needed for proper functioning of tox environments
        if os.path.exists(os.path.join(SRC_PATH, file_)):
            continue
        matchings.append(match)

    # Move files to source package temporarily for installation
    for match in matchings:
        from_, to_, name_ = match['from'], match['to'], match['name']
        shutil.move(os.path.join(from_, name_), to_)

    setuptools.setup()

    # Move folders and files back to original position
    for match in matchings:
        from_, to_, name_ = match['from'], match['to'], match['name']
        shutil.move(os.path.join(to_, name_), from_)
