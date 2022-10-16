#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2022-01-25

@author: R. Brian Akins - damiandudka0@gmail.com

Includes methods to run FREEDA from a frozen state inside a pyinstaller bundle.

"""

import os
import sys


def is_bundled():
    return getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS')


# PYINSTALLER: Update os path variable with location of frozen repository.
# Using https://stackoverflow.com/a/44352931
def resource_path(relative_path):
    """Get absolute path to given relative resource. Works for dev and for PyInstaller."""
    # If variable _MEIPASS exists, it stores the location in which the program is running. Otherwise, return input.
    base_path = getattr(sys, "_MEIPASS", "")
    return os.path.join(base_path, relative_path)

