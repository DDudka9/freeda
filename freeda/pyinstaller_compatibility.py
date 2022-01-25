#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2022-01-25

@author: brian

Includes methods to run FREEDA from a frozen state inside a pyinstaller bundle.

"""

import os
import sys


def is_bundled():
    return getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS')


def get_path(command):
    """If in a pyintstaller bundle, returns a path to the input in the 'include' folder. Otherwise, returns input."""
    if is_bundled():
        return os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "include", command)
    else:
        return command
