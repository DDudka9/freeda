#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copyright 2022 - Damian Dudka and R. Brian Akins - contact: damiandudka0@gmail.com

This file is part of FREEDA.

FREEDA is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

FREEDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with FREEDA.
If not, see <https://www.gnu.org/licenses/>.

"""


"""

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

