# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
A type is considered "primitive" if it is a non-collection type which can
be serialized and deserialized successfully.
"""
from __future__ import print_function, division, absolute_import

from functools import wraps
from six import string_types, integer_types

NoneType = type(None)

PRIMITIVE_TYPES = (bool, float, NoneType) + string_types + integer_types

def return_primitive(fn):
    """
    Decorator which wraps a single argument function to ignore any
    arguments of primitive type (simply returning them unmodified).
    """
    @wraps(fn)
    def wrapped_fn(x):
        if isinstance(x, PRIMITIVE_TYPES):
            return x
        return fn(x)
    return wrapped_fn
