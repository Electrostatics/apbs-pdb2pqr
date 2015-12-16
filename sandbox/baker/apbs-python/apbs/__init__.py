""" The goal of this module is to wrap as much APBS functionality as possible in Python. """
# TODO - add more docstrings with copyright, authors, etc.

import sys
if sys.version_info < (3, 4):
    raise OSError("Python 3.4 or greater is required")

import logging
DEFAULT_LOGGING_LEVEL = logging.INFO

__all__ = ["parser", "calculation"]
