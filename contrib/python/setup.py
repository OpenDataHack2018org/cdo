#!/usr/bin/env python
from distutils.core import setup, Extension
from setup import *

setup (name   = 'cdo',
  version     = '1.0.2',
  author      = "Ralf Mueller",
  author_email= "stark.dreamdetective@gmail.com",
  license     = "GPLv2",
  description = """pyhton bindings to CDO""",
  py_modules  = ["cdo"],
  )
