#!/usr/bin/env python
import unittest
from cdo import *
import testStreams

class TestArith(unittest.TestCase):

  def test_add(self):
    cdo = Cdo()
    cdo.debug = True
    testFiles = testStreams.createInputFiles(4,'r2x2',cdo,'-f nc')
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.add(input = ' '.join(files[0:2]))
      self.assertEqual("","")


  def test_sub(self):
    cdo = Cdo()
    cdo.debug = True
    testFiles = testStreams.createInputFiles(4,'r2x2',cdo,'-f nc')
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.sub(input = ' '.join(files[0:2]))
      self.assertEqual("","")


  def test_mul(self):
    cdo = Cdo()
    cdo.debug = True
    testFiles = testStreams.createInputFiles(4,'r2x2',cdo,'-f nc')
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.mul(input = ' '.join(files[0:2]))
      self.assertEqual("","")


  def test_div(self):
    cdo = Cdo()
    cdo.debug = True
    testFiles = testStreams.createInputFiles(4,'r2x2',cdo,'-f nc')
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.div(input = ' '.join(files[0:2]))
      self.assertEqual("","")


  def test_min(self):
    cdo = Cdo()
    cdo.debug = True
    testFiles = testStreams.createInputFiles(4,'r2x2',cdo,'-f nc')
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.min(input = ' '.join(files[0:2]))
      self.assertEqual("","")


  def test_max(self):
    cdo = Cdo()
    cdo.debug = True
    testFiles = testStreams.createInputFiles(4,'r2x2',cdo,'-f nc')
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.max(input = ' '.join(files[0:2]))
      self.assertEqual("","")


  def test_atan2(self):
    cdo = Cdo()
    cdo.debug = True
    testFiles = testStreams.createInputFiles(4,'r2x2',cdo,'-f nc')
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.atan2(input = ' '.join(files[0:2]))
      self.assertEqual("","")


if __name__ == '__main__':
  unittest.main()
