#!/usr/bin/env python
import unittest
from cdo import *
import testStreams

class TestArithc(unittest.TestCase):

  def test_addc(self):
    cdo = Cdo()
    testFiles = testStreams.createInputFiles(4,'r2x2',cdo,'-f nc')
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.addc('2.718281828459045',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_subc(self):
    cdo = Cdo()
    testFiles = testStreams.createInputFiles(4,'r2x2',cdo,'-f nc')
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.subc('2.718281828459045',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_mulc(self):
    cdo = Cdo()
    testFiles = testStreams.createInputFiles(4,'r2x2',cdo,'-f nc')
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.mulc('2.718281828459045',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_divc(self):
    cdo = Cdo()
    testFiles = testStreams.createInputFiles(4,'r2x2',cdo,'-f nc')
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.divc('2.718281828459045',input = ' '.join(files[0:1]))
      self.assertEqual("","")


if __name__ == '__main__':
  unittest.main()
