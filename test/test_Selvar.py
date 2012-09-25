#!/usr/bin/env python
import unittest
from cdo import *
import testStreams

cdo = Cdo()
testFiles = testStreams.createInputFiles(4,'r36x18',cdo,'-f nc')


class TestSelvar(unittest.TestCase):

  def test_selparam(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.selparam('130',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_delparam(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.delparam('130',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_selcode(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.selcode('130',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_delcode(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.delcode('130',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_selname(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.selname('T',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_delname(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.delname('T',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_selstdname(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.selstdname('T',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_sellevel(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.sellevel('0,10,20,50,100,200,500,1000,2000,5000,10000',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_sellevidx(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.sellevidx('1',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_selgrid(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.selgrid('1',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_selzaxis(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.selzaxis('1',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_selltype(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.selltype('105',input = ' '.join(files[0:1]))
      self.assertEqual("","")


  def test_seltabnum(self):
    cdo = Cdo()
    for i,value in enumerate(testFiles):
      files = testFiles[value]
      result = cdo.seltabnum('128',input = ' '.join(files[0:1]))
      self.assertEqual("","")


if __name__ == '__main__':
  unittest.main()
