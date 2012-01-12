import unittest
from cdo import *

class CdoTest(unittest.TestCase):
    def testCDO(self):
        cdo = Cdo()
        print(cdo.CDO)

    def testOps(self):
        cdo = Cdo()
        self.assertIn("sinfov",cdo.operators)
        self.assertIn("for",cdo.operators)
        self.assertIn("mask",cdo.operators)
        self.assertIn("studentt",cdo.operators)

    def testCall(self):
        cdo = Cdo()
        print cdo.sinfov(input='/home/ram/data/icon/oce.nc')

    def test_simple(self):
      cdo = Cdo()
      s = cdo.sinfov(input="-topo",options="-f nc")
      s = cdo.sinfov(input="-remapnn,r36x18 -topo",options="-f nc")
      f = 'ofile.nc'
      cdo.expr("'z=log(abs(topo+1))*9.81'",input="-topo",output = f,options="-f nc")
      s = cdo.infov(input=f)
      cdo.stdatm("0",output=f,options="-f nc")

    def test_returnArray(self):
      cdo = Cdo()
      cdo.setReturnArray()
      outfile = 'test.nc'
      press = cdo.stdatm("0",output=outfile,options="-f nc").var("P").get()
      self.assertEqual(1013.25,press.min())
      cdo.unsetReturnArray()
      press = cdo.stdatm("0",output=outfile,options="-f nc")
      self.assertEqual(press,outfile)
      press = cdo.stdatm("0",output=outfile,options="-f nc",returnArray=True).var("P").get()
      self.assertEqual(1013.25,press.min())
      print("press = "+press.min().__str__())


    def test_combine(self):
      cdo = Cdo()
      stdatm = cdo.stdatm("0",options = "-f nc",output="F")
      sum = cdo.fldsum(input = stdatm,output = "G")
      sum = cdo.fldsum(input = cdo.stdatm("0",options="-f nc",output="F"),output="G")
      sum = cdo.fldsum(input = cdo.stdatm("0",options="-f nc",output="F"),output="G",returnArray=True)

    def test_cdf(self):
      cdo = Cdo()
      self.assertNotIn("cdf",cdo.__dict__)
      cdo.setReturnArray()
      self.assertIn("cdf",cdo.__dict__)
      cdo.setReturnArray(False)
      sum = cdo.fldsum(input = cdo.stdatm("0",options="-f nc",output="F"),output="G",returnArray=True)
      self.assertEqual(1013.25,sum.var("P").get().min())

if __name__ == '__main__':
    unittest.main()
