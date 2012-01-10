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


if __name__ == '__main__':
    unittest.main()
