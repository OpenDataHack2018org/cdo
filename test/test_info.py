import unittest
from cdo import *


class TestInfo(unittest.TestCase):

    def test_info(self):
        cdo = Cdo()
        levels = cdo.showlevel(input = "-stdatm,0")
        info   = cdo.sinfo(input = "-stdatm,0")
        print(levels)
        self.assertEqual([0,0],map(float,levels))
        self.assertEqual("File format: GRIB",info[0])


if __name__ == '__main__':
    unittest.main()
