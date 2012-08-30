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

        sinfov = cdo.sinfov(input = "-stdatm,0,10,20,50,100,500,1000",options = "-f nc")
        self.assertEqual("P",sinfov[2].split(' ')[6])
        self.assertEqual("T",sinfov[3].split(' ')[6])

if __name__ == '__main__':
    unittest.main()
