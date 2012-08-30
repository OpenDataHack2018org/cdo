import unittest
from cdo import *


class TestInfo(unittest.TestCase):

    # check operators: (s)infov,showlevel
    def test_info(self):
        cdo  = Cdo()

        info = cdo.sinfo(input = "-stdatm,0")
        self.assertEqual("File format: GRIB",info[0])

        sinfov = cdo.sinfov(input = "-stdatm,0,10,20,50,100,500,1000",options = "-f nc")
        infov  = cdo.infov( input = "-stdatm,0,10,20,50,100,500,1000",options = "-f nc")
        self.assertEqual("P",sinfov[2].split(' ')[6])
        self.assertEqual("T",sinfov[3].split(' ')[6])
        self.assertEqual("1013.2", infov[1].split(' ')[-1])
        self.assertEqual("T", infov[-1].split(' ')[4])

        levels = cdo.showlevel(input = "-stdatm,0")
        self.assertEqual([0,0],map(float,levels))
        levels = cdo.showlevel(input = "-stdatm,10,20,30")[0].split(' ')
        self.assertEqual([10,20,30],map(float,levels))

if __name__ == '__main__':
    unittest.main()
