import unittest
from cdo import *


class TestInfo(unittest.TestCase):

    # check operators: (s)infov,showlevel
    def test_info(self):
        cdo  = Cdo()

        info = cdo.sinfo(input = "-stdatm,0")
        self.assertEqual("File format: GRIB",info[0])

        sinfov = cdo.sinfov(input = "-stdatm,0,10,20,50,100,500,1000",options = "-f nc")
        self.assertEqual("P",sinfov[2].split(' ')[-1])
        self.assertEqual("T",sinfov[3].split(' ')[-1])

        infov  = cdo.infov( input = "-stdatm,0,10,20,50,100,500,1000",options = "-f nc")
        self.assertEqual("1013.2", infov[1].split(' ')[-1])
        self.assertEqual("T", infov[-1].split(' ')[4])

        units  = cdo.showunit( input = "-stdatm,0", options = "-f nc")
        self.assertEqual(["hPa","K"],units[0].split(' '))

        levels = cdo.showlevel(input = "-stdatm,0")
        self.assertEqual([0,0],map(float,levels))
        levels = cdo.showlevel(input = "-stdatm,10,20,30")[0].split(' ')
        self.assertEqual([10,20,30],map(float,levels))

        code = cdo.showcode(input="-topo,r1x1")
        self.assertEqual("-1",code[0])
        code = cdo.showcode(input="-setcode,111 -topo,r1x1")
        self.assertEqual("111",code[0])

if __name__ == '__main__':
    unittest.main()
