import unittest
from cdo import *


class TestIntgridbil(unittest.TestCase):

    def test_intrgidbil(self):
        # dummy
        self.assertEqual(0,0)

    def test_global(self):
        # cdo cmd:
        # cdo -f nc intgridbil,global_.1 -topo result.nc
        cdo = Cdo()
        cdo.intgridbil('global_.1',input="-topo",output="result.nc",options="-f nc")
        # dummy, TODO cdo diff
        self.assertEqual(0,0)


if __name__ == '__main__':
    unittest.main()
