#!/usr/bin/env python
import unittest
from cdo import *


class TestDiff(unittest.TestCase):

    def test_diff(self):
        cdo = Cdo()
        diffv = cdo.diffn(input = "-const,4,r1x1 -const,2,r1x1")
        self.assertEqual(diffv[1].split(' ')[-1],"const")
        self.assertEqual(diffv[1].split(' ')[-3],"0.50000")
        diff  = cdo.diff(input = "-const,4,r1x1 -const,2,r1x1")
        self.assertEqual(diff[1].split(' ')[-3],"0.50000")


if __name__ == '__main__':
    unittest.main()
