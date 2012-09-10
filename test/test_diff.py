#!/usr/bin/env python
import unittest
from cdo import *


class TestDiff(unittest.TestCase):

    def test_diff(self):
        cdo = Cdo()
        diffv = cdo.diffn(input = "-random,r1x1 -random,r1x1")
        self.assertEqual(diffv[1].split(' ')[-1],"random")
        self.assertEqual(diffv[1].split(' ')[-3],"0.53060")
        diff  = cdo.diff(input = "-random,r1x1 -random,r1x1")
        self.assertEqual(diff[1].split(' ')[-3],"0.53060")


if __name__ == '__main__':
    unittest.main()
