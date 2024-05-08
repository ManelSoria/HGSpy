import unittest
from unittest.mock import patch

import numpy as np

from HGSpy.hgs_secant import hgs_secant


class HGSSecantTestCase(unittest.TestCase):

    @patch('HGSpy.hgs_secant.cr_start')
    @patch('HGSpy.hgs_secant.cr_stop')
    def test_hgs_secant(self, mock_crstop, mock_crstart):
        def fun(x,l):
            return x -20, l
        opt_sec = {
            'xmin': 0,
            'xmax': 30,
            'maxiter': 1000,
            'epsx': 0.1,
            'epsy': 0.001,
            'fchange': 1,
            'info': 0,
            'dTp': 1

        }
        res, *_ = hgs_secant(fun, [3], opt_sec, ilevel=0)

        self.assertEqual(res, 20)

    @patch('HGSpy.hgs_secant.cr_start')
    @patch('HGSpy.hgs_secant.cr_stop')
    def test_hgs_secant_squared(self, mock_crstop, mock_crstart):
        def fun(x, l):
            return 0.0000001*x**2+ 2*x- 90, l

        opt_sec = {
            'xmin': -500,
            'xmax': 1500,
            'maxiter': 10000,
            'epsx': 0.1,
            'epsy': 0.001,
            'fchange': 1,
            'info': 0,
            'dTp': 1

        }
        res, *_ = hgs_secant(fun, [3], opt_sec, ilevel=0)

        self.assertAlmostEqual(res, 44.9999, places=4)

    @patch('HGSpy.hgs_secant.cr_start')
    @patch('HGSpy.hgs_secant.cr_stop')
    def test_hgs_secant_squared_2(self, mock_crstop, mock_crstart):
        def fun(x, l):
            return -0.0000001 * x ** 2 - 2 * x + 90, l

        opt_sec = {
            'xmin': -500,
            'xmax': 1500,
            'maxiter': 10000,
            'epsx': 0.1,
            'epsy': 0.001,
            'fchange': 1,
            'info': 0,
            'dTp': 1

        }
        res, *_ = hgs_secant(fun, [3], opt_sec, ilevel=0)

        self.assertAlmostEqual(res, 44.9999, places=4)

    @patch('HGSpy.hgs_secant.cr_start')
    @patch('HGSpy.hgs_secant.cr_stop')
    def test_hgs_secant_no_sign_change(self, mock_crstop, mock_crstart):
        def fun(x, l):
            return x - 20, l

        opt_sec = {
            'xmin': 25,
            'xmax': 30,
            'maxiter': 1000,
            'epsx': 0.1,
            'epsy': 0.001,
            'fchange': 1,
            'info': 0,
            'dTp': 1

        }

        _, _, flag = hgs_secant(fun, [3], opt_sec, ilevel=0)
        self.assertEqual(flag, -2)

        @patch('HGSpy.hgs_secant.cr_start')
        @patch('HGSpy.hgs_secant.cr_stop')
        def test_hgs_secant_out_of_iter(self, mock_crstop, mock_crstart):
            def fun(x, l):
                return x - 20, l

            opt_sec = {
                'xmin': 0,
                'xmax': 40000,
                'maxiter': 10,
                'epsx': 0.1,
                'epsy': 0.001,
                'fchange': 1,
                'info': 0,
                'dTp': 1

            }

            _, _, flag = hgs_secant(fun, [3], opt_sec, ilevel=0)
            self.assertEqual(flag, -1)

if __name__ == '__main__':
    unittest.main()
