import copy
import unittest
from unittest.mock import patch
from HGSpy.hgs_solver import options, hastobezeroH_frozen, hastobezeroS_frozen, \
    hastobezeroH_shifting, hastobezeroS_shifting, hgs_solver
from HGSpy.hgs import HGSData
from test.mock_database import mock_data


class HGSSolveTestCase(unittest.TestCase):

    def setUp(self):
        self.data = mock_data

    def test_hastobezeroH_shifting(self):
        with patch('HGSpy.hgs_solver.hgs_prop_ids') as mock_prop, \
              patch('HGSpy.hgs_solver.hgs_eq_ids') as mock_eq:
            mock_eq.side_effect = [[1,2]]
            _ = hastobezeroH_shifting(300, 12, [1, 1], [0, 1], 0, options,
                                  HGSData().new(copy.deepcopy(self.data)))

    @patch('HGSpy.hgs_solver.hgs_prop_ids')
    def test_hastobezeroH_frozen(self, mock_hgs_prop):
        hastobezeroH_frozen(300, 12, [1, 1], [0, 1], 0, options,
                            HGSData().new(copy.deepcopy(self.data)))

    @patch('HGSpy.hgs_solver.hgs_prop_ids')
    @patch('HGSpy.hgs_solver.hgs_eq_ids')
    def test_hastobezeroS_shifting(self, mock_hgs_eq, mock_hgs_prop):
        mock_hgs_eq.side_effect = [[1, 2]]
        hastobezeroS_shifting(300, 12, [1, 1], [0, 1], 0, options,
                              HGSData().new(copy.deepcopy(self.data)))

    @patch('HGSpy.hgs_solver.hgs_prop_ids')
    def test_hastobezeroS_frozen(self, mock_hgs_prop):
        hastobezeroS_frozen(300, 12, [1, 1], [0, 1], 0, options,
                            HGSData().new(copy.deepcopy(self.data)))

    @patch('HGSpy.hgs_solver.hgs_secant')
    def test_hgs_solver_H(self, mock_hgs_secant):
        mock_hgs_secant.side_effect = [[1,2,3]]
        hgs_solver(
            [0, 1],
            [1, 1],
            'H',
            0,
            12,
            flow='shifting',
            solver='hgs_secant',
            Tstar=3000,
            opt_sec=options,
            hgs_data=HGSData().new(copy.deepcopy(self.data))
        )
    @patch('HGSpy.hgs_solver.hgs_secant')
    def test_hgs_solver_S(self, mock_hgs_secant):
        mock_hgs_secant.side_effect = [[1,2,3]]
        hgs_solver(
            [0, 1],
            [1, 1],
            'H',
            0,
            12,
            flow='shifting',
            solver='hgs_secant',
            Tstar=3000,
            opt_sec=options,
            hgs_data=HGSData().new(copy.deepcopy(self.data))
        )

    @patch('HGSpy.hgs_solver.hgs_secant')
    def test_hgs_solver_frozen(self, mock_hgs_secant):
        mock_hgs_secant.side_effect = [[1, 2, 3]]
        hgs_solver(
            [0, 1],
            [1, 1],
            'H',
            0,
            12,
            flow='frozen',
            solver='hgs_secant',
            Tstar=3000,
            opt_sec=options,
            hgs_data=HGSData().new(copy.deepcopy(self.data))
        )

    @patch('builtins.print')
    def test_hgs_solver_no_prop(self, mock_print):
        with self.assertRaises(SystemExit):
            hgs_solver(
                [0, 1],
                [1, 1],
                'G',
                0,
                12,
                flow='frozen',
                solver='hgs_secant',
                Tstar=3000,
                opt_sec=options,
                hgs_data=HGSData().new(copy.deepcopy(self.data))
            )



if __name__ == '__main__':
    unittest.main()
