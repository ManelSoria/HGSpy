import unittest
from unittest.mock import patch
from test.mock_database import mock_data, mock_data_mixture
import copy
from HGSpy.hgs import HGSData
from HGSpy.hgs_Tp import hgs_Tp, hgs_Tp_ids


class HGSTpTestCase(unittest.TestCase):
    @patch('HGSpy.hgs_Tp.hgs_prop_ids')
    @patch('HGSpy.hgs_Tp.hgs_solver')
    def test_hgs_Tp_ids(self, mock_hgs_solver, mock_hgs_prop):
        # Mock dependencies
        mock_hgs_solver.side_effect = [[1, 1, 1]]
        # Test with mocked dependencies
        hgs_Tp_ids([1, 2, 3], [0.1, 0.2, 0.3], 'H', [100, 100, 100], 1, flow='shifting',
                   solver='hgs_secant',
                   Tstar=3000, opt_eq={}, opt_sci={},
                   opt_sec={'xmin': 300, 'xmax': 4000, 'maxiter': 200,
                            'epsx': 0.1, 'epsy': 1, 'fchange': 5,
                            'tipo': 'Shifting', 'info': 0, 'dTp': 100},
                   hgs_data=HGSData().new(**copy.deepcopy(mock_data)))

    @patch('HGSpy.hgs_Tp.hgs_Tp_ids')
    def test_hgs_Tp(self, mock_hgs_Tp_ids):
        mock_hgsdata = HGSData().new(**copy.deepcopy(mock_data))
        mock_hgs_Tp_ids.side_effect = [[1, 1, 1]]
        # Test with mocked dependencies
        hgs_Tp(['species1', 'species2'], [2, 1], 'T', 1000, 1, hgs_data=mock_hgsdata)

    @patch('HGSpy.hgs_Tp.hgs_Tp_ids')
    def test_hgs_Tp_other(self, mock_hgs_Tp_ids):
        mock_hgsdata = HGSData().new(**copy.deepcopy(mock_data))
        mock_hgs_Tp_ids.side_effect = [[1, 1, 1]]
        # Test with mocked dependencies
        hgs_Tp(['species1', 'species2'], [2, 1], 'T', [1000, 1000], 1, hgs_data=mock_hgsdata)

    @patch('HGSpy.hgs_Tp.hgs_Tp_ids')
    def test_hgs_Tp_mixture(self, mock_hgs_Tp_ids):
        mock_hgsdata = HGSData().new(**copy.deepcopy(mock_data_mixture))

        mock_hgs_Tp_ids.side_effect = [[1, 1, 1]]
        # Test with mocked dependencies
        hgs_Tp(['mixture1'], [1], 'T', 1000, 1, hgs_data=mock_hgsdata)

    @patch('HGSpy.hgs_Tp.hgs_Tp_ids')
    def test_hgs_Tp_H(self, mock_hgs_Tp_ids):
        mock_hgsdata = HGSData().new(**copy.deepcopy(mock_data_mixture))

        mock_hgs_Tp_ids.side_effect = [[1, 1, 1]]
        # Test with mocked dependencies
        hgs_Tp(['mixture1'], [1], 'H', 1000, 1, hgs_data=mock_hgsdata)

    @patch('HGSpy.hgs_Tp.hgs_Tp_ids')
    def test_hgs_Tp_H_other(self, mock_hgs_Tp_ids):
        mock_hgsdata = HGSData().new(**copy.deepcopy(mock_data_mixture))

        mock_hgs_Tp_ids.side_effect = [[1, 1, 1]]
        # Test with mocked dependencies
        hgs_Tp(['mixture1'], [1], 'H', [1000], 1, hgs_data=mock_hgsdata)

    @patch('builtins.print')
    def test_hgs_Tperror(self, mock_print):
        mock_hgsdata = HGSData().new(**copy.deepcopy(mock_data_mixture))

        with self.assertRaises(SystemExit):
            # Test with mocked dependencies
            hgs_Tp(['mixture1'], [1], 'L', 1000, 1, hgs_data=mock_hgsdata)


if __name__ == '__main__':
    unittest.main()
