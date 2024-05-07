import unittest
from unittest.mock import patch
from test.mock_database import mock_data, mock_data_mixture
import copy
from HGSpy.hgs import HGSData
from HGSpy.hgs_isentropic import hastobeS_frozen, hastobeM_frozen, hastobeS_shifting, \
    hastobeM_shifting, hgs_isentropic, hgs_isentropic_ids
from HGSpy.hgs import HGSData
from test.mock_database import mock_data


class HGSisentropicTestCase(unittest.TestCase):

    def setUp(self):
        self.data = mock_data

    def test_hastobeS_shifting(self):
        with patch('HGSpy.hgs_isentropic.hgs_prop_ids') as mock_prop, \
                patch('HGSpy.hgs_isentropic.hgs_eq_ids') as mock_eq:
            mock_eq.side_effect = [[[1, 2], 2]]
            _ = hastobeS_shifting(300, 1000, 12, [1, 1], [0, 1], {},
                                  HGSData().new(copy.deepcopy(self.data)))

    @patch('HGSpy.hgs_isentropic.hgs_prop_ids')
    def test_hastobeS_frozen(self, mock_hgs_prop):
        hastobeS_frozen(300, 1000, 12, [1, 1], [0, 1], {},
                        HGSData().new(copy.deepcopy(self.data)))

    @patch('HGSpy.hgs_isentropic.hgs_secant')
    @patch('HGSpy.hgs_isentropic.hgs_prop_ids')
    @patch('HGSpy.hgs_isentropic.hgs_eq_ids')
    def test_hastobeM_shifting(self, mock_hgs_eq, mock_hgs_prop, mock_secant):
        mock_hgs_eq.side_effect = [[[1, 2], 2]]
        mock_hgs_prop.side_effect = [[1, 2]]
        mock_secant.side_effect = [[[1,1], [1, 1], 1]]
        hastobeM_shifting(300, 100, 23, 1000, 12, [1, 1], [0, 1], {},
                          HGSData().new(copy.deepcopy(self.data)))

    @patch('HGSpy.hgs_isentropic.hgs_secant')
    @patch('HGSpy.hgs_isentropic.hgs_prop_ids')
    @patch('HGSpy.hgs_isentropic.hgs_eq_ids')
    @patch('builtins.print')
    def test_hastobeM_shifting_flag(self, mock_print, mock_hgs_eq, mock_hgs_prop, mock_secant):
        mock_hgs_eq.side_effect = [[[1, 2], 2]]
        mock_secant.side_effect = [[[1, 1], [1, 1], 2]]
        with self.assertRaises(SystemExit):
            hastobeM_shifting(300, 100, 23, 1000, 12, [1, 1], [0, 1], {},
                              HGSData().new(copy.deepcopy(self.data)))

    @patch('HGSpy.hgs_isentropic.hgs_secant')
    @patch('HGSpy.hgs_isentropic.hgs_prop_ids')
    def test_hastobeM_frozen(self, mock_hgs_prop, mock_secant):
        mock_hgs_prop.side_effect = [[1,2]]
        mock_secant.side_effect =[[[1,1],[1,1],1]]
        hastobeM_frozen(300, 100, 23, 1000, 12, [1, 1], [0, 1], {},
                        HGSData().new(copy.deepcopy(self.data)))

    @patch('HGSpy.hgs_isentropic.hgs_secant')
    @patch('HGSpy.hgs_isentropic.hgs_prop_ids')
    @patch('builtins.print')
    def test_hastobeM_frozen_flag(self,  mock_print, mock_hgs_prop, mock_secant):
        mock_secant.side_effect = [[[1, 1], [1, 1], 2]]
        with self.assertRaises(SystemExit):
            hastobeM_frozen(300, 100, 23, 1000, 12, [1, 1], [0, 1], {},
                            HGSData().new(copy.deepcopy(self.data)))

    @patch('HGSpy.hgs_isentropic.hgs_prop_ids')
    @patch('HGSpy.hgs_isentropic.hgs_solver')
    def test_hgs_isentropic_ids(self, mock_hgs_solver, mock_hgs_prop):
        # Mock dependencies
        mock_hgs_solver.side_effect = [[1, 1, 1]]
        mock_hgs_prop.side_effect = [[12, 1, 3], [12, 1, 3]]
        # Test with mocked dependencies
        hgs_isentropic_ids([1, 2, 3], [0.1, 0.2, 0.3], [100, 100, 100], 1, 'P', 1000,
                           flow='shifting',
                           solver='hgs_secant',
                           Tstar=3000, opt_eq={}, opt_sci={},
                           opt_sec={'xmin': 300, 'xmax': 4000, 'maxiter': 200,
                                    'epsx': 0.1, 'epsy': 1, 'fchange': 5,
                                    'tipo': 'Shifting', 'info': 0, 'dT': 100},
                           hgs_data=HGSData().new(**copy.deepcopy(mock_data)))

    @patch('HGSpy.hgs_isentropic.hgs_secant')
    @patch('HGSpy.hgs_isentropic.hgs_prop_ids')
    @patch('HGSpy.hgs_isentropic.hgs_solver')
    def test_hgs_isentropic_ids_M(self, mock_hgs_solver, mock_hgs_prop, mock_hgs_secant):
        # Mock dependencies
        mock_hgs_secant.side_effect = [[1, [1,1,1], 1]]
        mock_hgs_solver.side_effect = [[1, 1, 1]]
        mock_hgs_prop.side_effect = [[12, 1, 3], [12, 1, 3]]
        # Test with mocked dependencies
        hgs_isentropic_ids([1, 2, 3], [0.1, 0.2, 0.3], [100, 100, 100], 1, 'M', 1000,
                           flow='shifting',
                           solver='hgs_secant',
                           Tstar=3000, opt_eq={}, opt_sci={},
                           opt_sec={'xmin': 300, 'xmax': 4000, 'maxiter': 200,
                                    'epsx': 0.1, 'epsy': 1, 'fchange': 5,
                                    'tipo': 'Shifting', 'info': 0, 'dTp': 100},
                           hgs_data=HGSData().new(**copy.deepcopy(mock_data)))

    @patch('HGSpy.hgs_isentropic.hgs_isentropic_ids')
    def test_hgs_isentropic(self, mock_hgs_isentropic_ids):
        mock_hgsdata = HGSData().new(**copy.deepcopy(mock_data))
        mock_hgs_isentropic_ids.side_effect = [[1, 1, 1, 1, 1]]
        # Test with mocked dependencies
        hgs_isentropic(['species1', 'species2'], [2, 1], 1000, 1, 'P', 100, hgs_data=mock_hgsdata)

    @patch('HGSpy.hgs_isentropic.hgs_isentropic_ids')
    def test_hgs_isentropic_other(self, mock_hgs_isentropic_ids):
        mock_hgsdata = HGSData().new(**copy.deepcopy(mock_data))
        mock_hgs_isentropic_ids.side_effect = [[1, 1, 1, 1, 1]]
        # Test with mocked dependencies
        hgs_isentropic(['species1', 'species2'], [2, 1], [1000, 1000], 1, 'P', 100,
                       hgs_data=mock_hgsdata)

    @patch('HGSpy.hgs_isentropic.hgs_isentropic_ids')
    def test_hgs_isentropic_mixture(self, mock_hgs_isentropic_ids):
        mock_hgsdata = HGSData().new(**copy.deepcopy(mock_data_mixture))

        mock_hgs_isentropic_ids.side_effect = [[1, 1, 1, 1, 1]]
        # Test with mocked dependencies
        hgs_isentropic(['mixture1'], [1], 1000, 1, 'P', 100, hgs_data=mock_hgsdata)

    @patch('HGSpy.hgs_isentropic.hgs_isentropic_ids')
    def test_hgs_isentropic_H(self, mock_hgs_isentropic_ids):
        mock_hgsdata = HGSData().new(**copy.deepcopy(mock_data_mixture))

        mock_hgs_isentropic_ids.side_effect = [[1, 1, 1, 1, 1]]
        # Test with mocked dependencies
        hgs_isentropic(['mixture1'], [1], 1000, 1, 'M', 100, hgs_data=mock_hgsdata)

    @patch('HGSpy.hgs_isentropic.hgs_isentropic_ids')
    def test_hgs_isentropic_mr(self, mock_hgs_isentropic_ids):
        mock_hgsdata = HGSData().new(**copy.deepcopy(mock_data_mixture))

        mock_hgs_isentropic_ids.side_effect = [[1, 1, 1, 1, 1]]
        # Test with mocked dependencies
        hgs_isentropic(['mixture1'], [1], [1000], 1, 'M', 100, hgs_data=mock_hgsdata)



if __name__ == '__main__':
    unittest.main()
