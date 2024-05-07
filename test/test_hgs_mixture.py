import unittest
import copy
from unittest.mock import patch
from test.mock_database import mock_data
from HGSpy.hgs_mixture import hgs_add_mixture
from HGSpy.hgs import HGSData

class TestAddMixtureData(unittest.TestCase):


    @patch('builtins.print')
    @patch('HGSpy.HGSData.id')
    def test_hgs_add_mixture_already_exist(self, mock_print, mock_hgs_id):
        test_data = HGSData().new(**copy.deepcopy(mock_data))
        mock_hgs_id.return_value = [1]
        with self.assertRaises(SystemExit):
            hgs_add_mixture("mock_mixture", ["species1", "species2"], [10,90],test_data)


    @patch('builtins.print')
    @patch('HGSpy.HGSData.id')
    def test_hgs_add_mixture_no_exist_species(self, mock_print, mock_hgs_id):
        test_data = HGSData().new(**copy.deepcopy(mock_data))
        mock_hgs_id.side_effect = [[None], [None,1]]
        with self.assertRaises(SystemExit):
            hgs_add_mixture("mock_mixture", ["species3", "species2"], [10, 90], test_data)

    @patch('builtins.print')
    @patch('HGSpy.HGSData.id')
    def test_hgs_add_mixture_no_exist_species(self, mock_print, mock_hgs_id):
        test_data = HGSData().new(**copy.deepcopy(mock_data))
        mock_hgs_id.side_effect = [[None], [0, 1]]
        with self.assertRaises(SystemExit):
            hgs_add_mixture("mock_mixture", ["species1", "species2"], [10, 90], test_data)


if __name__ == '__main__':
    unittest.main()
