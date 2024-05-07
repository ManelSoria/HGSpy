import unittest
from unittest.mock import patch
from test.mock_database import mock_data, mock_data_mixture
from HGSpy.hgs_id import hgs_id, find_hgs_id
from HGSpy.hgs import HGSData
import copy

class HGSIdTestCase(unittest.TestCase):


    @patch('HGSpy.hgs_id.find_hgs_id')
    def test_hgs_id(self, mock_find_hgs_id):
        data = HGSData.new(**copy.deepcopy(mock_data))
        _ = hgs_id(['species1', 'species2'], data, raise_error=True)
        mock_find_hgs_id.assert_called()
        self.assertEqual(mock_find_hgs_id.call_count, 2)

    @patch('HGSpy.hgs_id.find_hgs_id')
    def test_hgs_id_str(self, mock_find_hgs_id):
        data = HGSData.new(**copy.deepcopy(mock_data))
        _ = hgs_id('species1', data, raise_error=True)
        mock_find_hgs_id.assert_called()
        self.assertEqual(mock_find_hgs_id.call_count, 1)

    @patch('HGSpy.hgs_id.find_hgs_id')
    def test_hgs_id_mixture(self, mock_find_hgs_id):
        data = HGSData.new(**copy.deepcopy(mock_data_mixture), )
        _ = hgs_id(['mixture1'], data, raise_error=True)
        mock_find_hgs_id.assert_called()
        self.assertEqual(mock_find_hgs_id.call_count, 1)

    @patch('HGSpy.cr.cr_stop')
    @patch('builtins.print')
    def test_hgs_id_fail(self, mock_print, mock_cr):
        data = HGSData.new(**copy.deepcopy(mock_data_mixture))
        with self.assertRaises(SystemExit):
            _ = hgs_id([1], data, raise_error=True)

    def test_find_hgsdata(self):
        data = HGSData.new(**copy.deepcopy(mock_data_mixture))
        tmp = find_hgs_id('species1', data, True,[], True)
        self.assertEqual(tmp, 0)
        tmp = find_hgs_id('mixture1', data, True, [], True)
        self.assertEqual(tmp, 3)

    @patch('builtins.print')
    def test_find_hgsdata_2(self, mock_print):
        data = HGSData.new(**copy.deepcopy(mock_data))
        tmp = find_hgs_id('species2', data, False,[], True)
        self.assertEqual(tmp, 1)

    @patch('builtins.print')
    def test_find_hgsdata_on_fail(self, mock_print):
        data = HGSData.new(**copy.deepcopy(mock_data_mixture))
        with self.assertRaises(SystemExit):
            find_hgs_id('mixture4', data, True, [], True)

        tmp = find_hgs_id('mixture4', data, True, [], False)
        self.assertIsNone(tmp)

if __name__ == '__main__':
    unittest.main()
