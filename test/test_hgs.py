import unittest
from unittest.mock import patch
from HGSpy.hgs import HGSData
from io import StringIO
import sys
import copy

from test.mock_database import mock_data
class TestHGSData(unittest.TestCase):

    def setUp(self):
        # Create a sample data dictionary
        self.sample_data = copy.deepcopy(mock_data)

    def test_init(self):
        hgs_data = HGSData(self.sample_data)
        self.assertEqual(len(hgs_data), 3)

    def test_len(self):
        hgs_data = HGSData(self.sample_data)
        self.assertEqual(len(hgs_data), 3)

    def test_str(self):
        hgs_data = HGSData(self.sample_data)
        self.assertEqual(str(hgs_data), str(self.sample_data))

    def test_getitem(self):
        hgs_data = HGSData(self.sample_data)
        self.assertEqual(hgs_data['name'], ['species1', 'species2', 'species3'])

    def test_setitem(self):
        hgs_data = HGSData(self.sample_data)
        hgs_data['name'] = ['new_species1', 'new_species2']
        self.assertEqual(hgs_data['name'], ['new_species1', 'new_species2'])

    @patch('HGSpy.hgs_id.hgs_id')
    def test_id(self, mock_hgs_id):
        hgs_data = HGSData(self.sample_data)
        mock_hgs_id.return_value = [0, 1]
        self.assertEqual(hgs_data.id(['species1', 'species2']), [0, 1])

    @patch('HGSpy.hgs_add_mixture')
    @patch('HGSpy.hgs_subt_mixture')
    def test_add_subt_mixture(self, mock_hgs_subt_mixture , mock_hgs_add_mixture):
        hgs_data = HGSData(self.sample_data)

        if self.assertRaises(Exception,
            hgs_data.add('mixture_name', ['species1', 'species2'], [50, 50])):
            self.assertTrue(False)

        if self.assertRaises(Exception,
            hgs_data.subt('mixture_name')):
            self.assertTrue(False)


    @patch('HGSpy.hgs_rebuild')
    def test_rebuild(self, mock_hgs_rebuild):
        hgs_data = HGSData(self.sample_data)
        mock_hgs_rebuild.return_value = (['species1', 'species2'], [0.5, 0.5], [300, 300])
        if self.assertRaises(Exception,hgs_data.rebuild(['species1', 'species2'], [0.5, 0.5], [300, 300]),
                         (['species1', 'species2'], [0.5, 0.5], [300, 300])):
            self.assertTrue(False)

    @patch('HGSpy.hgs_print_info')
    def test_print_info(self, mock_hgs_print_info):

        captured_output = StringIO()

        sys.stdout = captured_output


        hgs_data = HGSData(self.sample_data)
        if self.assertRaises(Exception,hgs_data.print_info('species1')):
            self.assertTrue(False)
        captured_output.seek(0)
        output = captured_output.getvalue()

        self.assertIn("Species = <species1>   code = 0\n"
                      "- Composition:  atom1-2   atom2-1   \n"
                      "- Mm = 12.0000 \n"
                      " -----------------", output)


    def test_coefs(self):
        hgs_data = HGSData(self.sample_data)
        self.assertEqual(hgs_data.coefs(0, 500, all=True), ([1, 2, 3, 4, 5, 6,7], [7, 8, 9,10, 11, 12, 13]))

    def test_cp(self):
        hgs_data = HGSData(self.sample_data)
        self.assertAlmostEqual(hgs_data.cp(0, 500), 2602432881.4693513, places=4)

    def test_cv(self):
        hgs_data = HGSData(self.sample_data)
        self.assertAlmostEqual(hgs_data.cv(0, 500), 2602432881.4610367, places=4)

    def test_h(self):
        hgs_data = HGSData(self.sample_data)
        self.assertAlmostEqual(hgs_data.h(0, 500), 260347635896.8351, places=4)

    def test_s(self):
        hgs_data = HGSData(self.sample_data)
        self.assertAlmostEqual(hgs_data.s(0, 500, 1), 650956221.5934553, places=4)

    def test_g(self):
        hgs_data = HGSData(self.sample_data)
        self.assertAlmostEqual(hgs_data.g(0, 500, 1), -65130474899.89255, places=4)


if __name__ == '__main__':
    unittest.main()
