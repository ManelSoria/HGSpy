import copy
import unittest
from test.mock_database import mock_data_mixture
from HGSpy.hgs import HGSData
from HGSpy.hgs_find import hgs_find
from io import StringIO
import sys


class HGSFindTestCase(unittest.TestCase):
    def test_hgs_find_no_complete(self):
        data = HGSData().new(**copy.deepcopy(mock_data_mixture))
        captured_output = StringIO()

        sys.stdout = captured_output

        hgs_find('species', False, data)
        captured_output.seek(0)
        output = captured_output.getvalue()

        self.assertIn("Species that contain species\n"
                      "<0>  species1\n"
                      "<1>  species2\n"
                      "<2>  species3\n"
                      "Mixtures that contain species\n"
                      "None", output)

    def test_hgs_find_complete(self):
        data = HGSData().new(**copy.deepcopy(mock_data_mixture))
        captured_output = StringIO()

        sys.stdout = captured_output

        hgs_find('species', True, data)
        captured_output.seek(0)
        output = captured_output.getvalue()

        self.assertIn("Species that contain species\n"
                      "None\n"
                      "Mixtures that contain species\n"
                      "None", output)

    def test_hgs_find_complete_mixture(self):
        data = HGSData().new(**copy.deepcopy(mock_data_mixture))
        captured_output = StringIO()

        sys.stdout = captured_output

        hgs_find('mixture', True, data)
        captured_output.seek(0)
        output = captured_output.getvalue()

        self.assertIn("Species that contain mixture\n"
                      "None\n"
                      "Mixtures that contain mixture\n"
                      "<3>  mixture1", output)

    def test_hgs_find_another(self):
        data = HGSData().new(**copy.deepcopy(mock_data_mixture))
        captured_output = StringIO()

        sys.stdout = captured_output

        hgs_find('1', False, data)
        captured_output.seek(0)
        output = captured_output.getvalue()

        self.assertIn("Species that contain 1\n"
                      "<0>  species1\n"
                      "Mixtures that contain 1\n"
                      "<3>  mixture1", output)

    def test_hgs_find_complete_true(self):
        data = HGSData().new(**copy.deepcopy(mock_data_mixture))
        captured_output = StringIO()

        sys.stdout = captured_output

        hgs_find('species2', False, data)
        captured_output.seek(0)
        output = captured_output.getvalue()

        self.assertIn("Species that contain species2\n"
                      "<1>  species2\n"
                      "Mixtures that contain species2\n"
                      "None", output)

if __name__ == '__main__':
    unittest.main()
