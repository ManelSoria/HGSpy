import unittest
from unittest.mock import patch
from HGSpy.utils import set_options, get_options, raiseError, raiseWarning, truncate


class UtilsTestCase(unittest.TestCase):
    def test_set_options(self):
        set_options('warnings', False)
        self.assertFalse(get_options('warnings'))

    def test_get_options(self):
        set_options('errors', False)
        self.assertFalse(get_options('errors'))

    @patch('builtins.print')
    def test_raiseError(self, mock_print):
        # NOTE: The values is set to true due to  instance issues
        set_options('errors', True)
        with self.assertRaises(SystemExit):
            raiseError("Test error message")


    def test_truncate(self):
        self.assertEqual(truncate(3.14159, 2), 3.14)


if __name__ == '__main__':
    unittest.main()
