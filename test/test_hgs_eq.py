import unittest
from unittest.mock import patch, MagicMock
from HGSpy.hgs_eq import parameters_min, hgs_eq_ids, hgs_eq, options
from HGSpy.hgs import HGSData
from test.mock_database import mock_data
import numpy as np
import copy


class HGSeqTestCase(unittest.TestCase):
    def test_parameters_min(self):
        data = HGSData().new(**mock_data)
        bound, dictionary = parameters_min([0,1],[1,1], data)
        self.assertTrue(np.array_equal(bound.lb,[0]*2))
        self.assertTrue(np.array_equal(bound.ub, [np.inf] * 2))

        func = dictionary['fun']

        self.assertTrue(np.array_equal(func([3,2]), [7,3]))
        self.assertTrue(np.array_equal(func([1, 2]), [3, 1]))

    def test_hgs_eq_ids(self):

        class TmpClass:
            x =  MagicMock
            success =  True
            fun = MagicMock

        data = HGSData().new(**mock_data)

        with patch('HGSpy.hgs_eq.minimize', return_value=TmpClass()) as mock_minimize:
            x, fun = hgs_eq_ids([0, 1], [1, 2], [300] * 2, 12, options, data)


    def test_hgs_eq(self):

        class TmpClass:
            x =  MagicMock
            success =  True
            fun = MagicMock

        data = HGSData().new(**copy.deepcopy(mock_data))

        with patch('HGSpy.hgs_eq.minimize', return_value=TmpClass()) as mock_minimize:
             hgs_eq(['species1', 'species2', 'species3'], [1, 2, 3], 3000, 12, options, data)

        with patch('HGSpy.hgs_eq.minimize', return_value=TmpClass()) as mock_minimize:
            hgs_eq(['species1', 'species2', 'species3'], [1, 2, 3], [3000], 12, options, data)

        with patch('HGSpy.hgs_eq.minimize', return_value=TmpClass()) as mock_minimize:
            hgs_eq(['species1', 'species2', 'species3'], [1, 2, 3], [3000]*3, 12, options, data)

    @patch('builtins.print')
    def test_hgs_eq_fail(self, mock_print):

        class TmpClass:
            x =  MagicMock
            success =  True
            fun = MagicMock

        data = HGSData().new(**copy.deepcopy(mock_data))

        with self.assertRaises(SystemExit),  patch('HGSpy.hgs_eq.minimize', return_value=TmpClass()) as mock_minimize:
            x, fun = hgs_eq([0, 1], [1, 2, 3], 3000, 12, options, data)



if __name__ == '__main__':
    unittest.main()
