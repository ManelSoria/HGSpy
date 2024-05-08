import unittest
import copy
from unittest.mock import patch
from test.mock_database import mock_data
from HGSpy.hgs_prop import rg, gamma, sound, partial, hgs_prop_ids, hgs_single, hgs_prop
from parameterized import parameterized


class TestAddMixtureData(unittest.TestCase):

    def setUp(self):
        self.mock_database = copy.deepcopy(mock_data)

    def test_rg_float(self):
        self.assertAlmostEqual(rg(2.98), 2.7901, places=4)

    def test_rg_int(self):
        self.assertAlmostEqual(rg(10), 0.8314, places=4)

    def test_gamma_float(self):
        self.assertAlmostEqual(gamma(2.98, 1.1), 2.7091, places=4)

    def test_gamma_int(self):
        self.assertAlmostEqual(gamma(10, 1), 10.0, places=4)

    def test_gamma_return_float_on_int_input(self):
        self.assertTrue(type(gamma(10, 1)) == float)

    def test_sound(self):
        self.assertAlmostEqual(sound(2.98, 1.1, 1200), 1983.3305, places=4)

    def test_partial(self):
        self.assertAlmostEqual(partial([1], 1.2, [0], self.mock_database), [1.2], places=4)

    def test_partial_2(self):
        tmp = partial([1, 3], 1.2, [0, 1], self.mock_database)
        for val1, val2 in zip(tmp, [.3, .9]):
            self.assertAlmostEqual(val1, val2, places=4)

    @patch('builtins.print')
    def test_partial_species_no_gas(self, mock_print):
        with self.assertRaises(SystemExit):
            partial([1, 3], 1.2, [0, 2], self.mock_database)

    @patch('HGSpy.hgs.HGSData.coef')
    @patch('HGSpy.hgs.HGSData.g')
    @patch('HGSpy.hgs_prop.partial')
    @patch('HGSpy.hgs.HGSData.s')
    @patch('HGSpy.hgs.HGSData.h')
    @patch('HGSpy.hgs.HGSData.cv')
    @patch('HGSpy.hgs.HGSData.cp')
    @patch('HGSpy.hgs.HGSData')
    def test_hgs_props_ids_all(self, mock_hgsdata, mock_cp, mock_cv, mock_h, mock_s, mock_partial,
                               mock_g, mock_coef):
        mock_cp.side_effect = [10, 20, 10, 20, 10, 20]
        mock_cv.side_effect = [10, 20, 10, 20, 10, 20]
        mock_h.side_effect = [10, 20]
        mock_partial.side_effect = [[1, 2], [1, 2]]
        mock_s.side_effect = [10, 20]
        mock_g.side_effect = [10, 20]
        mock_coef.side_effect = [10, 20]
        _ = hgs_prop_ids([0, 1], [1, 1], [300, 300], [2, 2], [], mock_hgsdata)

        mock_cp.assert_called()
        self.assertEqual(mock_cp.call_count, 6)
        mock_cv.assert_called()
        self.assertEqual(mock_cv.call_count, 6)
        mock_h.assert_called()
        self.assertEqual(mock_h.call_count, 2)
        mock_partial.assert_called()
        self.assertEqual(mock_partial.call_count, 2)
        mock_s.assert_called()
        self.assertEqual(mock_s.call_count, 2)
        mock_g.assert_called()
        self.assertEqual(mock_g.call_count, 2)
        mock_coef.assert_not_called()

    @patch('HGSpy.hgs.HGSData.coef')
    @patch('HGSpy.hgs.HGSData.g')
    @patch('HGSpy.hgs_prop.partial')
    @patch('HGSpy.hgs.HGSData.s')
    @patch('HGSpy.hgs.HGSData.h')
    @patch('HGSpy.hgs.HGSData.cv')
    @patch('HGSpy.hgs.HGSData.cp')
    @patch('HGSpy.hgs.HGSData')
    def test_hgs_props_ids_mm(self, mock_hgsdata, mock_cp, mock_cv, mock_h, mock_s, mock_partial,
                              mock_g, mock_coef):
        mock_cp.side_effect = [10, 20, 10, 20, 10, 20]
        mock_cv.side_effect = [10, 20, 10, 20, 10, 20]
        mock_h.side_effect = [10, 20]
        mock_partial.side_effect = [[1, 2], [1, 2]]
        mock_s.side_effect = [10, 20]
        mock_g.side_effect = [10, 20]
        mock_coef.side_effect = [10, 20]
        _ = hgs_prop_ids([0, 1], [1, 1], [300, 300], [2, 2], ['mm'], mock_hgsdata)

        mock_cp.assert_not_called()
        mock_cv.assert_not_called()
        mock_h.assert_not_called()
        mock_partial.assert_not_called()
        mock_s.assert_not_called()
        mock_g.assert_not_called()
        mock_coef.assert_not_called()

    @patch('HGSpy.hgs.HGSData.coef')
    @patch('HGSpy.hgs.HGSData.g')
    @patch('HGSpy.hgs_prop.partial')
    @patch('HGSpy.hgs.HGSData.s')
    @patch('HGSpy.hgs.HGSData.h')
    @patch('HGSpy.hgs.HGSData.cv')
    @patch('HGSpy.hgs.HGSData.cp')
    @patch('HGSpy.hgs.HGSData')
    def test_hgs_props_ids_cp(self, mock_hgsdata, mock_cp, mock_cv, mock_h, mock_s,
                              mock_partial, mock_g, mock_coef):
        mock_cp.side_effect = [10, 20, 10, 20, 10, 20]
        mock_cv.side_effect = [10, 20, 10, 20, 10, 20]
        mock_h.side_effect = [10, 20]
        mock_partial.side_effect = [[1, 2], [1, 2]]
        mock_s.side_effect = [10, 20]
        mock_g.side_effect = [10, 20]
        mock_coef.side_effect = [10, 20]
        _ = hgs_prop_ids([0, 1], [1, 1], [300, 300], [2, 2], ['cp'], mock_hgsdata)

        mock_cp.assert_called()
        self.assertEqual(mock_cp.call_count, 2)
        mock_cv.assert_not_called()
        mock_h.assert_not_called()
        mock_partial.assert_not_called()
        mock_s.assert_not_called()
        mock_g.assert_not_called()
        mock_coef.assert_not_called()

    @patch('HGSpy.hgs.HGSData.coef')
    @patch('HGSpy.hgs.HGSData.g')
    @patch('HGSpy.hgs_prop.partial')
    @patch('HGSpy.hgs.HGSData.s')
    @patch('HGSpy.hgs.HGSData.h')
    @patch('HGSpy.hgs.HGSData.cv')
    @patch('HGSpy.hgs.HGSData.cp')
    @patch('HGSpy.hgs.HGSData')
    def test_hgs_props_ids_cv(self, mock_hgsdata, mock_cp, mock_cv, mock_h, mock_s,
                              mock_partial, mock_g, mock_coef):
        mock_cp.side_effect = [10, 20, 10, 20, 10, 20]
        mock_cv.side_effect = [10, 20, 10, 20, 10, 20]
        mock_h.side_effect = [10, 20]
        mock_partial.side_effect = [[1, 2], [1, 2]]
        mock_s.side_effect = [10, 20]
        mock_g.side_effect = [10, 20]
        mock_coef.side_effect = [10, 20]
        _ = hgs_prop_ids([0, 1], [1, 1], [300, 300], [2, 2], ['cv'], mock_hgsdata)

        mock_cp.assert_not_called()
        mock_cv.assert_called()
        self.assertEqual(mock_cv.call_count, 2)
        mock_h.assert_not_called()
        mock_partial.assert_not_called()
        mock_s.assert_not_called()
        mock_g.assert_not_called()
        mock_coef.assert_not_called()

    @patch('HGSpy.hgs.HGSData.coef')
    @patch('HGSpy.hgs.HGSData.g')
    @patch('HGSpy.hgs_prop.partial')
    @patch('HGSpy.hgs.HGSData.s')
    @patch('HGSpy.hgs.HGSData.h')
    @patch('HGSpy.hgs.HGSData.cv')
    @patch('HGSpy.hgs.HGSData.cp')
    @patch('HGSpy.hgs.HGSData')
    def test_hgs_props_ids_h(self, mock_hgsdata, mock_cp, mock_cv, mock_h, mock_s, mock_partial,
                             mock_g, mock_coef):
        mock_cp.side_effect = [10, 20, 10, 20, 10, 20]
        mock_cv.side_effect = [10, 20, 10, 20, 10, 20]
        mock_h.side_effect = [10, 20]
        mock_partial.side_effect = [[1, 2], [1, 2]]
        mock_s.side_effect = [10, 20]
        mock_g.side_effect = [10, 20]
        mock_coef.side_effect = [10, 20]
        _ = hgs_prop_ids([0, 1], [1, 1], [300, 300], [2, 2], ['h'], mock_hgsdata)

        mock_cp.assert_not_called()
        mock_cv.assert_not_called()
        mock_h.assert_called()
        self.assertEqual(mock_h.call_count, 2)
        mock_partial.assert_not_called()
        mock_s.assert_not_called()
        mock_g.assert_not_called()
        mock_coef.assert_not_called()

    @patch('HGSpy.hgs.HGSData.coef')
    @patch('HGSpy.hgs.HGSData.g')
    @patch('HGSpy.hgs_prop.partial')
    @patch('HGSpy.hgs.HGSData.s')
    @patch('HGSpy.hgs.HGSData.h')
    @patch('HGSpy.hgs.HGSData.cv')
    @patch('HGSpy.hgs.HGSData.cp')
    @patch('HGSpy.hgs.HGSData')
    def test_hgs_props_ids_s(self, mock_hgsdata, mock_cp, mock_cv, mock_h, mock_s, mock_partial,
                             mock_g, mock_coef):
        mock_cp.side_effect = [10, 20, 10, 20, 10, 20]
        mock_cv.side_effect = [10, 20, 10, 20, 10, 20]
        mock_h.side_effect = [10, 20]
        mock_partial.side_effect = [[1, 2], [1, 2]]
        mock_s.side_effect = [10, 20]
        mock_g.side_effect = [10, 20]
        mock_coef.side_effect = [10, 20]
        _ = hgs_prop_ids([0, 1], [1, 1], [300, 300], [2, 2], ['s'], mock_hgsdata)

        mock_cp.assert_not_called()
        mock_cv.assert_not_called()
        mock_h.assert_not_called()
        mock_partial.assert_called()
        self.assertEqual(mock_partial.call_count, 1)
        mock_s.assert_called()
        self.assertEqual(mock_s.call_count, 2)
        mock_g.assert_not_called()
        mock_coef.assert_not_called()

    @patch('HGSpy.hgs.HGSData.coef')
    @patch('HGSpy.hgs.HGSData.g')
    @patch('HGSpy.hgs_prop.partial')
    @patch('HGSpy.hgs.HGSData.s')
    @patch('HGSpy.hgs.HGSData.h')
    @patch('HGSpy.hgs.HGSData.cv')
    @patch('HGSpy.hgs.HGSData.cp')
    @patch('HGSpy.hgs.HGSData')
    def test_hgs_props_ids_g(self, mock_hgsdata, mock_cp, mock_cv, mock_h, mock_s, mock_partial,
                             mock_g, mock_coef):
        mock_cp.side_effect = [10, 20, 10, 20, 10, 20]
        mock_cv.side_effect = [10, 20, 10, 20, 10, 20]
        mock_h.side_effect = [10, 20]
        mock_partial.side_effect = [[1, 2], [1, 2]]
        mock_s.side_effect = [10, 20]
        mock_g.side_effect = [10, 20]
        mock_coef.side_effect = [10, 20]
        _ = hgs_prop_ids([0, 1], [1, 1], [300, 300], [2, 2], ['g'], mock_hgsdata)

        mock_cp.assert_not_called()
        mock_cv.assert_not_called()
        mock_h.assert_not_called()
        mock_partial.assert_called()
        self.assertEqual(mock_partial.call_count, 1)
        mock_s.assert_not_called()
        mock_g.assert_called()
        self.assertEqual(mock_g.call_count, 2)
        mock_coef.assert_not_called()

    @patch('HGSpy.hgs.HGSData.coef')
    @patch('HGSpy.hgs.HGSData.g')
    @patch('HGSpy.hgs_prop.partial')
    @patch('HGSpy.hgs.HGSData.s')
    @patch('HGSpy.hgs.HGSData.h')
    @patch('HGSpy.hgs.HGSData.cv')
    @patch('HGSpy.hgs.HGSData.cp')
    @patch('HGSpy.hgs.HGSData')
    def test_hgs_props_ids_rg(self, mock_hgsdata, mock_cp, mock_cv, mock_h, mock_s, mock_partial,
                              mock_g, mock_coef):
        mock_cp.side_effect = [10, 20, 10, 20, 10, 20]
        mock_cv.side_effect = [10, 20, 10, 20, 10, 20]
        mock_h.side_effect = [10, 20]
        mock_partial.side_effect = [[1, 2], [1, 2]]
        mock_s.side_effect = [10, 20]
        mock_g.side_effect = [10, 20]
        mock_coef.side_effect = [10, 20]
        _ = hgs_prop_ids([0, 1], [1, 1], [300, 300], [2, 2], ['rg'], mock_hgsdata)

        mock_cp.assert_not_called()
        mock_cv.assert_not_called()
        mock_h.assert_not_called()
        mock_partial.assert_not_called()
        mock_s.assert_not_called()
        mock_g.assert_not_called()
        mock_coef.assert_not_called()

    @patch('HGSpy.hgs.HGSData.coef')
    @patch('HGSpy.hgs.HGSData.g')
    @patch('HGSpy.hgs_prop.partial')
    @patch('HGSpy.hgs.HGSData.s')
    @patch('HGSpy.hgs.HGSData.h')
    @patch('HGSpy.hgs.HGSData.cv')
    @patch('HGSpy.hgs.HGSData.cp')
    @patch('HGSpy.hgs.HGSData')
    def test_hgs_props_ids_gamma(self, mock_hgsdata, mock_cp, mock_cv, mock_h, mock_s,
                                 mock_partial,
                                 mock_g, mock_coef):
        mock_cp.side_effect = [10, 20, 10, 20, 10, 20]
        mock_cv.side_effect = [10, 20, 10, 20, 10, 20]
        mock_h.side_effect = [10, 20]
        mock_partial.side_effect = [[1, 2], [1, 2]]
        mock_s.side_effect = [10, 20]
        mock_g.side_effect = [10, 20]
        mock_coef.side_effect = [10, 20]
        _ = hgs_prop_ids([0, 1], [1, 1], [300, 300], [2, 2], ['gamma'], mock_hgsdata)

        mock_cp.assert_called()
        self.assertEqual(mock_cp.call_count, 2)
        mock_cv.assert_called()
        self.assertEqual(mock_cv.call_count, 2)
        mock_h.assert_not_called()
        mock_partial.assert_not_called()
        mock_s.assert_not_called()
        mock_g.assert_not_called()
        mock_coef.assert_not_called()

    @patch('HGSpy.hgs.HGSData.coef')
    @patch('HGSpy.hgs.HGSData.g')
    @patch('HGSpy.hgs_prop.partial')
    @patch('HGSpy.hgs.HGSData.s')
    @patch('HGSpy.hgs.HGSData.h')
    @patch('HGSpy.hgs.HGSData.cv')
    @patch('HGSpy.hgs.HGSData.cp')
    @patch('HGSpy.hgs.HGSData')
    def test_hgs_props_ids_a(self, mock_hgsdata, mock_cp, mock_cv, mock_h, mock_s,
                             mock_partial,
                             mock_g, mock_coef):
        mock_cp.side_effect = [10, 20, 10, 20, 10, 20]
        mock_cv.side_effect = [10, 20, 10, 20, 10, 20]
        mock_h.side_effect = [10, 20]
        mock_partial.side_effect = [[1, 2], [1, 2]]
        mock_s.side_effect = [10, 20]
        mock_g.side_effect = [10, 20]
        mock_coef.side_effect = [10, 20]
        _ = hgs_prop_ids([0, 1], [1, 1], [300, 300], [2, 2], ['a'], mock_hgsdata)

        mock_cp.assert_called()
        self.assertEqual(mock_cp.call_count, 2)
        mock_cv.assert_called()
        self.assertEqual(mock_cv.call_count, 2)
        mock_h.assert_not_called()
        mock_partial.assert_not_called()
        mock_s.assert_not_called()
        mock_g.assert_not_called()
        mock_coef.assert_not_called()

    @patch('HGSpy.hgs.HGSData.coef')
    @patch('HGSpy.hgs.HGSData.g')
    @patch('HGSpy.hgs_prop.partial')
    @patch('HGSpy.hgs.HGSData.s')
    @patch('HGSpy.hgs.HGSData.h')
    @patch('HGSpy.hgs.HGSData.cv')
    @patch('HGSpy.hgs.HGSData.cp')
    @patch('HGSpy.hgs.HGSData')
    def test_hgs_props_ids_coef(self, mock_hgsdata, mock_cp, mock_cv, mock_h, mock_s,
                                mock_partial,
                                mock_g, mock_coef):
        mock_cp.side_effect = [10, 20, 10, 20, 10, 20]
        mock_cv.side_effect = [10, 20, 10, 20, 10, 20]
        mock_h.side_effect = [10, 20]
        mock_partial.side_effect = [[1, 2], [1, 2]]
        mock_s.side_effect = [10, 20]
        mock_g.side_effect = [10, 20]
        mock_coef.side_effect = [10, 20]
        _ = hgs_prop_ids([0, 1], [1, 1], [300, 300], [2, 2], ['coef'], mock_hgsdata)

        mock_cp.assert_not_called()
        mock_cv.assert_not_called()
        mock_h.assert_not_called()
        mock_partial.assert_not_called()
        mock_s.assert_not_called()
        mock_g.assert_not_called()
        mock_coef.assert_called()
        self.assertEqual(mock_coef.call_count, 2)

    @parameterized.expand([
        ('mm',), ('cp',), ('cv',), ('h',), ('s',), ('g',), ('rg',), ('gamma',), ('a',),
        ('coef',)
    ])
    @patch('HGSpy.hgs.HGSData.rebuild')
    @patch('HGSpy.hgs.HGSData.id')
    @patch('HGSpy.hgs.HGSData.__len__')
    @patch('HGSpy.hgs.HGSData')
    @patch('HGSpy.hgs_prop.hgs_prop_ids')
    def test_hgs_single_species(self, properties, mock_hgs_prop_ids, mock_hgsdata, mock_len,
                                mock_id, mock_rebuild):
        mock_id.side_effect = [[0, 1]]
        mock_len.side_effect = [2]
        mock_hgs_prop_ids.side_effect = [['success']]
        res = hgs_single(['species1', 'species2'], properties, [300], 12, mock_hgsdata)
        mock_hgs_prop_ids.assert_called()
        self.assertEqual(mock_hgs_prop_ids.call_count, 1)
        mock_id.assert_called()
        self.assertEqual(mock_id.call_count, 1)
        mock_rebuild.assert_not_called()
        self.assertEqual(mock_rebuild.call_count, 0)
        self.assertEqual(res, 'success')

    @parameterized.expand([
        ('mm',), ('cp',), ('cv',), ('h',), ('s',), ('g',), ('rg',), ('gamma',), ('a',),
        ('coef',)
    ])
    @patch('HGSpy.hgs.HGSData.rebuild')
    @patch('HGSpy.hgs.HGSData.id')
    @patch('HGSpy.hgs.HGSData.__len__')
    @patch('HGSpy.hgs.HGSData')
    @patch('HGSpy.hgs_prop.hgs_prop_ids')
    def test_hgs_single_mixture(self, properties, mock_hgs_prop_ids, mock_hgsdata, mock_len,
                                mock_id, mock_rebuild):
        mock_rebuild.side_effect = [(['species1', 'species2'], [0.5, 0.5], 300)]
        mock_id.side_effect = [[0, 3], [0, 1]]
        mock_len.side_effect = [2]
        mock_hgs_prop_ids.side_effect = [['success']]
        res = hgs_single(['species1', 'species2'], properties, [300], 12, mock_hgsdata)
        mock_hgs_prop_ids.assert_called()
        self.assertEqual(mock_hgs_prop_ids.call_count, 1)
        mock_id.assert_called()
        self.assertEqual(mock_id.call_count, 2)
        mock_rebuild.assert_called()
        self.assertEqual(mock_rebuild.call_count, 1)
        self.assertEqual(res, 'success')

    @patch('HGSpy.hgs.HGSData.rebuild')
    @patch('HGSpy.hgs.HGSData.id')
    @patch('HGSpy.hgs.HGSData.__len__')
    @patch('HGSpy.hgs.HGSData')
    @patch('HGSpy.hgs_prop.hgs_prop_ids')
    @patch('builtins.print')
    def test_hgs_single_wrong_property(self, mock_print, mock_hgs_prop_ids, mock_hgsdata, mock_len,
                                       mock_id,
                                       mock_rebuild):
        mock_rebuild.side_effect = [(['species1', 'species2'], [0.5, 0.5], 300)]
        mock_id.side_effect = [[0, 3], [0, 1]]
        mock_len.side_effect = [2]
        mock_hgs_prop_ids.side_effect = [['success']]
        with self.assertRaises(SystemExit):
            _ = hgs_single(['species1', 'species2'], 'not_exist_prop', [300], 12, mock_hgsdata)
        mock_hgs_prop_ids.assert_not_called()

    @parameterized.expand([
        ('species1', 1, 300, 12),
        ('species1', 1, [300], [12]),
        (['species1'], [1], 300, 12),
        (['species1'], [1], [300], [12]),
        (['species1', 'species2'], [1, 2], 300, 12),
        (['species1', 'species2'], [1, 2], [300], [12]),
        (['species1', 'species2'], [1, 2], [300, 300], [12, 12])
    ])
    @patch('HGSpy.hgs.HGSData.rebuild')
    @patch('HGSpy.hgs.HGSData.id')
    @patch('HGSpy.hgs.HGSData.__len__')
    @patch('HGSpy.hgs.HGSData')
    @patch('HGSpy.hgs_prop.hgs_prop_ids')
    @patch('builtins.print')
    def test_hgs_prop_species(self, config_species, config_mols, config_temp, config_pressure,
                              mock_print, mock_hgs_prop_ids, mock_hgsdata, mock_len,
                              mock_id, mock_rebuild):
        mock_id.side_effect = [[0]]
        mock_len.side_effect = [2]
        mock_hgs_prop_ids.side_effect = ['success']
        res = hgs_prop(config_species, config_mols, config_temp, config_pressure, 'mm',
                       hgs_data=mock_hgsdata)
        mock_hgs_prop_ids.assert_called()
        self.assertEqual(mock_hgs_prop_ids.call_count, 1)
        mock_id.assert_called()
        self.assertEqual(mock_id.call_count, 1)
        mock_rebuild.assert_not_called()
        self.assertEqual(mock_rebuild.call_count, 0)
        self.assertEqual(res, 'success')

    @parameterized.expand([
        ('species1', 1, 300, 12),
        ('species1', 1, [300], [12]),
        (['species1'], [1], 300, 12),
        (['species1'], [1], [300], [12]),
        (['species1', 'species2'], [1, 2], 300, 12),
        (['species1', 'species2'], [1, 2], [300], [12]),
        (['species1', 'species2'], [1, 2], [300, 300], [12, 12])
    ])
    @patch('HGSpy.hgs.HGSData.rebuild')
    @patch('HGSpy.hgs.HGSData.id')
    @patch('HGSpy.hgs.HGSData.__len__')
    @patch('HGSpy.hgs.HGSData')
    @patch('HGSpy.hgs_prop.hgs_prop_ids')
    @patch('builtins.print')
    def test_hgs_prop_mixture(self, config_species, config_mols, config_temp, config_pressure,
                              mock_print, mock_hgs_prop_ids, mock_hgsdata, mock_len,
                              mock_id, mock_rebuild):
        mock_rebuild.side_effect = [(['species1', 'species2'], [0.5, 0.5], 300)]
        mock_id.side_effect = [[0, 3], [0, 1]]
        mock_len.side_effect = [2]
        mock_hgs_prop_ids.side_effect = ['success']
        res = hgs_prop(config_species, config_mols, config_temp, config_pressure, 'mm',
                       hgs_data=mock_hgsdata)
        mock_hgs_prop_ids.assert_called()
        self.assertEqual(mock_hgs_prop_ids.call_count, 1)
        mock_id.assert_called()
        self.assertEqual(mock_id.call_count, 2)
        mock_rebuild.assert_called()
        self.assertEqual(mock_rebuild.call_count, 1)
        self.assertEqual(res, 'success')

    @patch('HGSpy.hgs.HGSData.rebuild')
    @patch('HGSpy.hgs.HGSData.id')
    @patch('HGSpy.hgs.HGSData.__len__')
    @patch('HGSpy.hgs.HGSData')
    @patch('HGSpy.hgs_prop.hgs_prop_ids')
    @patch('builtins.print')
    def test_hgs_prop_length_diff(self, mock_print, mock_hgs_prop_ids, mock_hgsdata, mock_len,
                                  mock_id, mock_rebuild):
        mock_rebuild.side_effect = [(['species1', 'species2'], [0.5, 0.5], 300)]
        mock_id.side_effect = [[0, 3], [0, 1]]
        mock_len.side_effect = [2]
        mock_hgs_prop_ids.side_effect = [['success']]
        with self.assertRaises(SystemExit):
            _ = hgs_prop(['species1', 'species2'], [1], 300, 12, 'mm',
                         hgs_data=mock_hgsdata)
        mock_hgs_prop_ids.assert_not_called()


if __name__ == '__main__':
    unittest.main()
