"""
author: @Udayan Sarin
Contains information pertaining to units and conversions that beam_io is compatible with
"""
import configparser as config
import os


class UnitControl:
    def __init__(self):
        config_file = os.path.join(os.getcwd(), 'unit_configurations.ini')
        self._cfg = config.ConfigParser()
        self._cfg.optionxform = str
        self._cfg.read(config_file)
        self._len_units = []
        self._load_units = []
        self._stress_units = []
        self._init_units()

    def _init_units(self):
        _pop_list = [self._len_units, self._load_units, self._stress_units]
        _targ_list = ['DistanceUnits', 'ForceUnits', 'StressUnits']
        for index, section_header in enumerate(_targ_list):
            for header, _ in self._cfg[section_header].items():
                _pop_list[index].append(header)

    def get_len_units(self):
        return self._len_units

    def get_len_conversion(self, unit):
        return float(self._cfg['DistanceUnits'][unit])

    def get_force_units(self):
        return self._load_units

    def get_force_conversion(self, unit):
        return float(self._cfg['ForceUnits'][unit])

    def get_stress_units(self):
        return self._stress_units

    def get_stress_conversion(self, unit):
        return float(self._cfg['StressUnits'][unit])

    @staticmethod
    def get_combined_unit(quantity, length, force):
        if (quantity == 'Moment Load') or (quantity == 'Torsion Load'):
            return f'{force}{length}'
        if (quantity == 'Vertical Load') or (quantity == 'Axial Load'):
            return f'{force}'
        if quantity == 'Distributed Load':
            return f'{force}/{length}'
        return False

    def get_combined_conversion(self, _load_type, length, force):
        if (_load_type == 'Vertical Load') or (_load_type == 'Axial Load'):
            return float(self.get_force_conversion(force))
        if (_load_type == 'Moment Load') or (_load_type == 'Torsion Load'):
            return float(self.get_force_conversion(force)*self.get_len_conversion(length))
        if _load_type == 'Distributed Load':
            return float(self.get_force_conversion(force)/self.get_len_conversion(length))