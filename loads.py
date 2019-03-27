"""
author: @Udayan Sarin
defines load types for beam_io. Lists load properties and features
"""
import configparser as config
import os


class SetupControl:
    """
    Define the types of connections and their reaction loads
    """
    def __init__(self):
        self._cfg = config.ConfigParser()
        self._cfg.optionxform = str
        self._cfg.read(os.path.join(os.getcwd(), 'setup_configurations.ini'))
        self._load_types = []
        self._support_types = []
        self._init_properties()

    def _init_properties(self):
        _pop_list = [self._load_types, self._support_types]
        _targ_list = ["LoadTypes", "SupportTypes"]
        for index, section_header in enumerate(_targ_list):
            for header, _ in self._cfg[section_header].items():
                _pop_list[index].append(header)

    def get_supports(self):
        return self._support_types

    def get_support_reactions(self, support):
        return self._cfg['SupportTypes'][support].split(' ')

    def get_loads(self):
        return self._load_types

    def get_definition_conditions(self, load):
        return int(self._cfg['LoadTypes'][load])
