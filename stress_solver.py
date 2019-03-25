"""
author: @Udayan Sarin
Includes functions and classes to evalute shear and axial stresses developed on a cross-section for given loading
conditions. Makes use of stress development due to bending moment, shear force and axial force.
"""
import numpy as np
import units


class StressSolver:
    """
    instances of this class possess the required functions to evaluate axial and shear stress for a given cross-section
    """
    def __init__(self, _rectangles, _properties, _bm, _sf, _af):
        self._rects = _rectangles
        self._prop = _properties
        self._bm = _bm
        self._sf = _sf
        self._af = _af
        self._len_cx = 0
        self._width_cx = 0
        self._units = units.UnitControl()
        # define the size of the cross-section
        for rect in self._rects:
            coor_y = rect['y_loc']+rect['h']
            coor_x = rect['x_loc']+rect['b']
            if coor_y > self._len_cx:
                self._len_cx = coor_y
            if coor_x > self._width_cx:
                self._width_cx = coor_x

    def get_stresses(self, len_unit, stress_unit):
        """
        calculate shear and axial stresses developed at a given location on the beam, for a given cross-section
        :return:
        """
        cx_w = np.arange(0, self._len_cx, ((self._width_cx + self._len_cx) / 1000))
        cx_h = np.arange(0, self._width_cx, ((self._width_cx + self._len_cx) / 100))
        shear_stress = [[0] * len(cx_h)] * len(cx_w)
        axial_stress = [[0] * len(cx_h)] * len(cx_w)

        for index in range(len(cx_w)):
            thick = 0
            position = cx_w[index]*self._units.get_len_conversion(len_unit)
            a_dash = 0
            ya_dash = 0
            for rect in self._rects:
                height = rect['h']
                width = rect['b']
                y_coor = rect['y_loc']
                if y_coor <= position < (y_coor + height):
                    thick = width
                    a_dash += ((position - y_coor) * width)
                    ya_dash += ((position + y_coor) / 2) * (position - y_coor) * width
                    break
                if position >= (y_coor + height):
                    thick = width
                    a_dash += (width * height)
                    ya_dash += (((height / 2) + y_coor) * width * height)
            if a_dash != 0:
                y_dash = ya_dash / a_dash
                y_dash = self._prop['yNA'] - y_dash
            else:
                y_dash = 0
            q = y_dash * a_dash
            shear_stress[len(cx_w) - index - 1] = [(self._sf * q / (self._prop['I'] * thick)) /
                                                   self._units.get_stress_conversion(stress_unit)] * len(cx_h)
            axial_stress[len(cx_w) - index - 1] = [((self._bm * (self._prop['yNA'] - position) / self._prop['I']) +
                                                    self._af / self._prop['A']) /
                                                   self._units.get_stress_conversion(stress_unit)]*len(cx_h)

        return shear_stress, axial_stress, [self._width_cx, self._len_cx]
