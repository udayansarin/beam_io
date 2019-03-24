"""
author: @Udayan Sarin
Includes functions and classes to evalute shear and axial stresses developed on a cross-section for given loading
conditions. Makes use of stress development due to bending moment, shear force and axial force.
"""
import numpy as np


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
        # define the size of the cross-section
        for rect in self._rects:
            coor_y = rect['y_loc']+rect['h']
            coor_x = rect['x_loc']+rect['b']
            if coor_y > self._len_cx:
                self._len_cx = coor_y
            if coor_x > self._width_cx:
                self._width_cx = coor_x

    def get_stresses(self):
        """
        calculate shear and axial stresses developed at a given location on the beam, for a given cross-section
        :return:
        """
        cx_n = np.arange(0, self._len_cx, ((self._width_cx+self._len_cx)/2000))
        cx_h = np.arange(0, self._width_cx, ((self._width_cx+self._len_cx)/2000))
        shear_stress = [[0] * len(cx_h)] * len(cx_n)
        axial_stress = [[0] * len(cx_h)] * len(cx_n)
        for index, position in enumerate(cx_n):
            thick = 0
            a_dash = 0
            ya_dash = 0
            for rect in self._rects:
                height = rect['h']
                width = rect['b']
                coor_y = rect['y_loc']
                if coor_y <= position < (coor_y + height):
                    thick = width
                    a_dash += ((position - coor_y) * width)
                    ya_dash += ((position + coor_y) / 2) * (position - coor_y) * width
                    break
                if position >= (coor_y + height):
                    thick = width
                    a_dash += (width * height)
                    ya_dash += (((height / 2) + coor_y) * width * height)
            if a_dash != 0:
                y_dsh = ya_dash/a_dash
                y_dash = self._prop['yNA'] - y_dsh
            else:
                y_dash = 0
            q = y_dash*a_dash
            shear_stress[len(cx_n)-index-1] = [self._sf*q/(self._prop['I']*thick)]*len(cx_h)
            axial_stress[len(cx_n)-index-1] = [(self._bm*(self._prop['yNA']-position)/self._prop['I'])+(
                self._af/self._prop['A'])]*len(cx_h)
            return shear_stress, axial_stress, [self._width_cx, self._len_cx]
