"""
author: @Udayan Sarin
Macaulay function/singularity function solver for beams. Applies knowledge from Shingley's Mechanical Design 10th Ed.
Principle - Solving structural loading problems where the load on the structure develops more unknowns than force and
moment balance equations. This requires that assumptions be made about deflection and slope of the structure due to
loading, to allow for the problem to be solved.

The script applies concepts of linear algebra and integration to develop shear force, bending moment and deflection
models for the structure being analyzed
"""
import numpy as np
import loads
import matplotlib.pyplot as plt


def _h(limit, current_position):
    """
    heaviside function
    :param limit: signal on limit
    :param current_position: test value
    :return: 1 if test value is greater than the limit
    """
    return 0 if np.sign(limit-current_position) > 0 else 1


class VerticalLoad:
    """
    defines properties for a vertical load
    """
    def __init__(self, loc, magnitude):
        """
        define properties of a vertical load
        :param loc: location in distance
        :param magnitude: in force
        """
        self._location = loc
        self._value = magnitude

    def singularity(self, x):
        """
        macaulay deflection result
        :param x: test location
        :return: deflection due to vertical load
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value/6*(x-self._location)**3

    def singularity_slope(self, x):
        """
        macaulay slope result
        :param x: test location
        :return: slope due to vertical load
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value/2*(x-self._location)**2

    def get_point_loading(self, mom_bal=False):
        """
        equivalent point loading due to a vertical load
        :param mom_bal: boolean, True if moment balance is being parsed
        :return: equivalent point loading
        """
        return self._value*(self._location if mom_bal else 1)

    def v_shear(self, x):
        """
        macaulay shear due to a vertical load
        :param x: test location
        :return: macaulay shear
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value

    def v_moment(self, x):
        """
        macaulay moment due to a vertical load
        :param x: test location
        :return: macaulay moment
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value*(x-self._location)


class AxialLoad:
    """
    define properties for an axial load
    """
    def __init__(self, loc):
        """
        additional functionality for axial loading in development
        :param loc:
        """
        self._location = loc
        print('hi')

    @staticmethod
    def get_point_loading(mom_bal=False):
        return 0


class Moment:
    """
    defines properties for a moment load
    """
    def __init__(self, loc, magnitude):
        """
        define properties of a moment load
        :param loc: location in distance
        :param magnitude: magnitude in force*distance
        """
        self._location = loc
        self._value = magnitude

    def singularity(self, x):
        """
        macaulay deflection result
        :param x: test location
        :return: deflection due to moment load
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value/2*(x-self._location)**2

    def singularity_slope(self, x):
        """
        macaulay slope result
        :param x: test location
        :return: slope due to moment load
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value*(x-self._location)

    def get_point_loading(self, mom_bal=False):
        """
        equivalent point loading due to a moment load
        :param mom_bal: boolean, True if moment balance is being parsed
        :return: equivalent point load
        """
        return self._value*(1 if mom_bal else 0)

    def m_moment(self, x):
        """
        macaulay moment due to a moment load
        :param x: test location
        :return: macaulay moment
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*-self._value

    @staticmethod
    def m_shear(self, x):
        if self._location == x:
            return 0
        return 0


class RectangleLoad:
    """
    defines properties for a rectangular distributed load
    """
    def __init__(self, loc, magnitude, end_loc):
        """
        define properties of a rectangular load
        :param loc: location in distance
        :param magnitude: magnitude in force/distance
        :param end_loc: end of the beam on which the force is generated
        """
        self._location = loc
        self._end = end_loc
        self._value = magnitude

    def singularity(self, x):
        """
        macaulay deflection result
        :param x: test location
        :return: deflection due to rectangle load
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*(self._value/24)*(x-self._location)**4

    def singularity_slope(self, x):
        """
        macaulay slope result
        :param x: test location
        :return: slope due to rectangle load
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value/6*(x-self._location)**3

    def point_eq(self):
        """
        equivalent force due to a rectangular load
        :return:
        """
        return self._value*(self._end - self._location)

    def get_point_loading(self, mom_bal=False):
        """
        equivalent point loading due to a rectangular load
        :param mom_bal: boolean, True if moment balance is being parsed
        :return: equivalent point load
        """
        return self.point_eq()*((self._location+((self._end - self._location)/2)) if mom_bal else 1)

    def v_shear(self, x):
        """
        macaulay shear due to a rectangular load
        :param x: test location
        :return: macaulay shear
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value*(x-self._location)

    def v_moment(self, x):
        """
        macaulay moment due to a rectangle load
        :param x: test location
        :return: macaulay moment
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value/2*(x-self._location)**2


class TriangleLoad:
    """
    defines properties for a triangular distributed load
    """
    def __init__(self, loc, slope, end_loc):
        """
        define properties of a triangular load
        :param loc: location in distance
        :param slope: slope in (force/distance)/distance
        :param end_loc: end of the beam on which the force is generated
        """
        self._location = loc
        self._end = end_loc
        self._value = slope

    def singularity(self, x):
        """
        macaulay deflection result
        :param x: test location
        :return: deflection due to triangle load
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*(self._value/120)*(x-self._location)**5

    def singularity_slope(self, x):
        """
        macaulay slope result
        :param x: test location
        :return: slope due to rectangle load
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value/24*(x-self._location)**4

    def point_eq(self):
        """
        equivalent force due to a triangular force
        :return:
        """
        return self._value*0.5*(self._end - self._location)**2

    def get_point_loading(self, mom_bal=False):
        """
        equivalent point loading due to a triangular load
        :param mom_bal: boolean, True if moment balance is being parsed
        :return: equivalent point load
        """
        return self.point_eq()*((self._location + (self._end - self._location)*2/3) if mom_bal else 1)

    def v_shear(self, x):
        """
        macaulay shear due to a triangular load
        :param x: test location
        :return: macaulay shear
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value/2*(x-self._location)**2

    def v_moment(self, x):
        """
        macaulay moment due to a triangular load
        :param x: test location
        :return: macaulay moment
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._value/6*(x-self._location)**3


class DistributedLoad:
    """
    Defines properties for an arbitrary distributed load
    """
    def __init__(self, loc1, loc2, val1, val2, b_len):
        """
        splits a distributed load into its sub-components
        :param loc1: start location in distance
        :param loc2: end location in distance
        :param val1: start value in force/distance
        :param val2: end value in force/distance
        :param b_len: length of the beam being parsed
        """
        self._location1 = loc1
        self._location2 = loc2
        self._beam_len = b_len
        self._value1 = val1
        self._value2 = val2
        self._slope = (val2-val1)/(loc2-loc1)
        self._model = [
            RectangleLoad(self._location1, self._value1, self._beam_len),
            RectangleLoad(self._location2, -self._value2, self._beam_len),
            TriangleLoad(self._location1, self._slope, self._beam_len),
            TriangleLoad(self._location2, -self._slope, self._beam_len)
        ]

    def singularity(self, x):
        """
        macaulay deflection result
        :param x: test location
        :return: deflection due to distributed load
        """
        sing_result = 0
        for sub_load in self._model:
            sing_result += sub_load.singularity(x)
        return sing_result

    def singularity_slope(self, x):
        """
        macaulay slope result
        :param x: test location
        :return: slope due to distributed load
        """
        slope_result = 0
        for sub_load in self._model:
            slope_result += sub_load.singularity_slope(x)
        return slope_result

    def get_point_loading(self, mom_bal=False):
        """
        equivalent point loading due to a distributed load
        :param mom_bal: boolean, True if moment balance is being parsed
        :return: equivalent point load
        """
        val_to_return = 0
        for sub_load in self._model:
            val_to_return += sub_load.get_point_loading(mom_bal)
        return val_to_return

    def v_shear(self, x):
        """
        macaulay shear due to a distributed load
        :param x: test location
        :return: macaulay shear
        """
        val_to_return = 0
        for sub_load in self._model:
            val_to_return += sub_load.v_shear(x)
        return val_to_return

    def v_moment(self, x):
        """
        macaulay moment due to a distributed load
        :param x: test location
        :return: macaulay moment
        """
        val_to_return = 0
        for sub_load in self._model:
            val_to_return += sub_load.v_moment(x)
        return val_to_return


class FixedSupport:
    """
    model a cantiliver/fixed support
    """
    def __init__(self, loc):
        """
        define properties of a fixed support
        :param loc: location of the support
        """
        self._location = loc
        # initial values of unity for reactions to allow for solution of matrices
        self._reaction_moment = 1
        self._reaction_vert = 1
        self._reaction_ax = 1
        self._reaction_tor = 1

    def m_get_point_loading(self, mom_bal=False):
        """
        get point loading due to moment reaction
        :param mom_bal: bool, indicates whether moment balance is being conducted
        :return: point loading as moment for moment balance or force balance, governed by mom_bal
        """
        return self._reaction_moment*(1 if mom_bal else 0)

    def v_get_point_loading(self, mom_bal=False):
        """
        get point loading due to shear reaction
        :param mom_bal: bool, indicates whether moment balance is being conducted
        :return: point loading as moment for moment balance or force balance, governed by mom_bal
        """
        return self._reaction_vert*(self._location if mom_bal else 1)

    def m_singularity(self, x):
        """
        macaulay deflection developed due to moment reaction
        :param x: test location
        :return: macaulay deflection
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_moment/2*(x-self._location)**2

    def m_singularity_slope(self, x):
        """
        macaulay slope developed due to moment reaction
        :param x: test location
        :return: macaulay slope
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_moment*(x-self._location)

    def v_singularity(self, x):
        """
        macaulay deflection due to shear reaction
        :param x: test location
        :return: macaulay deflection
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_vert/6*(x-self._location)**3

    def v_singularity_slope(self, x):
        """
        macaulay slope due to shear reaction
        :param x: test location
        :return: macaulay slope
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_vert/2*(x-self._location)**2

    def m_shear(self, x):
        """
        macaulay shear due to moment reaction
        :param x: test location
        :return: macaulay shear
        """
        if self._location == x:
            return 0
        return 0 #_h(self._location, x)*self._reaction_moment*(x-self._location)**-1 pg 90 in the TB

    def v_shear(self, x):
        """
        maculay shear due to shear reaction
        :param x: test location
        :return: macaulay shear
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_vert

    def m_moment(self, x):
        """'
        macaulay moment due to moment reaction
        :param x: test location
        :return: macaulay moment
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_moment

    def v_moment(self, x):
        """
        macaulay moment due to shear reaction
        :param x: test location
        :return: macaulay moment
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*(x-self._location)*self._reaction_vert

    def update_values(self, moment, shear):
        """
        update value for the reactions at the support
        :param moment: update value for the moment reaction
        :param shear: update value for the shear reaction
        :return:
        """
        self._reaction_moment = -moment
        self._reaction_vert = shear


class PinnedSupport:
    """
    model a pinned structural support
    """
    def __init__(self, loc):
        """
        define properties for a pinned structural support
        :param loc: location of the support
        """
        self._location = loc
        # initial value of unity for the reaction loads
        self._reaction_vert = 1
        self._reaction_ax = 1
        self._reaction_tor = 1

    def v_singularity(self, x):
        """
        macaulay deflection due to shear reaction
        :param x: test location
        :return: macaulay deflection
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_vert/6*(x-self._location)**3

    def v_singularity_slope(self, x):
        """
        macaulay slope due to shear reaction
        :param x: test location
        :return: macaulay slope
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_vert/2*(x-self._location)**2

    def v_shear(self, x):
        """
        macaulay shear due to shear reaction
        :param x: test location
        :return: macaulay shear
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_vert

    def v_moment(self, x):
        """
        macaulay moment due to shear reaction
        :param x: test location
        :return: macaulay moment
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*(x-self._location)*self._reaction_vert

    def get_point_loading(self, mom_bal=False):
        """
        get point loading due to shear reaction developed at a pinned support
        :param mom_bal: boolean indicating whether the point loading is for a moment or force balance
        :return: point loading due to shear reaction at a pinned support
        """
        return self._reaction_vert*(self._location if mom_bal else 1)

    def update_values(self, shear):
        """
        update values of the reaction loads developed at a pinned support
        :param shear: shear value to be updated
        :return:
        """
        self._reaction_vert = shear


class RollerSupport:
    """
    model a pinned structural support
    """
    def __init__(self, loc):
        """
        define properties of a pinned structural support
        :param loc: location of the support
        """
        self._location = loc
        self._reaction_vert = 1

    def v_singularity(self, x):
        """
        macaulay deflection due to shear reaction
        :param x: test location
        :return: macaulay deflection
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_vert/6*(x-self._location)**3

    def v_singularity_slope(self, x):
        """
        macaulay slope due to shear reaction
        :param x: test location
        :return: macaulay slope
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_vert/2*(x-self._location)**2

    def v_shear(self, x):
        """
        macaulay shear due to shear reaction
        :param x: test location
        :return: macaulay sh
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*self._reaction_vert

    def v_moment(self, x):
        """
        macaulay moment due to shear reaction
        :param x: test location
        :return: macaulay moment
        """
        if self._location == x:
            return 0
        return _h(self._location, x)*(x-self._location)*self._reaction_vert

    def get_point_loading(self, mom_bal=False):
        """
        get point loading due to shear reaction developed at a roller support
        :param mom_bal: boolean indicating whether the point loading is for a moment or force balance
        :return: point loading due to shear reaction at a roller support
        """
        return self._reaction_vert*(self._location if mom_bal else 1)

    def update_values(self, shear):
        """
        update values of the reaction loads developed at a roller support
        :param shear: shear value to be updated
        :return:
        """
        self._reaction_vert = shear


class SolveProblem:
    """
    define a generic beam solver
    """
    def __init__(self, beam_len, support_dict, load_dict, grid_def, special_properties=None):
        """
        define the features which need to be solved
        :param beam_len: length of the beam being solved
        :param support_dict: dict defining the location of supports
        supports = [{'support': 'Fixed', 'location': 0.0}, {'support': 'Fixed', 'location': 12.0}]
        :param load_dict: dict defining the location of loads
        loading = [{'type': 'Vertical Load', 'val': 6000.0, 'loc': 6.0}]
        :param grid_def: definition of the grid which is being solved
        :param special_properties: dict containing information pertaining to Young's modulus and cross section geometry
        """
        self._beam_len = beam_len
        self._LHS = []
        self._RHS = []
        self._supports = support_dict
        self._loads = load_dict
        self.simple_unknowns = None  # forces and moments
        self.unknowns = None  # including constants of integration
        self.C1 = None  # Integration const 1
        self.C2 = None  # Integration const 2
        self.equations = 0
        self._setup = loads.SetupControl()
        self._beam_body = list(np.linspace(0, self._beam_len, grid_def))
        self._shear_force = list(np.zeros(len(self._beam_body)))
        self._bending_moment = list(np.zeros(len(self._beam_body)))
        self._axial_force = list(np.zeros(len(self._beam_body)))
        self._deflection = list(np.zeros(len(self._beam_body)))
        self._identify_vertical_unknowns()
        self._create_equations()

    def _identify_vertical_unknowns(self):
        """
        identify the number of unknowns which exist in the structural loading problem
        :return:
        """
        # unknowns which can be solved using structural analysis
        self.simple_unknowns = 0
        for support in self._supports:
            reaction = self._setup.get_support_reactions(support['support'])
            self.simple_unknowns += int(reaction[0]) + int(reaction[1])
        # indeterminate unkowns, to be evaluated using singularity/macaulay properties
        self.unknowns = self.simple_unknowns + 2
        for _ in range(0, self.unknowns):
            self._LHS.append(list(np.zeros(self.unknowns)))
            self._RHS.append(0)

    @staticmethod
    def cls_router(name):
        """
        link keys elements in the support dict to corresponding property classes
        :param name:
        :return:
        """
        if name == 'Pinned':
            return PinnedSupport
        if name == 'Roller':
            return RollerSupport
        if name == 'Fixed':
            return FixedSupport
        if name == 'Vertical Load':
            return VerticalLoad
        if name == 'Axial Load':
            return AxialLoad
        if name == 'Moment Load':
            return Moment
        if name == 'Distributed Load':
            return DistributedLoad

    def _create_equations(self):
        """
        develop the equations to solve the problem
        :return:
        """
        for support in self._supports:
            support['obj'] = self.cls_router(support['support'])(loc=support['location'])

        for load in self._loads:
            if load['type'] == 'Distributed Load':
                load['obj'] = self.cls_router(load['type'])(load['loc'], load['loc2'], load['val'], load['val2'],
                                                            self._beam_len)
            elif load['type'] == 'Axial Load':
                continue
            else:
                load['obj'] = self.cls_router(load['type'])(load['loc'], load['val'])
        self._basic_equations(moment_balance=False)  # force balance
        self._basic_equations(moment_balance=True)  # moment balance
        self._get_slope_equations()  # slope equation from evaluating locations of zero slope
        self._get_deflection_equations()  # slope equations from evaluating locations of zero deflection

        self._parse_results()  # solve system of linear equations and develop result datasets
        self._make_plots()  # plot results

    def _get_slope_equations(self):
        """
        identify locations of zero deflection and evaluate macaulay deflection for the beam at target locations
        :return:
        """
        # check if the equation is needed - if a fixed support exists
        for parent_support in self._supports:
            if not parent_support['support'] == "Fixed":
                continue
            indet_loc = parent_support['location']
            var_ad = 0
            for support in self._supports:
                if support['support'] == 'Fixed':
                    # account for moment and shear
                    self._LHS[self.equations][var_ad] = support['obj'].m_singularity_slope(indet_loc)
                    self._LHS[self.equations][var_ad] = support['obj'].v_singularity_slope(indet_loc)
                    var_ad += 2
                else:
                    self._LHS[self.equations][var_ad] = support['obj'].v_singularity_slope(indet_loc)
                    var_ad += 1
            self._LHS[self.equations][self.simple_unknowns] = 1
            rhs_slope_add = 0
            for load in self._loads:
                if load['type'] == 'Axial Load':
                    continue
                else:
                    rhs_slope_add += load['obj'].singularity_slope(indet_loc)
            self._RHS[self.equations] = -rhs_slope_add
            self.equations += 1

    def _get_deflection_equations(self):
        """
        identify locations where macaulay deflection can be evaluated - i.e supports as beam deflection is zero
        :return:
        """
        for parent_support in self._supports:
            indet_loc = parent_support['location']
            var_ad = 0
            for support in self._supports:
                if support['support'] == 'Fixed':
                    # account for moment and shear
                    self._LHS[self.equations][var_ad] = -support['obj'].m_singularity(indet_loc)
                    self._LHS[self.equations][var_ad + 1] = support['obj'].v_singularity(indet_loc)
                    var_ad += 2
                else:
                    self._LHS[self.equations][var_ad] = support['obj'].v_singularity(indet_loc)
                    var_ad += 1
            self._LHS[self.equations][self.simple_unknowns] = indet_loc
            self._LHS[self.equations][self.simple_unknowns + 1] = 1
            rhs_add_deflection = 0
            for load in self._loads:
                if load['type'] == 'Axial Load':
                    continue
                else:
                    rhs_add_deflection += load['obj'].singularity(indet_loc)
            self._RHS[self.equations] = -rhs_add_deflection
            self.equations += 1

    def _basic_equations(self, moment_balance=False):
        """
        perform force and moment balance on the structure, develop corresponding equations
        :param moment_balance: boolean indicating whether the analysis is a moment balance
        :return:
        """
        var_count = 0
        for support in self._supports:
            if support['support'] == 'Fixed':
                self._LHS[self.equations][var_count] = support['obj'].m_get_point_loading(mom_bal=moment_balance)
                self._LHS[self.equations][var_count + 1] = support['obj'].v_get_point_loading(mom_bal=moment_balance)
                var_count += 2
            else:
                self._LHS[self.equations][var_count] = support['obj'].get_point_loading(mom_bal=moment_balance)
                var_count += 1
        load_bal = 0
        for load in self._loads:
            if load['type'] == 'Axial Load':
                continue
            load_bal += load['obj'].get_point_loading(mom_bal=moment_balance)
        self._RHS[self.equations] = -load_bal
        self.equations += 1

    def _parse_results(self):
        """
        solve system of linea equations to evaluate reactions, update support values and determine developed constants
        of integration when macaulay moments are converted to deflection
        :return:
        """
        print("this is the matrix system we are solving")
        print(self._LHS, self._RHS)
        solution_roots = np.linalg.solve(self._LHS, self._RHS)
        print("this is the solution matrix")
        parser_count = 0
        a = list(solution_roots)
        self.C2 = a.pop()
        self.C1 = a.pop()
        print(solution_roots)
        print(self.C1)
        print(self.C2)
        for support in self._supports:
            if support['support'] == 'Fixed':
                support['obj'].update_values(moment=solution_roots[parser_count], shear=solution_roots[parser_count+1])
                parser_count += 2
            else:
                support['obj'].update_values(shear=solution_roots[parser_count])
                parser_count += 1

        for index, x_loc in enumerate(self._beam_body):
            shear_force = 0
            bending_moment = 0
            deflection = 0
            for support in self._supports:
                if support['support'] == 'Fixed':
                    shear_force += support['obj'].m_shear(x_loc) + support['obj'].v_shear(x_loc)
                    bending_moment += support['obj'].m_moment(x_loc) + support['obj'].v_moment(x_loc)
                    deflection += support['obj'].m_singularity(x_loc) + support['obj'].v_singularity(x_loc)
                else:
                    shear_force += support['obj'].v_shear(x_loc)
                    bending_moment += support['obj'].v_moment(x_loc)
                    deflection += support['obj'].v_singularity(x_loc)
            for load in self._loads:
                if load['type'] == 'Moment Load':
                    shear_force += load['obj'].m_shear(x_loc)
                    bending_moment += load['obj'].m_moment(x_loc)
                    deflection += load['obj'].singularity(x_loc)
                elif load['type'] == 'Axial Load':
                    continue
                else:
                    shear_force += load['obj'].v_shear(x_loc)
                    bending_moment += load['obj'].v_moment(x_loc)
                    deflection += load['obj'].singularity(x_loc)
            deflection += self.C1*x_loc + self.C2
            self._shear_force[index] = shear_force
            self._bending_moment[index] = bending_moment
            self._deflection[index] = deflection

    def _make_plots(self):
        """
        develop plots to display matplotlib data
        :return:
        """
        plt.plot(self._beam_body, self._shear_force)
        plt.show()
        plt.plot(self._beam_body, self._bending_moment)
        plt.show()
        plt.plot(self._beam_body, self._deflection)
        plt.show()


def dummy_run():

    # advanced example from cloud beam
    length = 5
    supports = [{'support': 'Fixed', 'location': 0.0}]
    loadsi = [{'type': 'Vertical Load', 'val': -20000.0, 'loc': 4.5},
              {'type': 'Moment Load', 'val': 40000.0, 'loc': 3.0},
              {'val2': 2000.0, 'loc2': 4.0, 'type': 'Distributed Load', 'val': 5000.0, 'loc': 1.0}]

    #a = SolveProblem(length, supports, loadsi, 1000)

    # simple indeterminate from text book
    length = 10.0
    supports = [{'support': 'Fixed', 'location': 0.0}, {'support': 'Pinned', 'location': 10.0}]
    loadsi = [{'type': 'Vertical Load', 'val': -300.0, 'loc': 5.0}]
    #a = SolveProblem(length, supports, loadsi, 1000)

    # Simple Beam
    length = 5.0
    supports = [{'support': 'Roller', 'location': 0.0}, {'support': 'Roller', 'location': 5.0}]
    loadsi = [{'val2': -5000.0, 'loc2': 5.0, 'type': 'Distributed Load', 'val': -5000.0, 'loc': 0.0}]
    a = SolveProblem(length, supports, loadsi, 1000)

    # Simple indeterminate problem 2
    length = 12.0
    supports = [{'support': 'Fixed', 'location': 0.0}, {'support': 'Fixed', 'location': 12.0}]
    loadsi = [{'type': 'Vertical Load', 'val': 6000.0, 'loc': 6.0}]
    #a = SolveProblem(length, supports, loadsi, 1000)


#dummy_run()
# Test cases: fixed with distributed load: PASS
# Test cases: fixed with point load: PASS
# Test cases: fixed with point, distributed and moment: PASS
# Test cases: roller/pinned with point load: PASS
# Test cases: roller/pinned with distributed load: PASS