from tkinter import *
from colorama import Fore
from scipy.integrate import quad, dblquad
import numpy as np
from sympy import Symbol, diff, lambdify
import matplotlib.pyplot as plt
import matplotlib.patches as patches

LOGO = "A:\Programs I\PyCharm Projects\LoadPlotter\I-beam.png.ico"


class SolveProblem:
    def __init__(self, beam, loads, axial_loads, area, I_ar, E_mod, adv_dict, resolution_const=5000):
        """
        format for BEAM: {'Length': 'distance', 'Supports': [['#type', 'location'], ['#type', 'location']]}
        format for LOADS: {'LOAD#': ['type', 'value', 'location']}
        for distributed loads {'LOAD#': ['DISTRIBURED', [val0, val1], [loc0, loc1]]}
        types of loads: moment, force, distributed
        types of supports: fixed, pinned
        """
        print("Starting")
        self.AREA = area
        self.MOMENT_AREA = I_ar
        self.E_mod = E_mod
        if adv_dict:
            self.yNA = adv_dict['yNA']
            self.rectangles = adv_dict['Rectangle']
            self.loc_solve = adv_dict['Loc Solve']

        print("Moment of Area and Young's Modulus:")
        print(self.MOMENT_AREA, E_mod)
        self.MOMENT_SUPPORT_FLAG = False
        self.BEAM = beam
        self.LOADS = loads
        self.AXIAL_LOADS = axial_loads
        print("Beam, Loads and Axial Loads")
        print(self.BEAM)
        print(self.LOADS)
        print(self.AXIAL_LOADS)
        self.axial_loading_flag = False
        if self.AXIAL_LOADS:
            self.axial_loading_flag = True
        length = self.BEAM['Length']
        self.reaction_window = None
        self.resolution = length/resolution_const
        self.X_COOR = np.arange(0, (self.BEAM['Length']), self.resolution)
        self.NUM_REAC_FORCES = 0
        self.NUM_REAC_MOMENTS = 0
        self.SHEAR_FORCE = [0]*len(self.X_COOR)
        self.BENDING_MOMENT = [0]*len(self.X_COOR)
        self.LOAD_SLOPE = [0]*len(self.X_COOR)
        self.LOAD_BEND = [0]*len(self.X_COOR)
        self.AXIAL_FORCE = [0]*len(self.X_COOR)
        self.MOMENT_SUPPORTS = []
        self.FORCE_SUPPORTS = []
        self.NUM_UNKNOWNS = self.get_unknowns()
        self.LHS = []
        self.RHS = []
        self.REACTION_SOLUTIONS = {}
        print('0%')
        print('Forming force balance equation...')
        self.get_force_eq()
        print('10%')
        print('Forming moment balance equation...')
        self.get_moment_eq()
        print('Getting Load Features..')
        self.get_load_features()
        print('20%')
        print('Forming deflection and slope equations...')
        if self.MOMENT_SUPPORT_FLAG:
            self.get_slope_equations()
        self.get_defl_equation()
        print('70%')
        print('Solving equations...')
        print(self.LHS, self.RHS)
        self.solutions = np.linalg.solve(self.LHS, self.RHS)
        print(self.solutions)
        print('90%')
        print('Creating plots...')
        self.correct_shear_moment_bending()
        self.correct_axial_force()
        print('100%')
        print('Done!')
        plt.plot(self.X_COOR, self.SHEAR_FORCE, 'g')
        plt.title('Shear')
        plt.grid()
        plt.ylabel("Shear Force in N")
        plt.xlabel("Distance in m")
        plt.axhline(0, color='black')
        plt.figure()
        print("Now", self.AXIAL_LOADS)
        if self.axial_loading_flag:
            plt.plot(self.X_COOR, self.AXIAL_FORCE, 'm')
            plt.title('Axial Force')
            plt.grid()
            plt.ylabel("Axial Force in N")
            plt.xlabel("Distance in m")
            plt.axhline(0, color='black')
            plt.figure()
        plt.plot(self.X_COOR, self.BENDING_MOMENT, 'r')
        plt.title('Moment')
        plt.grid()
        plt.ylabel("Bending Moment in Nm")
        plt.xlabel("Distance in m")
        plt.axhline(0, color='black')
        self.format_reactions()
        if not (self.MOMENT_AREA == 1 and self.AREA == 1 and self.E_mod == 1):
            for index in range(len(self.LOAD_BEND)):
                self.LOAD_BEND[index] = self.LOAD_BEND[index]/(self.MOMENT_AREA*self.E_mod)
            plt.figure()
            plt.plot(self.X_COOR, self.LOAD_BEND, 'c')
            plt.title('Deflection')
            plt.grid()
            plt.ylabel("Beam Deflection in m")
            plt.xlabel("Distance in m")
            plt.axhline(0, color='black')
            if adv_dict:
                self.get_cross_section_stresses()
            else:
                plt.show()
        else:
            plt.show()

    def get_cross_section_stresses(self):
        location = 0
        for index_flag in range(len(self.X_COOR)):
            if self.loc_solve == self.X_COOR[index_flag]:
                location = index_flag
                break
        M = self.BENDING_MOMENT[location]
        V = self.SHEAR_FORCE[location]
        Fa = self.AXIAL_FORCE[location]
        length_cx = 0
        width_cx = 0
        for id, rect in self.rectangles.items():
            length_cx += rect['h']
            if rect['b'] > width_cx:
                width_cx = rect['b']
        cxN = np.arange(0, length_cx, ((width_cx + length_cx) / 2000))
        cxH = np.arange(0, width_cx, ((width_cx + length_cx) / 2000))
        shear_stress = [[0] * len(cxH)] * len(cxN)
        horz_stress = [[0] * len(cxH)] * len(cxN)

        for index in range(len(cxN)):
            Q = 0
            thick = 0
            position = cxN[index]
            A_dash = 0
            yA_dash = 0
            for ID, rect in self.rectangles.items():
                height = rect['h']
                width = rect['b']
                y_coor = rect['y']

                if position < (y_coor + height) and position >= y_coor:
                    thick = width
                    A_dash += ((position - y_coor) * width)
                    yA_dash += ((position + y_coor) / 2) * (position - y_coor) * width
                    break
                if position >= (y_coor + height):
                    thick = width
                    A_dash += (width * height)
                    yA_dash += (((height / 2) + y_coor) * width * height)
            if A_dash != 0:
                y_dsh = yA_dash / A_dash
                y_dash = self.yNA - y_dsh
            else:
                y_dash = 0
            Q = y_dash * A_dash
            shear_stress[len(cxN)-index-1] = [V * Q / (self.MOMENT_AREA * thick)] * len(cxH)
            horz_stress[len(cxN)-index-1] = [(M * (self.yNA - position) / self.MOMENT_AREA) + Fa / self.AREA] * len(cxH)

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        for ID, rect in self.rectangles.items():
            y_coor = rect['y']
            x_coor = rect['x']
            width = rect['b']
            height = rect['h']
            ax1.add_patch(patches.Rectangle((x_coor, y_coor), width, height, fill=False))
        plt.imshow(horz_stress, cmap='cool', extent=[0, width_cx, 0, length_cx])
        plt.colorbar()
        plt.title("Axial Stress at %sm in Pa" %(str(self.loc_solve)))
        plt.ylabel("Distance in m")
        plt.xlabel("Distance in m")
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, aspect='equal')
        for ID, rect in self.rectangles.items():
            y_coor = rect['y']
            x_coor = rect['x']
            width = rect['b']
            height = rect['h']
            ax2.add_patch(patches.Rectangle((x_coor, y_coor), width, height, fill=False))
        plt.imshow(shear_stress, cmap='cool', extent=[0, width_cx, 0, length_cx])
        plt.colorbar()
        plt.title("Shear Stress at %sm in Pa" % (str(self.loc_solve)))
        plt.ylabel("Distance in m")
        plt.xlabel("Distance in m")
        plt.show()

    def format_reactions(self):
        self.reaction_window = Tk()
        self.reaction_window.iconbitmap(LOGO)
        self.reaction_window.title("Beam.io/Reactions")
        results_str = ""
        for ID, support in self.REACTION_SOLUTIONS.items():
            if "FIXED" in ID:
                addition_str = "For the fixed support at {LOC}m:\nM = {MOMENT}Nm, Ry = {VERT}N, Rx = {HORZ}N\n".format(
                    LOC=str(support["Location"]),
                    MOMENT=str(support["Moment"]),
                    VERT=str(support["Force_y"]),
                    HORZ=str(support["Force_x"])
                )
                results_str += addition_str
            elif "ROLLER" in ID:
                addition_str = "For the roller support at {LOC}m:\nRy = {VERT}N\n".format(
                    LOC=str(support["Location"]),
                    VERT=str(support["Force_y"]),
                )
                results_str += addition_str
            elif "PINNED" in ID:
                addition_str = "For the pinned support at {LOC}m:\nRy = {VERT}N, Rx = {HORZ}N\n".format(
                    LOC=str(support["Location"]),
                    VERT=str(support["Force_y"]),
                    HORZ=str(support["Force_x"])
                )
                results_str += addition_str
        Label(self.reaction_window, text=results_str, fg='white', bg='black').pack()

    def get_axial_reactions(self, left_supp, right_supp, value, location):
        LHS = [[1, 1], [(location-left_supp), (location-right_supp)]]
        RHS = [value, 0]
        left_val, right_val = np.linalg.solve(LHS, RHS)
        return left_val, right_val

    def correct_axial_force(self):
        AXIAL_REACTIONS = []
        AXIAL_REACTIONS_LOC = []
        for support in self.BEAM['Supports']:
            if "ROLLER" not in support[0]:
                AXIAL_REACTIONS.append([support[0], {'Location': support[1], 'Value': 0}])
                AXIAL_REACTIONS_LOC.append(support[1])

        for ident, values in self.AXIAL_LOADS.items():
            value = values[0]
            location = values[1]
            for index in range(len(AXIAL_REACTIONS_LOC)):
                if location == AXIAL_REACTIONS_LOC[index]:
                    for identifier, attributes in AXIAL_REACTIONS:
                        if attributes['Location'] == location:
                            attributes['Value'] += (-1 * value)
                            break
                    continue
                elif location > AXIAL_REACTIONS_LOC[index]:
                    left_supp = AXIAL_REACTIONS_LOC[index]
                    try:
                        right_supp = AXIAL_REACTIONS_LOC[index + 1]
                        if location >= right_supp:
                            continue
                    except IndexError:
                        for identifier, attributes in AXIAL_REACTIONS:
                            if attributes['Location'] == left_supp:
                                attributes['Value'] += (-1 * value)
                                break
                        continue
                    left_reac, right_reac = self.get_axial_reactions(left_supp, right_supp, value, location)
                    for identifier, attributes in AXIAL_REACTIONS:
                        if attributes['Location'] == left_supp:
                            attributes['Value'] += (-1 * left_reac)
                            break
                    for identifier, attributes in AXIAL_REACTIONS:
                        if attributes['Location'] == right_supp:
                            attributes['Value'] += (-1 * right_reac)
                            break
                    continue
                elif location < AXIAL_REACTIONS_LOC[index]:
                    for identifier, attributes in AXIAL_REACTIONS:
                        if attributes['Location'] == AXIAL_REACTIONS[index]:
                            attributes['Value'] += (-1 * value)
                            break
                    continue
        count = 1
        placeholder_axial_loads = {}
        placeholder_axial_loads.update(self.AXIAL_LOADS)
        for identifier, attributes in AXIAL_REACTIONS:
            variable_name = "REAC_LOAD" + str(count)
            addition_list = [attributes['Value'], attributes['Location']]
            placeholder_axial_loads.update({variable_name: addition_list})
            count += 1

        for identifier, attributes in placeholder_axial_loads.items():
            for i in range(len(self.X_COOR)):
                if attributes[1] >= self.X_COOR[i]:
                    self.AXIAL_FORCE[i] += attributes[0]

        for identifier, attributes in placeholder_axial_loads.items():
            for ID, support_location in self.BEAM['Supports']:
                if support_location == attributes[1]:
                    (self.REACTION_SOLUTIONS[ID]).update({'Force_x': attributes[0]})

    def correct_shear_moment_bending(self):
        supports = self.BEAM['Supports']
        x_sym = Symbol('x')
        A = self.solutions[self.NUM_REAC_FORCES+self.NUM_REAC_MOMENTS]
        B = self.solutions[self.NUM_REAC_FORCES+self.NUM_REAC_MOMENTS+1]
        for support in supports:
            if 'FIXED' in support[0]:
                index_m = int((support[0])[0]) + self.NUM_REAC_FORCES - 1
                index_f = self.NUM_REAC_FORCES - int((support[0])[0])
                value_m = self.solutions[index_m]
                value_f = self.solutions[index_f]
                loc = support[1]
                moment = lambda x: self.heavside(x, loc) * (-1*value_m + value_f*(x-loc))
                moment_d = lambda x, y: self.heavside(x, loc) * (-1*value_m + value_f*(x-loc))
                shear = diff((-1*value_m + value_f * (x_sym - loc)), x_sym)
                shear = lambdify(x_sym, shear, 'numpy')
                count = 0
                for x in self.X_COOR:
                    self.SHEAR_FORCE[count] += self.heavside(x, loc) * shear(x)
                    sl = quad(moment, loc, x)
                    defl = dblquad(moment_d, loc, x, lambda x: loc, lambda x: x)
                    self.BENDING_MOMENT[count] += moment(x)
                    self.LOAD_SLOPE[count] += sl[0]
                    self.LOAD_BEND[count] = (self.LOAD_BEND[count] + defl[0])
                    count += 1
                results_addition = {"Location": support[1], "Moment": value_m, "Force_y": value_f}
                self.REACTION_SOLUTIONS.update({support[0]: results_addition})
            else:
                index = int((support[0])[0]) - 1
                loc = support[1]
                val = self.solutions[index]
                moment = lambda x: self.heavside(x, loc) * (val * (x - loc))
                moment_d = lambda x, y: self.heavside(x, loc) * (val * (x - loc))
                shear = diff((val * (x_sym - loc)), x_sym)
                shear = lambdify(x_sym, shear, 'numpy')
                count = 0
                for x in self.X_COOR:
                    self.SHEAR_FORCE[count] += self.heavside(x, loc)*shear(x)
                    self.BENDING_MOMENT[count] += moment(x)
                    sl = quad(moment, loc, x)
                    defl = dblquad(moment_d, loc, x, lambda x: loc, lambda x: x)
                    self.LOAD_SLOPE[count] += sl[0]
                    self.LOAD_BEND[count] = (self.LOAD_BEND[count] + defl[0])
                    count += 1
                results_addition = {"Location": loc, "Force_y": val}
                self.REACTION_SOLUTIONS.update({support[0]: results_addition})
        count = 0
        for x in self.X_COOR:
            self.LOAD_SLOPE[count] += A
            self.LOAD_BEND[count] += (A*x + B)
            count += 1

    def get_defl_equation(self):
        supports = self.BEAM['Supports']
        for location in self.FORCE_SUPPORTS:
            defl_equation = [0]*self.NUM_UNKNOWNS
            for support in supports:
                if 'FIXED' in support[0]:
                    index_m = int((support[0])[0]) + self.NUM_REAC_FORCES - 1
                    index_f = self.NUM_REAC_FORCES - int((support[0])[0])
                    loc = support[1]
                    moment_m = lambda x, y: self.heavside(x, loc)
                    moment_f = lambda x, y: self.heavside(x, loc)*(x-loc)
                    defl_m = dblquad(moment_m, loc, location, lambda x: loc, lambda x: x)
                    defl_f = dblquad(moment_f, loc, location, lambda x: loc, lambda x: x)
                    defl_equation[index_m] = -1*defl_m[0]
                    defl_equation[index_f] = defl_f[0]

                else:
                    index = int((support[0])[0]) - 1
                    loc = support[1]
                    moment = lambda x, y: self.heavside(x, loc)*(x-loc)
                    defl = dblquad(moment, loc, location, lambda x: loc, lambda x: x)
                    defl_equation[index] = defl[0]
            defl_equation[self.NUM_REAC_MOMENTS + self.NUM_REAC_FORCES] = location
            defl_equation[self.NUM_REAC_MOMENTS + self.NUM_REAC_FORCES+1] = 1
            self.LHS.append(defl_equation)
            tmp_index = location * (1/self.resolution)
            tmp_index = int(tmp_index)
            if tmp_index != 0:
                tmp_index -=1
            self.RHS.append(-1 * self.LOAD_BEND[tmp_index])

    def get_slope_equations(self):
            supports = self.BEAM['Supports']
            for location in self.MOMENT_SUPPORTS:
                slope_equation = [0]*self.NUM_UNKNOWNS
                for support in supports:
                    if 'FIXED' in support[0]:
                        index_m = int((support[0])[0]) + self.NUM_REAC_FORCES - 1
                        index_f = self.NUM_REAC_FORCES - int((support[0])[0])
                        loc = support[1]
                        moment_m = lambda x: self.heavside(x, loc)
                        moment_f = lambda x: self.heavside(x, loc)*(x-loc)
                        sl_m = quad(moment_m, loc, location)
                        sl_f = quad(moment_f, loc, location)
                        slope_equation[index_m] = -1*sl_m[0]
                        slope_equation[index_f] = sl_f[0]
                    else:
                        index = int((support[0])[0]) - 1
                        loc = support[1]
                        moment = lambda x: self.heavside(x, loc)*(x-loc)
                        sl = quad(moment, loc, location)
                        slope_equation[index] = sl[0]
                slope_equation[self.NUM_REAC_MOMENTS+self.NUM_REAC_FORCES] = 1
                self.LHS.append(slope_equation)
                tmp_index = int(location*(1/self.resolution))
                if tmp_index !=0:
                    tmp_index -=1
                self.RHS.append(-1*self.LOAD_SLOPE[tmp_index])

    def get_load_features(self):
        x_sym = Symbol('x')
        for load, specs in self.LOADS.items():
            if 'FORCE' in specs[0]:
                val = specs[1]
                loc = specs[2]
                moment = lambda x: self.heavside(x, loc)*(val*(x-loc))
                moment_d = lambda x, y: self.heavside(x, loc)*(val*(x-loc))
                shear = diff((val * (x_sym - loc)), x_sym)
                shear = lambdify(x_sym, shear, 'numpy')
                count = 0
                for x in self.X_COOR:
                    self.SHEAR_FORCE[count] += self.heavside(x, loc)*shear(x)
                    self.BENDING_MOMENT[count] += moment(x)
                    if self.heavside(x, loc) == 1:
                        sl = quad(moment, loc, x)
                        defl = dblquad(moment_d, loc, x, lambda x: loc, lambda x: x)
                        self.LOAD_SLOPE[count] += sl[0]
                        self.LOAD_BEND[count] += defl[0]
                    else:
                        self.LOAD_SLOPE[count] += 0
                        self.LOAD_BEND[count] += 0
                    count += 1

            elif 'MOMENT' in specs[0]:
                val = specs[1]
                loc = specs[2]
                moment = lambda x: self.heavside(x, loc)*-1*val
                moment_d = lambda x, y: self.heavside(x, loc) * -1 * val
                count = 0
                for x in self.X_COOR:
                    self.BENDING_MOMENT[count] += moment(x)
                    if self.heavside(x, loc) == 1:
                        sl = quad(moment, loc, x)
                        defl = dblquad(moment_d, loc, x, lambda x: loc, lambda x: x)
                        self.LOAD_SLOPE[count] += sl[0]
                        self.LOAD_BEND[count] += defl[0]
                    else:
                        self.LOAD_SLOPE[count] += 0
                        self.LOAD_BEND[count] += 0
                    count += 1

            elif 'DISTRIBUTED' in specs[0]:
                val = specs[1]
                loc = specs[2]
                val_0 = float(val[0])
                val_1 = float(val[1])
                loc_0 = float(loc[0])
                loc_1 = float(loc[1])
                slope = (val_1-val_0)/(loc_1-loc_0)
                coeff = val_0 - slope * loc_0
                if slope != 0.0:
                    coeff = float(coeff) + (slope*loc_0)
                else:
                    coeff = float(coeff)
                moment = lambda x: self.heavside(x, loc_0)*((coeff/2)*(x-loc_0)**2) - self.heavside(x, loc_1) * \
                                                                                      ((coeff/2)*(x-loc_1)**2)
                moment_d = lambda x, y: self.heavside(x, loc_0) * ((coeff / 2) * (x - loc_0) ** 2) - self.heavside(x,
                                                                                                loc_1) * ((coeff / 2) *
                                                                                                (x - loc_1) ** 2)
                shear_1 = diff((coeff/2)*((x_sym-loc_0)**2))
                shear_1 = lambdify(x_sym, shear_1, 'numpy')
                shear_2 = diff((coeff / 2) * ((x_sym - loc_1) ** 2))
                shear_2 = lambdify(x_sym, shear_2, 'numpy')
                count = 0
                for x in self.X_COOR:
                    sl = quad(moment, loc_0, x)
                    defl = dblquad(moment_d, loc_0, x, lambda x: loc_0, lambda x: x)
                    self.LOAD_SLOPE[count] += sl[0]
                    self.LOAD_BEND[count] += defl[0]
                    self.SHEAR_FORCE[count] += (self.heavside(x, loc_0) * shear_1(x)) - (self.heavside(x, loc_1) * shear_2(x))
                    self.BENDING_MOMENT[count] += moment(x)
                    count += 1

                if slope != 0.0:
                    moment = lambda x: self.heavside(x, loc_0) * ((slope / 6) * (x - loc_0) ** 3) - \
                                       (self.heavside(x,loc_1) * ((slope / 6) * (x - loc_1) ** 3)) - \
                                       (self.heavside(x, loc_1) * (((slope*(loc_1-loc_0))/2) * (x - loc_1)**2))
                    moment_d = lambda x, y: self.heavside(x, loc_0) * ((slope / 6) * (x - loc_0) ** 3) - \
                                       (self.heavside(x, loc_1) * ((slope / 6) * (x - loc_1) ** 3)) - \
                                       (self.heavside(x, loc_1) * (((slope * (loc_1 - loc_0)) / 2) * (x - loc_1) ** 2))
                    shear_1 = diff((slope / 6) * ((x_sym - loc_0) ** 3))
                    shear_1 = lambdify(x_sym, shear_1, 'numpy')
                    shear_2 = diff((slope / 6) * ((x_sym - loc_1) ** 3))
                    shear_2 = lambdify(x_sym, shear_2, 'numpy')
                    shear_3 = diff(((slope*(loc_1-loc_0))/2) * ((x_sym - loc_1)**2))
                    shear_3 = lambdify(x_sym, shear_3, 'numpy')
                    count = 0
                    for x in self.X_COOR:
                        sl = quad(moment, loc_0, x)
                        defl = dblquad(moment_d, loc_0, x, lambda x: loc_0, lambda x: x)
                        self.LOAD_SLOPE[count] += sl[0]
                        self.LOAD_BEND[count] += defl[0]
                        self.SHEAR_FORCE[count] += (self.heavside(x, loc_0) * shear_1(x)) - (self.heavside(x, loc_1) *
                                                                    shear_2(x)) - (self.heavside(x, loc_1) * shear_3(x))
                        self.BENDING_MOMENT[count] += moment(x)
                        count += 1

    def get_moment_eq(self):
        MOMENT_EQUATION = [0]*self.NUM_UNKNOWNS
        sum_external_moment = 0
        for load, specs in self.LOADS.items():
            if specs[0] == 'MOMENT':
                sum_external_moment += specs[1]
            if specs[0] == 'DISTRIBUTED':
                eq_val, eq_loc = self.get_equivalent_PL(specs)
                sum_external_moment += eq_val*eq_loc
            elif specs[0] == 'FORCE':
                sum_external_moment += specs[1]*specs[2]
        for support in self.BEAM['Supports']:
            if 'FIXED' in support[0]:
                index_m = int((support[0])[0]) - 1
                index_f = self.NUM_REAC_FORCES - int((support[0])[0])
                index_m += self.NUM_REAC_FORCES
                MOMENT_EQUATION[index_m] = 1
                MOMENT_EQUATION[index_f] = support[1]
            else:
                index = int((support[0])[0]) - 1
                MOMENT_EQUATION[index] = support[1]
        self.LHS.append(MOMENT_EQUATION)
        self.RHS.append(-1*sum_external_moment)

    def get_force_eq(self):
        FORCE_EQUATION = [0]*self.NUM_UNKNOWNS
        sum_external_force = 0
        for load, specs in self.LOADS.items():
            if specs[0] == 'DISTRIBUTED':
                eq_val, eq_loc = self.get_equivalent_PL(specs)
                sum_external_force += eq_val
            elif specs[0] == 'FORCE':
                 sum_external_force += specs[1]
        for i in range(self.NUM_REAC_FORCES):
            FORCE_EQUATION[i] = 1
        self.LHS.append(FORCE_EQUATION)
        self.RHS.append(-1*sum_external_force)

    def get_equivalent_PL(self, specs):
        slope = (((specs[1])[1])-((specs[1])[0]))/(((specs[2])[1])-((specs[2])[0]))
        c = (specs[1])[0] - slope * ((specs[2])[0])
        avg_func = lambda x: x*(slope*x + c)
        load_func = lambda x: slope*x + c
        I = quad(load_func, ((specs[2])[0]), ((specs[2])[1]))
        I2 = quad(avg_func, ((specs[2])[0]), ((specs[2])[1]))
        return I[0], (I2[0]) / (I[0])

    def get_unknowns(self):
        supports = self.BEAM['Supports']
        unknowns = 0
        for support in supports:
            if 'FIXED' in support[0]:
                self.MOMENT_SUPPORT_FLAG = True
                self.MOMENT_SUPPORTS.append(support[1])
                self.FORCE_SUPPORTS.append(support[1])
                self.NUM_REAC_MOMENTS += 1
                self.NUM_REAC_FORCES += 1
                unknowns += 2
            elif 'ROLLER' in support[0]:
                self.NUM_REAC_FORCES += 1
                self.FORCE_SUPPORTS.append(support[1])
                unknowns += 1
            elif 'PINNED' in support[0]:
                self.NUM_REAC_FORCES += 1
                self.FORCE_SUPPORTS.append(support[1])
                unknowns += 1
        unknowns += 2
        return unknowns

    def heavside(self, x, a, sign=True):
        if sign:
            if (x-a) >= 0:
                return 1
        if not sign:
            if (x-a) > 0:
                return 1
        return 0


class SetupProblem:
    def __init__(self):
        self.tk = Tk()
        self.rectangles = None
        self.rectangle_log = ""
        self.rectangle_label = None
        self.rec_count = 1
        self.A = 1
        self.I = 1
        self.E = 1
        self.loc_solver = 0
        self.yNA = 0
        self.adv_dict = None
        self.load_window = None
        self.load_types = None
        self.CUST_RECEIVED = False
        self.ADV_RECEIVED = False
        self.tk.iconbitmap(LOGO)
        self.tk.title("Beam.io")
        self.BEAM = {'Length': 0, 'Supports': []}
        self.LOADS = {}
        self.AXIAL_LOADS = {}
        self.num_loads = 1
        self.num_axial_loads = 1
        self.fixed_sup_count = 1
        self.vertical_sup_count = 1
        self.load_counts = ""
        self.supp_counts = ""
        self.status_label = Label(self.tk, text='Status: Setup the beam...')
        self.status_label.grid(row=0, columns=1)
        Label(self.tk, relief=FLAT, width=25).grid(row=0, column=2)
        self.help_button = Button(self.tk, relief=GROOVE, bg='black', fg='white', text='Help', width=20,
                                  command=lambda: self.get_help())
        self.help_button.grid(row=0, column=3)
        self.beam_len_label = Label(self.tk, text="Enter the Length of the beam in metres: ")
        self.beam_len_label.grid(row=1, column=0)
        self.beam_len_entry = Entry(self.tk)
        self.beam_len_entry.grid(row=1, column=1)
        Label(self.tk, relief=FLAT, width=25).grid(row=2, column=0)
        self.supp_clear = Button(self.tk, relief=GROOVE, text='Clear Supports', width=20, fg='red',
                                 command = lambda: self.clear_supp())
        self.supp_clear.grid(row=3, column=0)
        Label(self.tk, width=25).grid(row=3, column=1)
        self.load_clear = Button(self.tk, relief=GROOVE, text='Clear Loads', width=20, fg='red',
                                 command= lambda : self.clear_loads())
        self.load_clear.grid(row=3, column=2)
        self.add_supp = Button(self.tk, relief=RAISED, text='Add a Support', command= lambda: self.get_supp())
        self.add_supp.grid(row=4, column=0)
        self.add_load = Button(self.tk, relief=RAISED, text='Add a Load', command= lambda:self.get_load())
        self.add_load.grid(row=4, column=2)
        self.support_log = Label(self.tk, text=self.supp_counts)
        self.support_log.grid(row=5, column=0)
        self.load_log = Label(self.tk, text=self.load_counts)
        self.load_log.grid(row=5, column=2)

        Label(self.tk, text="For Advanced Problems:").grid(row=6, column=0)
        self.advanced_options_btn = Button(self.tk, text="Enter Advanced Prop.s", relief=RAISED,
                                           command=lambda: self.get_adv_prop())
        self.advanced_options_btn.grid(row=7, column=0)
        Label(self.tk, text='OR').grid(row=7, column=1)
        self.custom_beam = Button(self.tk, text="Enter Beam X-Section", relief=RAISED, command=lambda: self.get_Xsec())
        self.custom_beam.grid(row=7, column=2)

        self.solve_btn = Button(self.tk, text="Solve Beam", relief=RAISED, bg='green', command=lambda: self.solver())
        self.solve_btn.grid(row=8, column=0)

        self.tk.mainloop()

    def solver(self):
        try:
            length = self.beam_len_entry.get()
            if length == "":
                self.status_label.configure(text="Beam length not defined...", fg='red')
                raise ValueError
            length = float(length)
            if length <= 0:
                self.status_label.configure(text="Invalid Beam Length...", fg='red')
                raise ValueError
            support_bool = False
            semi_support_bool = False
            for beam_el in self.BEAM['Supports']:
                supp_el = beam_el[0]
                if "FIXED" in supp_el:
                    support_bool = True
                    break
                else:
                    if semi_support_bool:
                        support_bool = True
                        break
                    semi_support_bool = True

            if not support_bool:
                self.status_label.configure(text="Beam not fully supproted...", fg='red')
                raise ValueError
        except ValueError:
            print(Fore.RED + "Definition Error Preventing Solution!" + Fore.RESET)
            return
        self.BEAM['Length'] = length
        if not self.ADV_RECEIVED:
            self.A = 1
            self.I = 1
            self.E = 1
            self.adv_dict = None
        else:
            if self.adv_dict:
                check = self.adv_dict['Loc Solve']
                if float(check) > float(self.BEAM['Length']):
                    self.status_label.configure(text="Your location for evaluating\nstresses "
                                                     "is beyond the beam...", fg='red')
                    return
        self.status_label.configure(text="Solving the beam!", fg='green')
        SolveProblem(self.BEAM, self.LOADS, self.AXIAL_LOADS, self.A, self.I, self.E, self.adv_dict)

    def get_Xsec(self):
        custom_xsec_win = Tk()
        custom_xsec_win.iconbitmap(LOGO)
        custom_xsec_win.title("Beam.io/Custom Beam")
        Label(custom_xsec_win, text="Enter E in GPa:").grid(row=0, column=0)
        self.rectangles = {}
        E_entry = Entry(custom_xsec_win)
        E_entry.grid(row=0, column=1)
        Label(custom_xsec_win, text="Location along the beam length where\n"
                                    "you want stress plots (in m):").grid(row=1, column=0)
        loc_entry = Entry(custom_xsec_win)
        loc_entry.grid(row=1, column=1)
        Label(custom_xsec_win, text="Construct your cross section with rectangles.\n"\
                                    "Ensure there is no overlap!").grid(row=2, column=0)
        clear_x_sec_btn = Button(custom_xsec_win, text="Clear Cross Section", relief=RAISED,
                                 command = lambda: self.clear_x_section())
        clear_x_sec_btn.grid(row=3, column=0)
        self.rectangle_label = Label(custom_xsec_win, text=self.rectangle_log)
        self.rectangle_label.grid(row=4, column=0)
        add_rect_btn = Button(custom_xsec_win, text="Add Rectangle", relief=RAISED,
                              command = lambda: self.add_x_section())
        add_rect_btn.grid(row=5, column=0)
        view_cross_sec = Button(custom_xsec_win, text="View Cross Section", relief=RAISED,
                                command= lambda: self.show_x_sec())
        view_cross_sec.grid(row=5, column=1)
        submit_cross_sec = Button(custom_xsec_win, text="Submit", relief=RAISED,
                                command=lambda: self.submit_cross_section(custom_xsec_win, E_entry, loc_entry))
        submit_cross_sec.grid(row=6, column=0)

    def submit_cross_section(self, win, E_entry, loc_entry):
        if self.rectangles:
            try:
                self.E = E_entry.get()
                self.loc_solver = loc_entry.get()
                self.E = float(self.E)*(10**9)
                self.loc_solver = float(self.loc_solver)
                if self.E <= 0:
                    self.status_label.configure(text="E must be positive..", fg = 'red')
                    raise ValueError
                if self.loc_solver <= 0:
                    self.status_label.configure(text="Location for stress solution must be positive..", fg='red')
                    raise ValueError
            except ValueError:
                return
            A = 0
            yA = 0
            I = 0
            for ID, rect in self.rectangles.items():
                y_coor = rect['y']
                width = rect['b']
                height = rect['h']
                A += width * height
                yA += ((height / 2) + y_coor) * width * height

            yNA = yA / A

            for ID, rect in self.rectangles.items():
                y_coor = rect['y']
                width = rect['b']
                height = rect['h']

                I1 = (1 / 12) * width * (height ** 3)
                I2 = width * height * (yNA - y_coor - (height / 2)) ** 2
                I += (I1 + I2)

            self.A = A
            self.yNA = yNA
            self.I = I
            self.adv_dict = {"Rectangle": self.rectangles, "yNA": self.yNA, "Loc Solve": self.loc_solver}
            self.ADV_RECEIVED = True
            self.rectangle_label.configure(text="")
            self.rectangle_log = ""
            win.destroy()
            self.status_label.configure(text="Received custom cross section and properties")
            self.rectangles = {}
            self.yNA = 0
            self.loc_solver = 0

    def show_x_sec(self):
        if self.rectangles:
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            for ID, rect in self.rectangles.items():
                y_coor = rect['y']
                x_coor = rect['x']
                width = rect['b']
                height = rect['h']
                ax1.add_patch(patches.Rectangle((x_coor, y_coor), width, height, fill=False))
            plt.xlabel("Distance in m")
            plt.ylabel("Distance in m")
            plt.title("Beam Cross Section")
            plt.show()

    def clear_x_section(self):
        self.rectangles = {}
        self.rec_count = 1
        self.rectangle_log = ""
        self.rectangle_label.configure(text=self.rectangle_log)
        return

    def add_x_section(self):
        add_rect_window = Tk()
        add_rect_window.iconbitmap(LOGO)
        add_rect_window.title("Beam.io/Add Rectangle")
        Label(add_rect_window, text='Enter the x coordinate of the bottom\n'
                                    'left corner of the rectangle in cm:').grid(row=0, column=0)
        x_loc_entry = Entry(add_rect_window)
        x_loc_entry.grid(row=0, column=1)
        Label(add_rect_window, text='Enter the y coordinate of the bottom\n'
                                    'left corner of the rectangle in cm:').grid(row=1, column=0)
        y_loc_entry = Entry(add_rect_window)
        y_loc_entry.grid(row=1, column=1)
        Label(add_rect_window, text='Enter the width of the rectangle in cm:').grid(row=2, column=0)
        width_entry = Entry(add_rect_window)
        width_entry.grid(row=2, column=1)
        Label(add_rect_window, text='Enter the height of the rectangle in cm:').grid(row=3, column=0)
        height_entry = Entry(add_rect_window)
        height_entry.grid(row=3, column=1)
        submit_rect = Button(add_rect_window, text="Submit",
                             command=lambda: self.accept_rect(add_rect_window, x_loc_entry, y_loc_entry, width_entry,
                                                              height_entry))
        submit_rect.grid(row=4, column=0)

    def accept_rect(self, win, x_loc, y_loc, width, height):
        try:
            x_loc = x_loc.get()
            y_loc = y_loc.get()
            width = width.get()
            height = height.get()
            x_loc = float(x_loc)
            y_loc = float(y_loc)
            width = float(width)
            height = float(height)
        except ValueError:
            self.status_label.configrue(text='Invalid Rectangle Geometery...', fg='red')
            return
        self.rectangle_log += "\nRectangle at ({X}, {Y})cm\nwith width {WID}cm and height {HEI}cm".format(
            X=str(x_loc),
            Y=str(y_loc),
            WID=str(width),
            HEI=str(height)
        )
        x_loc = float(x_loc) / 100
        y_loc = float(y_loc) / 100
        width = float(width) / 100
        height = float(height) / 100
        ID_rec = "rec"+str(self.rec_count)
        self.rec_count += 1
        addition = {ID_rec: {'y': y_loc, 'b': width, 'h': height, 'x': x_loc}}
        self.rectangles.update(addition)
        self.rectangle_label.configure(text=self.rectangle_log)
        win.destroy()


    def get_adv_prop(self):
        advanced_prop_window = Tk()
        advanced_prop_window.iconbitmap(LOGO)
        advanced_prop_window.title("Beam.io/Advanced Properties")
        Label(advanced_prop_window, text="Enter E in GPa:").grid(row=0, column=0)
        E_entry = Entry(advanced_prop_window)
        E_entry.grid(row=0, column=1)
        Label(advanced_prop_window, text="Enter A in cm^2:").grid(row=1, column=0)
        A_entry = Entry(advanced_prop_window)
        A_entry.grid(row=1, column=1)
        Label(advanced_prop_window, text="Enter I in cm^4:").grid(row=2, column=0)
        I_entry = Entry(advanced_prop_window)
        I_entry.grid(row=2, column=1)
        submit_adv_btn = Button(advanced_prop_window, text="Submit", command=lambda: self.submit_adv_prop(E_entry,
                                                                                               A_entry, I_entry,
                                                                                               advanced_prop_window))
        submit_adv_btn.grid(row=3, column=0)
        advanced_prop_window.mainloop()

    def submit_adv_prop(self, E_entry, A_entry, I_entry, win):
        try:
            self.E = E_entry.get()
            self.E = float(self.E)
            if self.E <=0:
                self.status_label.configure(text="E must be a positive number!", fg='red')
                raise ValueError
            self.I = I_entry.get()
            self.I = float(self.I)
            if self.I <= 0:
                self.status_label.configure(text="I must be a positive number!", fg='red')
                raise ValueError
            self.A = A_entry.get()
            self.A = float(self.A)
            if self.A <= 0:
                self.status_label.configure(text="A must be a positive number!", fg='red')
                raise ValueError
        except ValueError:
            print("Error in receiving advanced beam properties...")
            return
        self.E = self.E * (10**9)
        self.I = self.I * (10**-8)
        self.A = self.A * (10**-4)
        self.CUST_RECEIVED = False
        self.adv_dict = None
        self.ADV_RECEIVED = True
        win.destroy()

    def clear_loads(self):
        self.LOADS = {}
        self.AXIAL_LOADS = {}
        self.num_loads = 1
        self.num_axial_loads = 1
        self.load_counts = ""
        self.load_log.configure(text=self.load_counts)

    def accept_load1(self, *args):
        load_type = self.load_types.get()
        Label(self.load_window, text=load_type, width=20, relief=FLAT, bg='black', fg = 'white').grid(row=0, column=0)
        if load_type == 'DISTRIBUTED':
            Label(self.load_window, text="Enter the start location in m: ").grid(row=1, column=0)
            start_loc = Entry(self.load_window)
            start_loc.grid(row=1, column=1)
            Label(self.load_window, text="Enter the end location in m: ").grid(row=2, column=0)
            end_loc = Entry(self.load_window)
            end_loc.grid(row=2, column=1)
            Label(self.load_window, text="Enter the start value in N: ").grid(row=3, column=0)
            start_val = Entry(self.load_window)
            start_val.grid(row=3, column=1)
            Label(self.load_window, text="Enter the end value in N: ").grid(row=4, column=0)
            end_val = Entry(self.load_window)
            end_val.grid(row=4, column=1)
            submit_btm = Button(self.load_window, text="Submit Load", command = lambda :self.accept_load2(
                self.load_window, load_type, start_val, start_loc, end_val, end_loc))
            submit_btm.grid(row=5, column=0)
        elif load_type == 'MOMENT':
            Label(self.load_window, text="Enter the moment location in m: ").grid(row=1, column=0)
            start_loc = Entry(self.load_window)
            start_loc.grid(row=1, column=1)
            Label(self.load_window, text="Enter the moment value in Nm: ").grid(row=2, column=0)
            start_val = Entry(self.load_window)
            start_val.grid(row=2, column=1)
            submit_btm = Button(self.load_window, text="Submit Load", command=lambda: self.accept_load2(
                self.load_window, load_type, start_val, start_loc))
            submit_btm.grid(row=5, column=0)

        elif load_type == 'FORCE':
            Label(self.load_window, text="Enter the force location in m: ").grid(row=1, column=0)
            start_loc = Entry(self.load_window)
            start_loc.grid(row=1, column=1)
            Label(self.load_window, text="Enter the force value in N: ").grid(row=2, column=0)
            start_val = Entry(self.load_window)
            start_val.grid(row=2, column=1)
            submit_btm = Button(self.load_window, text="Submit Load", command=lambda: self.accept_load2(
                self.load_window, load_type, start_val, start_loc))
            submit_btm.grid(row=5, column=0)

        elif load_type == 'AXIAL':
            Label(self.load_window, text="Enter the force location in m: ").grid(row=1, column=0)
            position = Entry(self.load_window)
            position.grid(row=1, column=1)
            Label(self.load_window, text="Enter the force value in N: ").grid(row=2, column=0)
            value = Entry(self.load_window)
            value.grid(row=2, column=1)
            submit_btm = Button(self.load_window, text="Submit Load", command=lambda: self.accept_load3(
                self.load_window, load_type, value, position))
            submit_btm.grid(row=5, column=0)

    def get_load(self):
        self.load_window = Tk()
        self.load_window.title("Beam.io/Load")
        self.load_window.iconbitmap(LOGO)
        self.load_types = StringVar(self.load_window)
        self.load_types.set("Select the Load Type")
        load_options = OptionMenu(self.load_window, self.load_types, "DISTRIBUTED", "FORCE", "MOMENT", "AXIAL")
        load_options.grid(row=0, column=0)
        self.load_types.trace("w", self.accept_load1)
        self.load_window.mainloop()

    def accept_load2(self, window, load_type, start_val, start_loc, end_val=None, end_loc=None):
        addition = []
        addition.append(str(load_type))
        length = self.beam_len_entry.get()
        try:
            if length == "":
                self.status_label.configure(text="Beam length not defined...", fg='red')
                raise ValueError
            length = float(length)
            if length <= 0:
                self.status_label.configure(text="Invalid Beam Length...", fg='red')
                raise ValueError
            start_loc = start_loc.get()
            start_val = start_val.get()
            start_val = float(start_val)
            start_loc = float(start_loc)
            if end_val and end_loc:
                end_val = end_val.get()
                end_loc = end_loc.get()
                end_val = float(end_val)
                end_loc = float(end_loc)
            if (start_loc<0 or start_loc>length):
                self.status_label.configure(text="Loads must be on the beam...", fg='red')
                raise ValueError
            if end_loc != None:
                if (end_loc<0 or end_loc>length):
                    self.status_label.configure(text="Loads must be on the beam...", fg='red')
                    raise ValueError
            if end_val and end_loc:
                addition.append([start_val, end_val])
                addition.append([start_loc, end_loc])
                self.load_counts += "\n"+str(load_type)+" at " + str(start_loc) + " to " + str(end_loc) + " with " \
                                                                                                    "" + str(
                    start_val) + " to " + str(end_val) + " in N"
                self.load_log.configure(text=self.load_counts)
                load_num = "LOAD" + str(self.num_loads)
                self.num_loads += 1
                self.LOADS.update({load_num:addition})
                window.destroy()
            else:
                addition.append(start_val)
                addition.append(start_loc)
                self.load_counts += "\n" + str(load_type) + " at " + str(start_loc) + "m and " + str(start_val)
                self.load_log.configure(text=self.load_counts)
                load_num = "LOAD" + str(self.num_loads)
                self.num_loads += 1
                self.LOADS.update({load_num: addition})
                window.destroy()
        except ValueError:
            print(Fore.RED + 'Invalid Loads Parameters...' + Fore.RESET)
            window.destroy()
            return

    def accept_load3(self, window, load_type, value, position):
        addition = []
        length = self.beam_len_entry.get()
        try:
            if length == "":
                self.status_label.configure(text="Beam length not defined...", fg='red')
                raise ValueError
            length = float(length)
            if length <= 0:
                self.status_label.configure(text="Invalid Beam Length...", fg='red')
                raise ValueError
            value = value.get()
            position = position.get()
            value = float(value)
            position = float(position)
            if (position < 0 or position > length):
                self.status_label.configure(text="Loads must be on the beam...", fg='red')
                raise ValueError
            addition.append(value)
            addition.append(position)
            self.load_counts += "\n" + str(load_type) + " at " + str(position) + "m and " + str(value)
            self.load_log.configure(text=self.load_counts)
            load_num = "LOAD" + str(self.num_axial_loads)
            self.num_axial_loads += 1
            self.AXIAL_LOADS.update({load_num: addition})
            window.destroy()
        except ValueError:
            print(Fore.RED + 'Invalid Loads Parameters...' + Fore.RESET)
            window.destroy()
            return


    def clear_supp(self):
        self.BEAM['Supports'] = []
        self.supp_counts = ""
        self.support_log.configure(text="")
        self.fixed_sup_count = 1
        self.vertical_sup_count = 1

    def get_supp(self):
        support_window = Tk()
        support_window.title("Beam.io/Support")
        support_window.iconbitmap(LOGO)
        supp_types = StringVar(support_window)
        supp_types.set("Select the Type of Support")
        support_options = OptionMenu(support_window, supp_types, "FIXED", "ROLLER", "PINNED")
        support_options.grid(row=0, column=0)
        Label(support_window, text="Enter the location of the support in m: ").grid(row=1, column=0)
        location_entry = Entry(support_window)
        location_entry.grid(row=1, column=1)
        submit_support = Button(support_window, relief=RAISED, text="Submit Support",
                                command = lambda: self.accept_supp(supp_types, location_entry, support_window))
        submit_support.grid(row=2, column=0)
        support_window.mainloop()

    def accept_supp(self, type, loc, win):
        try:
            length = self.beam_len_entry.get()
            if length == "":
                self.status_label.configure(text="Beam length not defined...", fg='red')
                raise ValueError
            length = float(length)
            if length <= 0:
                self.status_label.configure(text="Invalid Beam Length...", fg='red')
                raise ValueError
            sup_type = type.get()
            sup_loc = loc.get()
            sup_loc = float(sup_loc)
            if sup_loc < 0 or sup_loc > length:
                self.status_label.configure(text="Support must be on the beam...", fg='red')
                raise ValueError
            if sup_type == "FIXED":
                ident = str(self.fixed_sup_count)+sup_type
                addition = [ident, sup_loc]
                self.BEAM['Supports'].append(addition)
                self.fixed_sup_count += 1
                self.supp_counts += "\n" + str(sup_type) + " at " + str(sup_loc) + "m"
                self.support_log.configure(text=self.supp_counts)
            elif sup_type == "PINNED":
                ident = str(self.vertical_sup_count)+sup_type
                addition = [ident, sup_loc]
                self.BEAM['Supports'].append(addition)
                self.vertical_sup_count += 1
                self.supp_counts += "\n" + str(sup_type) + " at " + str(sup_loc) + "m"
                self.support_log.configure(text=self.supp_counts)
            elif sup_type == "ROLLER":
                ident = str(self.vertical_sup_count)+sup_type
                addition = [ident, sup_loc]
                self.BEAM['Supports'].append(addition)
                self.vertical_sup_count += 1
                self.supp_counts += "\n" + str(sup_type) + " at " + str(sup_loc) + "m"
                self.support_log.configure(text=self.supp_counts)
            self.status_label.configure(text="Setting up the beam...", fg='black')
            win.destroy()
        except ValueError:
            print(Fore.RED + "Error in setting up beam..." + Fore.RESET)
            win.destroy()
            return

    def get_help(self):
        help_window = Tk()
        help_window.iconbitmap(LOGO)
        help_window.title("Beam.io/help")
        help_text = 'Beam.io is a python based software for solving shear force, bending moment and deflection for ' \
                    'simply supported, cantilevered and statically indeterminate problems.\nDefine your loads and ' \
                    'supports by clicking the Add a Support and Add a Load buttons.\nClick Clear Loads and Clear' \
                    ' Supports if you wish to reset either of the beam parameters.\nClick Solve Beam when you want to' \
                    ' solve the problem.'
        help_label = Label(help_window, text=help_text)
        help_label.pack()
        help_window.mainloop()

if __name__ == '__main__':
    beam_problem = SetupProblem()
