"""
author: @Udayan Sarin
Develops an interface for beam_io. The given script accepts user inputs for beam geometry, material and dimensions.
Contains a tab to present plots, stresses and load result information
"""
import os
from tkinter import *
from tkinter import ttk
from PIL import ImageTk, Image

import units
import loads
import beam_solver as bs
import stress_solver as ss

import matplotlib
import matplotlib.patches as patches
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
matplotlib.use('TkAgg')


class GUIPlotter:
    """
    Generic class to create a matplotlib tkinter plot using canvas tools
    """
    @staticmethod
    def make_plot(_x, _y, _window, x_axis=None, y_axis=None, plot_title=None):
        """
        create a matplotlib plot instance in a wkinter frame/window
        :param _x: x list to plot
        :param _y: y list to plot
        :param _window: tkinter widget to develop the plot in
        :param x_axis: x axis label
        :param y_axis: y axis label
        :param plot_title: plot title
        :return:
        """
        fig = Figure(figsize=(9, 4))
        a = fig.add_subplot(111)
        a.set_title(plot_title)
        a.set_xlabel(x_axis)
        a.set_ylabel(y_axis)
        a.grid()
        a.plot(_x, _y, color='red')

        canvas = FigureCanvasTkAgg(fig, master=_window)
        return canvas

    @staticmethod
    def plot_cross_section(_rectangles, _window, x_axis=None, y_axis=None, plot_title=None, stress_val=None,
                           cross_section_range=None):
        if not stress_val:
            font_ = 8
        else:
            font_ = 12
        max_side = 0
        if not stress_val:
            fig = Figure(figsize=(3, 3))
        else:
            fig = Figure(figsize=(9, 4))
        a = fig.add_subplot(111)
        for rect in _rectangles:
            y_coor = rect['y_loc']
            x_coor = rect['x_loc']
            width = rect['b']
            height = rect['h']
            if (y_coor + height) > max_side:
                max_side = y_coor + 1.5 * height
            if (x_coor + width) > max_side:
                max_side = x_coor + 1.5 * width
            a.add_patch(patches.Rectangle((x_coor, y_coor), width, height, fill=False, edgecolor='r'))
        if stress_val:
            img = a.imshow(stress_val, cmap='cool', extent=[0, cross_section_range[0], 0, cross_section_range[1]])
            fig.colorbar(img, ax=a)
        a.set_title(plot_title, fontsize=font_)
        a.set_xlabel(x_axis, fontsize=font_-2)
        a.set_ylabel(y_axis, fontsize=font_-2)
        a.tick_params(axis='both', which='major', labelsize=font_-2)
        a.tick_params(axis='both', which='minor', labelsize=font_-2)
        a.set_xlim(0, max_side)
        a.set_ylim(0, max_side)
        a.grid()
        canvas = FigureCanvasTkAgg(fig, master=_window)
        return canvas


class Interface:
    """
    Acts as the central class for controlling the interface for beam_io.
    """

    def __init__(self):
        self._resolution = 1000
        self._supports = []
        self._beam_defined = False
        self._loads = []
        self._rectangles = None
        self._cross_section_unit = None
        self._advanced_properties = {}
        self._len_units = None
        self._force_units = None
        self._beam_len = None
        self._stress_units = None
        self._logo = os.path.join(os.getcwd(), 'I-beam.png.ico')
        self._units = units.UnitControl()
        self._setup = loads.SetupControl()
        self._results_dict = None
        self._widgets = {}
        self._root = Tk()
        self._root.title("beam.io")
        self._root.iconbitmap(self._logo)
        self._root.geometry('900x800')
        self._root.resizable(height=False, width=False)
        self._populate_window()
        self._root.mainloop()

    def _solve_beam(self):
        suff_support = False
        no_roller = True
        lc = []
        try:
            if len(self._loads) == 0:
                raise RuntimeError("No loading provided\n")
            for index, support in enumerate(self._supports):
                if support['location'] in lc:
                    raise RuntimeError("Two supports share same location\n")
                lc.append(support['location'])
                if support['support'] == "Fixed":
                    suff_support = True
                    break
                if support['support'] == 'Roller':
                    no_roller = False
                if not no_roller:
                    if support['support'] == 'Pinned':
                        suff_support = True
                        break
            if not suff_support and (len(self._supports) < 2):
                raise RuntimeError("Insufficient defining features to solve problem\n")
        except RuntimeError as e:
            self.output(e)
            return
        _sol = bs.SolveProblem(self._beam_len, self._supports, self._loads, self._resolution)
        self._results_dict = _sol.get_results()
        self._setup_results_tab()
        self._setup_advanced_tab()
        return

    def _setup_advanced_tab(self):
        """
        setup a tkinter window tab to accept and display advanced problem information - input: material properties
        output: deflection, stresses
        :return:
        """
        self._rectangles = []
        self._widgets['TabControl'].add(self._widgets['AdvancedTab'], text='Advanced Properties')
        self._widgets['MechInputFrame'] = LabelFrame(self._widgets['AdvancedTab'])
        self._widgets['MechInputFrame'].pack()
        self._widgets['MechOutputFrame'] = LabelFrame(self._widgets['AdvancedTab'])
        self._widgets['MechOutputFrame'].pack()
        self._widgets['ExplicitPropsSection'] = LabelFrame(self._widgets['MechInputFrame'])
        self._widgets['ImplicitPropsFrame'] = LabelFrame(self._widgets['MechInputFrame'])
        self._widgets['ImplicitCrossSection'] = LabelFrame(self._widgets['MechInputFrame'])
        self._widgets['ExplicitPropsSection'].pack(side=LEFT)
        self._widgets['ImplicitCrossSection'].pack(side=RIGHT)
        self._widgets['ImplicitPropsFrame'].pack(side=RIGHT)
        self._widgets['ExplicitPropsFrame'] = LabelFrame(self._widgets['ExplicitPropsSection'])
        self._widgets['ExplicitPropsFrame'].pack()
        self._widgets['AdvancedStatus'] = Text(self._widgets['ExplicitPropsSection'], height=6, width=30)
        self._widgets['AdvancedStatus'].pack()
        self.output('waiting..', adv_window=True)
        self._defined_implicit_frame()

        self._widgets['AdvancedPlotFrame'] = LabelFrame(self._widgets['MechOutputFrame'])
        self._widgets['AdvancedPlotFrame'].pack()

        Label(self._widgets['ExplicitPropsFrame'], text="Enter Young's Modulus").grid(row=0, column=0)
        self._widgets['YoungsEntry'] = Entry(self._widgets['ExplicitPropsFrame'], width=5)
        self._widgets['YoungsUnit'] = StringVar(self._widgets['ExplicitPropsFrame'])
        self._widgets['YoungsUnit'].set('--')
        self._widgets['YoungsMenu'] = OptionMenu(self._widgets['ExplicitPropsFrame'], self._widgets['YoungsUnit'],
                                                 *['GPa', 'Mpsi'])
        self._widgets['YoungsEntry'].grid(row=0, column=1)
        self._widgets['YoungsMenu'].grid(row=0, column=2)
        Label(self._widgets['ExplicitPropsFrame'], text="Enter Moment of Area").grid(row=1, column=0)
        self._widgets['MomAreaEntry'] = Entry(self._widgets['ExplicitPropsFrame'], width=5)
        self._widgets['MomAreaUnit'] = StringVar(self._widgets['ExplicitPropsFrame'])
        self._widgets['MomAreaUnit'].set('--')
        self._widgets['MomAreaMenu'] = OptionMenu(self._widgets['ExplicitPropsFrame'], self._widgets['MomAreaUnit'],
                                                  *['cm^4', 'in^4'])
        self._widgets['MomAreaEntry'].grid(row=1, column=1)
        self._widgets['MomAreaMenu'].grid(row=1, column=2)

        self._widgets['AdvancedPlotControl'] = Frame(self._widgets['AdvancedPlotFrame'])
        self._widgets['AdvancedPlotControl'].pack()
        self._widgets['LockAdvancedProp'] = Button(self._widgets['AdvancedPlotControl'], text="Submit Properties",
                                                   command=lambda: self._get_advanced_properties())
        self._widgets['LockAdvancedProp'].config(font=("Helvetica", 15))
        self._widgets['LockAdvancedProp'].pack(side=LEFT, padx=30)
        self._widgets['AdvancedPropTargets'] = Frame(self._widgets['AdvancedPlotControl'])
        self._widgets['AdvancedPropTargets'].pack(side=RIGHT)
        self._widgets['AdvancedPlotDisplay'] = Frame(self._widgets['MechOutputFrame'])
        self._widgets['AdvancedPlotDisplay'].pack()
        return

    def _pop_adv_plot_control(self):
        self._widgets['AdvancedPlotVar'] = StringVar(self._widgets['AdvancedPropTargets'])
        self._widgets['AdvancedPlotVar'].set("Select Plot")
        self._widgets['AdvancedPlotMenu'] = OptionMenu(self._widgets['AdvancedPropTargets'],
                                                       self._widgets['AdvancedPlotVar'],
                                                       *['Deflection', 'Shear Stress', 'Axial Stress', 'von Mises'])
        self._widgets['AdvancedPlotMenu'].grid(row=0, column=0)
        self._widgets['AdvancedPlotButton'] = Button(self._widgets['AdvancedPropTargets'], text="Plot",
                                                     command=lambda: self._get_advanced_plots())
        self._widgets['AdvancedPlotButton'].grid(row=0, column=1)
        Label(self._widgets['AdvancedPropTargets'],
              text=f'Select Stress Location: {self._len_units}').grid(row=1, column=0)
        _scaled_len = self._beam_len / self._units.get_len_conversion(self._len_units)
        self._widgets['StressLocationSelector'] = Scale(self._widgets['AdvancedPropTargets'], from_=0, to=_scaled_len,
                                                        resolution=_scaled_len / self._resolution, orient=HORIZONTAL,
                                                        command=lambda x: self._update_results(
                                                            self._widgets['StressLocationSelector']))
        self._widgets['StressLocationSelector'].grid(row=2, column=0)
        return

    def _defined_implicit_frame(self):
        self._widgets['ImplicitComFrame'] = LabelFrame(self._widgets['ImplicitPropsFrame'])
        self._widgets['ImplicitSketchFrame'] = LabelFrame(self._widgets['ImplicitPropsFrame'])
        self._widgets['ImplicitComFrame'].pack()
        self._widgets['ImplicitSketchFrame'].pack()
        Label(self._widgets['ImplicitComFrame'], text="Draw Custom Cross Section").pack()

        self._widgets['ImplicitControl'] = Frame(self._widgets['ImplicitSketchFrame'])
        self._widgets['ImplicitControl'].pack()
        self._widgets['ImplicitConfirm'] = Frame(self._widgets['ImplicitSketchFrame'])
        self._widgets['ImplicitConfirm'].pack()
        Label(self._widgets['ImplicitControl'], text='Select Position Unit: ').grid(row=0, column=0)
        self._widgets['CrossSectionVar'] = StringVar(self._widgets['ImplicitControl'])
        self._widgets['CrossSectionVar'].set("--")
        self._widgets['CrossSectionMenu'] = OptionMenu(self._widgets['ImplicitControl'],
                                                       self._widgets['CrossSectionVar'], *self._units.get_len_units())
        self._widgets['CrossSectionMenu'].grid(row=0, column=1)
        Label(self._widgets['ImplicitControl'], text='Enter Rectangle Width: ').grid(row=1, column=0)
        self._widgets['RectangleWidth'] = Entry(self._widgets['ImplicitControl'], width=5)
        self._widgets['RectangleWidth'].grid(row=1, column=1)
        Label(self._widgets['ImplicitControl'], text='Enter Rectangle Height: ').grid(row=2, column=0)
        self._widgets['RectangleHeight'] = Entry(self._widgets['ImplicitControl'], width=5)
        self._widgets['RectangleHeight'].grid(row=2, column=1)
        Label(self._widgets['ImplicitControl'], text='Enter Rectangle X: ').grid(row=3, column=0)
        self._widgets['RectangleX'] = Entry(self._widgets['ImplicitControl'], width=5)
        self._widgets['RectangleX'].grid(row=3, column=1)
        Label(self._widgets['ImplicitControl'], text='Enter Rectangle Y: ').grid(row=4, column=0)
        self._widgets['RectangleY'] = Entry(self._widgets['ImplicitControl'], width=5)
        self._widgets['RectangleY'].grid(row=4, column=1)
        self._widgets['SubmitSupport'] = Button(self._widgets['ImplicitConfirm'], text="Submit",
                                                command=lambda: self._add_rectangle())
        self._widgets['ClearRecentSupport'] = Button(self._widgets['ImplicitConfirm'], text='Undo',
                                                     command=lambda: self._clear_rectangle())
        self._widgets['SubmitSupport'].config(font=("Helvetica", 15))
        self._widgets['ClearRecentSupport'].config(font=("Helvetica", 15))
        self._widgets['SubmitSupport'].grid(row=0, column=0, padx=10, pady=5)
        self._widgets['ClearRecentSupport'].grid(row=0, column=1, padx=10, pady=6)

    def _add_rectangle(self):
        try:
            if not (self._cross_section_unit == self._widgets['CrossSectionVar'].get()):
                _current = self._widgets['CrossSectionVar'].get()
                if _current == '--':
                    raise RuntimeError("Please select a unit")
                if self._cross_section_unit:
                    self.output("Cross-section units have been changed, clearing old design..", adv_window=True)
                    self._rectangles = []
                self._cross_section_unit = _current
            width = float(self._widgets['RectangleWidth'].get())
            height = float(self._widgets['RectangleHeight'].get())
            x = float(self._widgets['RectangleX'].get())
            y = float(self._widgets['RectangleY'].get())
        except ValueError:
            self.output('Invalid location/dimension for cross-section member..', adv_window=True)
            return
        except RuntimeError as e:
            self.output(f"Runtime error in cross-section design: {e}")
            return
        self._rectangles.append({
            'b': width,
            'h': height,
            'x_loc': x,
            'y_loc': y
        })
        self._widgets['RectangleWidth'].delete(0, 'end')
        self._widgets['RectangleHeight'].delete(0, 'end')
        self._widgets['RectangleX'].delete(0, 'end')
        self._widgets['RectangleY'].delete(0, 'end')
        print(self._rectangles)
        self._update_cross_section_plot()
        return

    def _clear_rectangle(self):
        self._rectangles.pop()
        self._update_cross_section_plot()
        print(self._rectangles)

    def _update_cross_section_plot(self):
        try:
            self._widgets['CrossSecPlot'].get_tk_widget().pack_forget()
        except KeyError:
            pass
        self._widgets['CrossSecPlot'] = GUIPlotter.plot_cross_section(self._rectangles,
                                                                      self._widgets['ImplicitCrossSection'],
                                                                      f'{self._cross_section_unit}',
                                                                      f'{self._cross_section_unit}',
                                                                      "Beam Cross Section")
        self._widgets['CrossSecPlot'].get_tk_widget().pack()
        self._widgets['CrossSecPlot'].draw()

    def _get_advanced_plots(self):
        target_ = self._widgets['AdvancedPlotVar'].get()
        if target_ == "Deflection":
            self._get_deflection_plots()
        if ("Stress" in target_) or ("von" in target_):
            self._get_stress_plots(target_)
        else:
            self.output("Please select a plot!", adv_window=True)
        return

    def _get_advanced_properties(self):
        try:
            u_e = self._widgets['YoungsUnit'].get()
            _e = float(self._widgets['YoungsEntry'].get())
            if u_e == '--':
                raise ValueError("Please select a valid E unit")
            if "Pa" in u_e:
                self._stress_units = 'MPa'
            else:
                self._stress_units = 'kpsi'
            self._advanced_properties['E'] = _e*self._units.get_stress_conversion(u_e)
        except ValueError as e:
            self.output(e, adv_window=True)
            return
        if not self._rectangles:
            self.output("User has not developed a custom cross-section", adv_window=True)
            try:
                u_i = self._widgets['MomAreaUnit'].get()
                _i = float(self._widgets['MomAreaEntry'].get())
                if u_i == '--':
                    raise ValueError("Please enter a valid I unit")
                self._advanced_properties['I'] = _i*(self._units.get_len_conversion(u_i.split('^')[0]))**4
            except ValueError as e:
                self.output(e, adv_window=True)
                return
        else:
            _a = 0
            _ya = 0
            _i = 0
            for rect in self._rectangles:
                _a += rect['b'] * rect['h']
                _ya += ((rect['h']/2)+rect['y_loc'])*rect['b']*rect['h']
            _yna = _ya / _a

            for rect in self._rectangles:
                i1 = (1/12)*rect['b']*(rect['h']**3)
                i2 = rect['b']*rect['h']*(_yna-rect['y_loc']-(rect['h']/2))**2
                _i += (i1 + i2)
            self._advanced_properties['I'] = _i*(self._units.get_len_conversion(self._cross_section_unit))**4
            self._advanced_properties['A'] = _a*(self._units.get_len_conversion(self._cross_section_unit))**2
            self._advanced_properties['yNA'] = _yna*(self._units.get_len_conversion(self._cross_section_unit))
        self._pop_adv_plot_control()
        return

    def _get_stress_plots(self, target_plot):
        try:
            self._widgets['AdvPlotView'].get_tk_widget().pack_forget()
        except KeyError:
            pass
        if not self._rectangles:
            self.output("Cannot solve stresses without cross-section", adv_window=True)
            return
        _loc = self._widgets['StressLocationSelector'].get()
        _beam = self._results_dict['Beam']
        loc_diff = [abs(_loc - target_val) for target_val in _beam]
        target_ind = loc_diff.index(min(loc_diff))
        _axial_force = self._results_dict['AxialForce'][target_ind]
        _shear_force = self._results_dict['ShearForce'][target_ind]
        _bending_moment = self._results_dict['BendingMoment'][target_ind]
        stress_sol = ss.StressSolver(self._rectangles, self._advanced_properties, _bending_moment, _shear_force,
                                     _axial_force)
        shear_s, axial_s, sec_size = stress_sol.get_stresses()
        for row, shear_x in enumerate(shear_s):
            for index, shear in enumerate(shear_x):
                shear_s[row][index] = shear/self._units.get_stress_conversion(self._stress_units)
        for row, axial_x in enumerate(axial_s):
            for index, axial in enumerate(axial_x):
                axial_s[row][index] = axial/self._units.get_stress_conversion(self._stress_units)
        if target_plot == "Axial Stress":
            self._widgets['AdvPlotView'] = GUIPlotter.plot_cross_section(
                self._rectangles, self._widgets['AdvancedPlotDisplay'],self._cross_section_unit,
                self._cross_section_unit, f"Axial Stress in {self._stress_units}", axial_s, sec_size)
        elif target_plot == "Shear Stress":
            self._widgets['AdvPlotView'] = GUIPlotter.plot_cross_section(
                self._rectangles, self._widgets['AdvancedPlotDisplay'], self._cross_section_unit,
                self._cross_section_unit, f"Shear Stress in {self._stress_units})", shear_s, sec_size)
        self._widgets['AdvPlotView'].get_tk_widget().pack()
        self._widgets['AdvPlotView'].draw()
        self.output(f'{target_plot} plot displayed', adv_window=True)
        print_list1 = []
        print_list2 = []
        return

    def _get_deflection_plots(self):
        x_plot = [x/self._units.get_len_conversion(self._len_units) for x in self._results_dict['Beam']]
        y_plot = [y*self._advanced_properties['E']*self._advanced_properties['I']/(self._units.get_len_conversion(
            self._len_units)) for y in self._results_dict['Deflection']]
        try:
            self._widgets['AdvPlotView'].get_tk_widget().pack_forget()
        except KeyError:
            pass
        self._widgets['AdvPlotView'] = \
            GUIPlotter.make_plot(x_plot, y_plot, self._widgets['AdvancedPlotDisplay'],
                                 f"Beam Length {self._len_units}", f"Deflection {self._len_units}", "Deflection")
        self._widgets['AdvPlotView'].get_tk_widget().pack()
        self._widgets['AdvPlotView'].draw()
        self.output('Deflection plot displayed', adv_window=True)
        return

    def _populate_window(self):
        """
        root method to create interface geometry and layout
        :return:
        """
        # Define upper level tab design for beam_io interface
        self._widgets['TabControl'] = ttk.Notebook(self._root)
        self._widgets['IntroTab'] = ttk.Frame(self._widgets['TabControl'])
        self._widgets['TabControl'].add(self._widgets['IntroTab'], text='Read Me')
        self._widgets['TabControl'].pack(expand=1, fill='both')
        self._widgets['ProblemTab'] = ttk.Frame(self._widgets['TabControl'])
        self._widgets['TabControl'].add(self._widgets['ProblemTab'], text='Start Solving!')
        self._widgets['AdvancedTab'] = ttk.Frame(self._widgets['TabControl'])
        self._widgets['ResultsTab'] = ttk.Frame(self._widgets['TabControl'])
        self._setup_problem_tab()
        self._setup_intro_tab()
        return

    def _setup_intro_tab(self):
        target_dir = os.path.join(os.getcwd(), 'app_data')
        with open(os.path.join(target_dir, "disclaimer.txt"), 'r') as f:
            dis_text = f.read()
        self._widgets['Disclaimer'] = Message(self._widgets['IntroTab'], width=805)
        self._widgets['Disclaimer'].pack()
        self._widgets['Disclaimer'].config(text=dis_text)
        self._widgets['Disclaimer2'] = ImageTk.PhotoImage(Image.open(os.path.join(target_dir, "disclaimer.png")))
        self._widgets['Disclaimer2Disp'] = Label(self._widgets['IntroTab'], image=self._widgets['Disclaimer2'])
        self._widgets['Disclaimer2Disp'].pack()


    def _display_reactions(self):
        """
        populate the reactions window with the reaction loads calculated for each support in the structure
        :return:
        """
        reactions = self._results_dict['Reactions']
        _row = 0
        _column = 0
        for reaction in reactions:
            reac_str = ""
            for key, line in reaction.items():
                if 'Sup' in key:
                    _type = key.split('Sup')[0]
                    _loc = float(line)/self._units.get_len_conversion(self._len_units)
                    reac_str += f"{_type} support at {_loc}{self._len_units}\n"
                    continue
                if "moment" in key:
                    conv_ = self._units.get_combined_conversion("Moment Load", self._len_units, self._force_units)
                    units_ = self._units.get_combined_unit("Moment Load", self._len_units, self._force_units)
                else:
                    conv_ = self._units.get_force_conversion(self._force_units)
                    units_ = self._force_units
                reac_str += f'{key} reaction: {float(line)/conv_}{units_}\n'
            txt = Text(self._widgets['ReactionsFrame'], height=6, width=30)
            txt.grid(row=_row, column=_column)
            _column += 1
            if _column > 3:
                _row += 1
                _column = 0
            txt.insert(INSERT, reac_str)
        return

    def _setup_results_tab(self):
        """
        populate the force and loading results tab with tkinter widgets to display the reaction loads
        :return:
        """
        self._widgets['TabControl'].add(self._widgets['ResultsTab'], text='Results')
        self._widgets['ReactionsFrame'] = LabelFrame(self._widgets['ResultsTab'])
        self._widgets['ReactionsFrame'].pack()
        self._widgets['SliderFrame'] = LabelFrame(self._widgets['ResultsTab'])
        self._widgets['SliderFrame'].pack()
        self._widgets['PlotsFrame'] = LabelFrame(self._widgets['ResultsTab'])
        self._widgets['PlotsFrame'].pack()
        self._display_reactions()
        Label(self._widgets['SliderFrame'], text='Please Select Desired Location').grid(row=0, column=0)
        _scaled_len = self._beam_len/self._units.get_len_conversion(self._len_units)
        self._widgets['LocationSelector'] = Scale(self._widgets['SliderFrame'], from_=0, to=_scaled_len,
                                                  resolution=_scaled_len/self._resolution, orient=HORIZONTAL,
                                                  command=lambda x: self._update_results(
                                                      self._widgets['LocationSelector']))
        self._widgets['LocationSelector'].grid(row=2, column=0)
        Label(self._widgets['SliderFrame'], text="Results at Selected Location").grid(row=0, column=1)
        self._widgets['AxialResults'] = Label(self._widgets['SliderFrame'], width=50)
        self._widgets['AxialResults'].grid(row=1, column=1)
        self._widgets['ShearResults'] = Label(self._widgets['SliderFrame'], width=50)
        self._widgets['ShearResults'].grid(row=2, column=1)
        self._widgets['MomentResults'] = Label(self._widgets['SliderFrame'], width=50)
        self._widgets['MomentResults'].grid(row=3, column=1)
        self._create_plotter()
        return

    def _create_plotter(self):
        """
        populate plotter frame with tkinter widgets
        :return:
        """
        self._widgets['PlotSelectionFrame'] = Frame(self._widgets['PlotsFrame'])
        self._widgets['PlotDisplay'] = Frame(self._widgets['PlotsFrame'])
        self._widgets['PlotSelectionFrame'].pack()
        self._widgets['PlotDisplay'].pack()
        self._widgets['PlotSelectionVar'] = StringVar(self._widgets['PlotSelectionFrame'])
        self._widgets['PlotSelectionVar'].set("Select Plot")
        self._widgets['PlotSelectionMenu'] = OptionMenu(self._widgets['PlotSelectionFrame'],
                                                        self._widgets['PlotSelectionVar'],
                                                        *["Bending Moment", "Shear Force", "Axial Force"])
        self._widgets['PlotSelectionMenu'].pack(side=LEFT)
        self._widgets['PlotterButton'] = Button(self._widgets['PlotSelectionFrame'], text="Plot",
                                                command=lambda: self._display_image())
        self._widgets['PlotterButton'].pack(side=LEFT)
        return

    def _display_image(self):
        """
        create a tkinter Figure of the desired plot on the basis of user input
        :return:
        """
        target_file = self._widgets["PlotSelectionVar"].get()
        if "Select Plot" == target_file:
            return
        target_ = target_file.replace(' ', '')
        x_to_plot = [x/self._units.get_len_conversion(self._len_units) for x in self._results_dict['Beam']]
        if "Moment" in target_:
            units_ = self._units.get_combined_unit('Moment Load', self._len_units, self._force_units)
            conv_ = self._units.get_combined_conversion('Moment Load', self._len_units, self._force_units)
        else:
            units_ = self._force_units
            conv_ = self._units.get_force_conversion(self._force_units)
        y_to_plot = [y/conv_ for y in self._results_dict[target_]]
        try:
            self._widgets['PlotView'].get_tk_widget().pack_forget()
        except KeyError:
            pass
        self._widgets['PlotView'] = \
            GUIPlotter.make_plot(x_to_plot, y_to_plot, self._widgets['PlotDisplay'],
                                 f"Beam Length {self._len_units}", f"{target_file} {units_}", target_file)
        self._widgets['PlotView'].get_tk_widget().pack()
        self._widgets['PlotView'].draw()
        return

    def _update_results(self, loc_slider):
        """
        provide the user with loading - shear, bending moment and axial loading values for a given value on the length
        of theh structure on the basis of a Tkinter slider input
        :param loc_slider: tkinter slider widget
        :return:
        """
        target_location = loc_slider.get()*self._units.get_len_conversion(self._len_units)
        _ax = self._results_dict['AxialForce']
        _sf = self._results_dict['ShearForce']
        _bm = self._results_dict['BendingMoment']
        _beam = self._results_dict['Beam']
        loc_diff = [abs(target_location - target_val) for target_val in _beam]
        target_ind = loc_diff.index(min(loc_diff))
        axial_load = _ax[target_ind]/self._units.get_force_conversion(self._force_units)
        vertical_load = _sf[target_ind]/self._units.get_force_conversion(self._force_units)
        bending_moment = _bm[target_ind]/self._units.get_combined_conversion('Moment Load', self._len_units,
                                                                             self._force_units)
        axial_text = f"Axial Loading: {axial_load} {self._force_units}"
        vertical_text = f"Vertical Loading: {vertical_load} {self._force_units}"
        moment_text = f"Moment Loading: {bending_moment} {self._force_units}{self._len_units}"
        self._widgets['AxialResults'].config(text=axial_text)
        self._widgets['ShearResults'].config(text=vertical_text)
        self._widgets['MomentResults'].config(text=moment_text)
        return

    def _setup_problem_tab(self):
        """
        develop the window to be populated with widgets to accept problem specifics
        :return:
        """
        # Populate sub-frames in the problem setup tab
        self._widgets['BeamFrame'] = LabelFrame(self._widgets['ProblemTab'], padx=30)
        self._widgets['BeamFrame'].grid(row=0, column=0)
        self._widgets['ControlFrame'] = LabelFrame(self._widgets['ProblemTab'])
        self._widgets['ControlFrame'].grid(row=0, column=1)
        self._widgets['CommandWindow'] = LabelFrame(self._widgets['ProblemTab'])
        self._widgets['CommandWindow'].grid(row=0, column=2)
        self._widgets['SupportFrame1'] = LabelFrame(self._widgets['ProblemTab'])
        self._widgets['SupportFrame1'].grid(row=1, column=0)
        self._widgets['SupportFrame2'] = LabelFrame(self._widgets['ProblemTab'])
        self._widgets['SupportFrame2'].grid(row=1, column=1)
        self._widgets['SupportFrame3'] = LabelFrame(self._widgets['ProblemTab'])
        self._widgets['SupportFrame3'].grid(row=1, column=2)
        Frame(self._widgets['ProblemTab']).grid(row=2, column=0)
        self._widgets['LoadFrame1'] = LabelFrame(self._widgets['ProblemTab'])
        self._widgets['LoadFrame1'].grid(row=3, column=0)
        self._widgets['LoadFrame2'] = LabelFrame(self._widgets['ProblemTab'])
        self._widgets['LoadFrame2'].grid(row=3, column=1)
        self._widgets['LoadFrame3'] = LabelFrame(self._widgets['ProblemTab'])
        self._widgets['LoadFrame3'].grid(row=3, column=2)

        # add widgets to the problem setup tab
        self._widgets['BeamUnits'] = StringVar(self._widgets['BeamFrame'])
        self._widgets['BeamUnits'].set('--')
        self._widgets['ForceUnits'] = StringVar(self._widgets['BeamFrame'])
        self._widgets['ForceUnits'].set('--')
        self._widgets['BeamLength'] = Entry(self._widgets['BeamFrame'], width=5)
        self._widgets['BeamUnitSelection'] = OptionMenu(self._widgets['BeamFrame'], self._widgets['BeamUnits'],
                                                        *self._units.get_len_units())
        self._widgets['ForceUnitSelection'] = OptionMenu(self._widgets['BeamFrame'], self._widgets['ForceUnits'],
                                                         *self._units.get_force_units())
        self._widgets['SubmitBeamButton'] = Button(self._widgets['BeamFrame'], text='Enter Beam',
                                                   command=lambda: self._define_beam())
        Label(self._widgets['BeamFrame'], text="Beam Length:").grid(row=0, column=0)
        self._widgets['BeamLength'].grid(row=0, column=1)
        Label(self._widgets['BeamFrame'], text='Select Length Units:').grid(row=1, column=0)
        self._widgets['BeamUnitSelection'].grid(row=1, column=1)
        Label(self._widgets['BeamFrame'], text='Select Force Units:').grid(row=2, column=0)
        self._widgets['ForceUnitSelection'].grid(row=2, column=1)
        self._widgets['SubmitBeamButton'].grid(row=3, column=0)

        self._widgets['CmdScroll'] = Scrollbar(self._widgets['CommandWindow'])
        self._widgets['CommandLine'] = Text(self._widgets['CommandWindow'], height=7, width=50)
        self._widgets['CmdScroll'].pack(side=RIGHT)
        self._widgets['CommandLine'].pack(side=LEFT)
        self._widgets['CmdScroll'].config(command=self._widgets['CommandLine'].yview)
        self._widgets['CommandLine'].config(yscrollcommand=self._widgets['CmdScroll'].set)

        self._widgets['SolvePerc'] = Label(self._widgets['ControlFrame'], text='SOLVED:\n 0%')
        self._widgets['SolvePerc'].config(font=("Arial", 30))
        self._widgets['SolvePerc'].pack(pady=9)

        self._widgets['SupportType'] = StringVar(self._widgets['SupportFrame1'])
        self._widgets['SupportType'].set("Select Support")
        self._widgets['SupportMenu'] = OptionMenu(self._widgets['SupportFrame1'], self._widgets['SupportType'],
                                                  *self._setup.get_supports())
        self._widgets['SupportMenu'].grid(row=0, column=0, padx=52)
        Label(self._widgets['SupportFrame1'], text='\t').grid(row=1, column=0)
        Label(self._widgets['SupportFrame1'], text='\t').grid(row=2, column=0)
        Label(self._widgets['SupportFrame1'], text='\t').grid(row=3, column=0)
        self._widgets['Sub1SupportFrame2'] = Frame(self._widgets['SupportFrame2'])
        self._widgets['Sub2SupportFrame2'] = Frame(self._widgets['SupportFrame2'])
        self._widgets['Sub1SupportFrame2'].pack()
        self._widgets['Sub2SupportFrame2'].pack()
        Label(self._widgets['Sub1SupportFrame2'], text='Enter Support').grid(row=0, column=0)
        Label(self._widgets['Sub1SupportFrame2'], text='Location').grid(row=0, column=1)
        self._widgets['SupportLocation'] = Entry(self._widgets['Sub1SupportFrame2'], width=10)
        self._widgets['SupportLocation'].grid(row=1, column=0)
        self._widgets['SupportLocUnit'] = Label(self._widgets['Sub1SupportFrame2'], text='')
        self._widgets['SupportLocUnit'].grid(row=1, column=1)

        self._widgets['SubmitSupport'] = Button(self._widgets['Sub2SupportFrame2'], text="Submit",
                                                command=lambda: self._determine_support())
        self._widgets['ClearRecentSupport'] = Button(self._widgets['Sub2SupportFrame2'], text='Undo',
                                                     command=lambda: self._clear_last_support())
        self._widgets['SubmitSupport'].config(font=("Helvetica", 15))
        self._widgets['ClearRecentSupport'].config(font=("Helvetica", 15))
        self._widgets['SubmitSupport'].grid(row=0, column=0, padx=10, pady=5)
        self._widgets['ClearRecentSupport'].grid(row=0, column=1, padx=10, pady=6)
        self._widgets['SupportScroll'] = Scrollbar(self._widgets['SupportFrame3'])
        self._widgets['SupportLine'] = Text(self._widgets['SupportFrame3'], height=7, width=50)
        self._widgets['SupportScroll'].pack(side=RIGHT)
        self._widgets['SupportLine'].pack(side=LEFT)
        self._widgets['SupportScroll'].config(command=self._widgets['SupportLine'].yview)
        self._widgets['SupportLine'].config(yscrollcommand=self._widgets['SupportScroll'].set)

        self._widgets['Sub1LoadFrame2'] = Frame(self._widgets['LoadFrame2'])
        self._widgets['Sub2LoadFrame2'] = Frame(self._widgets['LoadFrame2'])
        self._widgets['Sub1LoadFrame2'].pack()
        self._widgets['Sub2LoadFrame2'].pack()
        self._widgets['LoadType'] = StringVar(self._widgets['LoadFrame1'])
        self._widgets['LoadType'].set("Select Load")
        self._widgets['LoadType'].trace('w', self._populate_load)
        self._widgets['LoadMenu'] = OptionMenu(self._widgets['LoadFrame1'], self._widgets['LoadType'],
                                               *self._setup.get_loads())
        self._widgets['LoadMenu'].grid(row=0, column=0, padx=60)
        Label(self._widgets['LoadFrame1'], text='\t').grid(row=1, column=0, pady=63)
        self._widgets['LoadRequirements'] = Label(self._widgets['LoadFrame1'], text='\t').grid(row=2, column=0)
        Label(self._widgets['Sub1LoadFrame2'], text='Enter Load').grid(row=0, column=0)
        Label(self._widgets['Sub1LoadFrame2'], text='Config:').grid(row=0, column=1)
        Label(self._widgets['Sub1LoadFrame2'], text='  ').grid(row=0, column=2)
        self._widgets['LoadStartLabel'] = Label(self._widgets['Sub1LoadFrame2'], text='\t').grid(row=1, column=0)
        self._widgets['LoadValue1'] = Label(self._widgets['Sub1LoadFrame2'], text='\t')
        self._widgets['LoadValue1'].grid(row=2, column=0)
        self._widgets['LoadValueEntry1'] = Entry(self._widgets['Sub1LoadFrame2'], width=7)
        self._widgets['LoadValueUnits1'] = Label(self._widgets['Sub1LoadFrame2'], text='')
        self._widgets['LoadValueUnits1'].grid(row=2, column=2)
        self._widgets['LoadLoc1'] = Label(self._widgets['Sub1LoadFrame2'], text='\t')
        self._widgets['LoadLoc1'].grid(row=3, column=0)
        self._widgets['LoadLocEntry1'] = Entry(self._widgets['Sub1LoadFrame2'], width=7)
        self._widgets['LoadLocUnits1'] = Label(self._widgets['Sub1LoadFrame2'], text='')
        self._widgets['LoadLocUnits1'].grid(row=3, column=2)
        self._widgets['LoadEndLabel'] = Label(self._widgets['Sub1LoadFrame2'], text='\t')
        self._widgets['LoadEndLabel'].grid(row=4, column=0)
        self._widgets['LoadValue2'] = Label(self._widgets['Sub1LoadFrame2'], text='\t')
        self._widgets['LoadValue2'].grid(row=5, column=0)
        self._widgets['LoadValueEntry2'] = Entry(self._widgets['Sub1LoadFrame2'], width=7)
        self._widgets['LoadValueUnits2'] = Label(self._widgets['Sub1LoadFrame2'], text='')
        self._widgets['LoadValueUnits2'].grid(row=5, column=2)
        self._widgets['LoadLoc2'] = Label(self._widgets['Sub1LoadFrame2'], text='\t')
        self._widgets['LoadLoc2'].grid(row=6, column=0)
        self._widgets['LoadLocEntry2'] = Entry(self._widgets['Sub1LoadFrame2'], width=7)
        self._widgets['LoadLocUnits2'] = Label(self._widgets['Sub1LoadFrame2'], text='')
        self._widgets['LoadLocUnits2'].grid(row=6, column=2)

        self._widgets['SubmitLoad'] = Button(self._widgets['Sub2LoadFrame2'], text="Submit",
                                             command=lambda: self._determine_load())
        self._widgets['ClearRecentLoad'] = Button(self._widgets['Sub2LoadFrame2'], text='Undo',
                                                  command=lambda: self._clear_last_load())
        self._widgets['SubmitLoad'].config(font=("Helvetica", 15))
        self._widgets['ClearRecentLoad'].config(font=("Helvetica", 15))
        self._widgets['SubmitLoad'].grid(row=0, column=0, padx=10, pady=5)
        self._widgets['ClearRecentLoad'].grid(row=0, column=1, padx=10, pady=6)

        self._widgets['LoadScroll'] = Scrollbar(self._widgets['LoadFrame3'])
        self._widgets['LoadLine'] = Text(self._widgets['LoadFrame3'], height=14, width=50)
        self._widgets['LoadScroll'].pack(side=RIGHT)
        self._widgets['LoadLine'].pack(side=LEFT)
        self._widgets['LoadScroll'].config(command=self._widgets['LoadLine'].yview)
        self._widgets['LoadLine'].config(yscrollcommand=self._widgets['LoadScroll'].set)
        self._widgets['SolveButton1'] = Button(self._widgets['ProblemTab'], text="SOLVE", command=lambda:
                                               self._solve_beam(), width=18)
        self._widgets['SolveButton1'].grid(row=4, column=2, pady=15)
        self._widgets['SolveButton1'].config(font=("Helvetica", 30))
        return

    def _define_beam(self, *args):
        """
        create the structural definition of the beam being solved - length and characteristic units desired by the user
        :param args:
        :return:
        """
        self._loads = []
        self._supports = []
        self._len_units = self._widgets['BeamUnits'].get()
        self._force_units = self._widgets['ForceUnits'].get()
        try:
            if not (self._len_units in self._units.get_len_units()):
                raise RuntimeError("Please select a length unit")
            if not (self._force_units in self._units.get_force_units()):
                raise RuntimeError("Please select a force unit")
            self._beam_len = float(self._widgets['BeamLength'].get())*self._units.get_len_conversion(self._len_units)
        except ValueError as e:
            self.output(e)
            self.output("Invalid Beam Length")
            return
        except RuntimeError as e:
            self.output(e)
            return
        self._widgets['ForceUnits'].trace('w', self._define_beam)
        self._widgets['BeamUnits'].trace('w', self._define_beam)
        self._widgets['SupportLocUnit'].config(text=self._len_units)
        self.output("Beam has been defined!")
        self._widgets['SupportType'].set("Select Support")
        self._widgets['LoadType'].set("Select Load")
        self._beam_defined = True
        self._update_display()
        return

    def _clear_last_support(self):
        """
        delete the most recent support from the support list
        :return: whether the deletion was successful
        """
        try:
            _cleared = self._supports.pop()
        except IndexError:
            self.output("No supports exist!")
            return False
        self.output(f"Cleared:\nType: {_cleared['support']}\nLocation: {_cleared['location']}{self._len_units}")
        self._update_display()
        return True

    def _determine_support(self):
        """
        parse entries to accept or decline a new support
        :return: whether the support was accepted
        """
        try:
            if self._beam_len is None:
                raise RuntimeError("Define beam first!")
            _support_type = self._widgets['SupportType'].get()
            if _support_type == "Select Support":
                raise RuntimeError("Please select a support type")
            _support_location = float(self._widgets['SupportLocation'].get())*self._units.get_len_conversion(
                self._len_units)
            if (_support_location > self._beam_len) or _support_location < 0:
                raise RuntimeError("Please place support on the beam")
        except ValueError as e:
            self.output(e)
            self.output("Please enter a valid support location")
            return False
        except RuntimeError as e:
            self.output(e)
            return False
        self.output("Support has been accepted")
        self.output(f"Type: {_support_type}\n"
                    f"Location: {_support_location/self._units.get_len_conversion(self._len_units)}{self._len_units}")
        self._supports.append({"support": _support_type, "location": _support_location})
        self._widgets['SupportLocation'].delete(0, 'end')
        self._widgets['SupportType'].set("Select Support")
        self._update_display()
        return True

    def _populate_load(self, *args):
        """
        populate load window to accept information about the selected load
        :param args:
        :return:
        """
        try:
            if self._beam_len is None:
                raise RuntimeError("Define beam first!")
            _load_type = self._widgets['LoadType'].get()
            if _load_type == "Select Load":
                if not self._beam_defined:
                    return False
                raise RuntimeError("Please select a load type")
            combined_unit = self._units.get_combined_unit(_load_type, self._len_units, self._force_units)
            def_conditions = self._setup.get_definition_conditions(_load_type)
            if not combined_unit:
                raise RuntimeError("Something went wrong while selecting load")
        except RuntimeError as e:
            self.output(e)
            return
        ident_1 = "Enter"
        ident_2 = None
        if def_conditions == 2:
            ident_1 = "Start"
            ident_2 = "End"
        self._widgets['LoadValue1'].config(text=f'{ident_1} val:')
        self._widgets['LoadLoc1'].config(text=f'{ident_1} loc:')
        self._widgets['LoadValueEntry1'].grid(row=2, column=1)
        self._widgets['LoadLocEntry1'].grid(row=3, column=1)
        self._widgets['LoadLocUnits1'].config(text=f'{self._len_units}')
        self._widgets['LoadValueUnits1'].config(text=f'{combined_unit}')
        if ident_2:
            self._widgets['LoadValue2'].config(text=f'{ident_2} val:')
            self._widgets['LoadLoc2'].config(text=f'{ident_2} loc:')
            self._widgets['LoadValueEntry2'].grid(row=5, column=1)
            self._widgets['LoadLocEntry2'].grid(row=6, column=1)
            self._widgets['LoadLocUnits2'].config(text=f'{self._len_units}')
            self._widgets['LoadValueUnits2'].config(text=f'{combined_unit}')
        else:
            self._widgets['LoadValue2'].config(text='')
            self._widgets['LoadLoc2'].config(text='')
            self._widgets['LoadValueEntry2'].grid_remove()
            self._widgets['LoadLocEntry2'].grid_remove()
            self._widgets['LoadLocUnits2'].config(text='')
            self._widgets['LoadValueUnits2'].config(text='')
        return

    def _determine_load(self):
        """
        obtain and validate information about the definition of the input loading
        :return:
        """
        add_load = {}
        try:
            if self._beam_len is None:
                raise RuntimeError("Define beam first!")
            _load_type = self._widgets['LoadType'].get()
            if _load_type == "Select Support":
                raise RuntimeError("Please select a load type")
            _load_loc = float(self._widgets['LoadLocEntry1'].get())*self._units.get_len_conversion(
                self._len_units)
            if (_load_loc > self._beam_len) or _load_loc < 0:
                raise RuntimeError("Please place support on the beam")
            _load_value = float(self._widgets['LoadValueEntry1'].get())*self._units.get_combined_conversion(
                _load_type, self._len_units, self._force_units)
            def_conditions = self._setup.get_definition_conditions(_load_type)
            if def_conditions == 2:
                _load_loc2 = float(self._widgets['LoadLocEntry2'].get())*self._units.get_len_conversion(
                    self._len_units)
                if (_load_loc2 > self._beam_len) or _load_loc2 < 0:
                    raise RuntimeError("Please place support on the beam")
                _load_value2 = float(self._widgets['LoadValueEntry2'].get())*self._units.get_combined_conversion(
                    _load_type, self._len_units, self._force_units)
                add_load['val2'] = _load_value2
                add_load['loc2'] = _load_loc2
                self._widgets['LoadValueEntry2'].delete(0, 'end')
                self._widgets['LoadLocEntry2'].delete(0, 'end')
            self._widgets['LoadLocEntry1'].delete(0, 'end')
            self._widgets['LoadValueEntry1'].delete(0, 'end')
            add_load['type'] = _load_type
            add_load['val'] = _load_value
            add_load['loc'] = _load_loc
        except ValueError as e:
            self.output(e)
            self.output("Please enter a valid support location/value")
            return False
        except RuntimeError as e:
            self.output(e)
            return False
        self._loads.append(add_load)
        self._update_display()
        return

    def _clear_last_load(self):
        """
        undo the last load entry, update the most recent entry to the second last entry on the time of the function call
        :return:
        """
        try:
            _cleared = self._loads.pop()
        except IndexError:
            self.output("No loads exist!")
            return False
        out_state = f"Cleared:\n{_cleared['type']}\n"
        out_loc = f"Location {float(_cleared['loc'])/self._units.get_len_conversion(self._len_units)}"
        conversion = self._units.get_combined_conversion(_cleared['type'], self._len_units, self._force_units)
        out_load = \
            f"Value {float(_cleared['val'])/conversion}"
        try:
            out_loc += f" to {float(_cleared['loc2']/self._units.get_len_conversion(self._len_units))}"
            out_load = f" to {float(_cleared['val2'])/conversion}"
        except KeyError:
            pass

        out_loc += f"{self._len_units}\n"
        out_load += f"{self._units.get_combined_unit(_cleared['type'], self._len_units, self._force_units)}\n"

        self.output(out_state + out_loc + out_load)
        self._update_display()
        return True

    def _update_display(self):
        """
        update loading and support information displays when a new beam is defined or a new feature is added
        :return:
        """
        self._widgets['SupportLine'].delete("1.0", END)
        self._widgets['LoadLine'].delete("1.0", END)
        add_text = ""
        for support in self._supports:
            add_text += f'{support["support"]} support at ' \
                        f'{float(support["location"]/self._units.get_len_conversion(self._len_units))}{self._len_units}'
            add_text += '\n'
        self._widgets['SupportLine'].insert(INSERT, add_text)
        add_load_text = ""
        for load in self._loads:
            conversion = self._units.get_combined_conversion(load['type'], self._len_units, self._force_units)
            load_loc = f'{load["type"]}\n{float(load["loc"])/self._units.get_len_conversion(self._len_units)}'
            load_val = f'{float(load["val"])/conversion}'
            try:
                load_loc += f'to {float(load["loc2"])/self._units.get_len_conversion(self._len_units)}'
                load_val += f'to {float(load["val2"])/conversion}'
            except KeyError:
                pass
            load_loc += f"{self._len_units}\n"
            load_val += f"{self._units.get_combined_unit(load['type'], self._len_units, self._force_units)}\n"
            add_load_text += load_loc + load_val
        self._widgets['LoadLine'].insert(INSERT, add_load_text)
        return

    def output(self, string, adv_window=False):
        """
        print string and display it on the commandline provided on the master tkinter window
        :param string: string to be output
        :param adv_window: indicate whether the output goes on the regular or advanced frame command output
        :return:
        """
        print(str(string))
        if not adv_window:
            self._widgets['CommandLine'].insert(INSERT, f'\n{string}')
        else:
            self._widgets['AdvancedStatus'].insert(INSERT, f'\n{string}')
        return


if __name__ == '__main__':
    it = Interface()
