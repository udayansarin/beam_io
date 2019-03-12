"""
author: @Udayan Sarin
Develops an interface for beam_io. The given script accepts user inputs for beam geometry, material and dimensions.
Contains a tab to present plots, stresses and load result information
"""
import os
from tkinter import *
from tkinter import ttk
from PIL import ImageTk, Image
from colorama import Fore

import units
import loads
import beam_solver as bs


class Interface:
    """
    Acts as the central class for controlling the interface for beam_io.
    """

    def __init__(self):
        self._supports = []
        self._loads = []
        self._special_properties = {}
        self._len_units = None
        self._force_units = None
        self._beam_len = None
        self._logo = os.path.join(os.getcwd(), 'I-beam.png.ico')
        self._units = units.UnitControl()
        self._setup = loads.SetupControl()
        self._widgets = {}
        self._root = Tk()
        self._root.title("beam.io")
        self._root.iconbitmap(self._logo)
        self._root.geometry('900x700')
        self._root.resizable(height=False, width=False)
        self._populate_window()
        self._root.mainloop()

    def _solve_beam(self):
        print(self._beam_len)
        print(self._supports)
        print(self._loads)
        bs.SolveProblem(self._beam_len, self._supports, self._loads, 1000)

    def _populate_window(self):
        """
        root method to create interface geometry and layout
        :return:
        """
        # Define upper level tab design for beam_io interface
        self._widgets['TabControl'] = ttk.Notebook(self._root)
        self._widgets['ProblemTab'] = ttk.Frame(self._widgets['TabControl'])
        self._widgets['TabControl'].add(self._widgets['ProblemTab'], text='Basic Geometry')
        self._widgets['AdvancedTab'] = ttk.Frame(self._widgets['TabControl'])
        self._widgets['TabControl'].add(self._widgets['AdvancedTab'], text='Advanced Properties')
        self._widgets['ResultsTab'] = ttk.Frame(self._widgets['TabControl'])
        self._widgets['TabControl'].add(self._widgets['ResultsTab'], text='Results')
        self._widgets['HelpTab'] = ttk.Frame(self._widgets['TabControl'])
        self._widgets['TabControl'].add(self._widgets['HelpTab'], text='Help')
        self._widgets['TabControl'].pack(expand=1, fill='both')

        self._setup_problem_tab()
        self._setup_advanced_tab()
        return

    def _setup_advanced_tab(self):
        self._widgets['StructuralPropertiesFrame'] = LabelFrame(self._widgets['AdvancedTab'])
        self._widgets['CrossSectionFrame'] = LabelFrame(self._widgets['AdvancedTab'])
        self._widgets['StructuralPropertiesFrame'].pack()
        self._widgets['CrossSectionFrame'].pack()
        self._widgets['MessageLabel'] = Label(self._widgets['StructuralPropertiesFrame'],
                                              text='Select Advanced Properties')
        self._widgets['MessageLabel'].grid(row=0, column=0)
        Label(self._widgets['StructuralPropertiesFrame'], text='Enter Elas. Modulus:').grid(row=1, column=0)

    def _setup_problem_tab(self):
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

    def _define_beam(self, *args):
        self._loads = []
        self._supports = []
        self._len_units = self._widgets['BeamUnits'].get()
        self._force_units = self._widgets['ForceUnits'].get()
        try:
            self._beam_len = float(self._widgets['BeamLength'].get())*self._units.get_len_conversion(self._len_units)
            if not (self._len_units in self._units.get_len_units()):
                raise RuntimeError("Please select a length unit")
            if not (self._force_units in self._units.get_force_units()):
                raise RuntimeError("Please select a force unit")
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
        self._update_display()

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
            if _load_type == "Select Support":
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

    def output(self, string):
        """
        print string and display it on the commandline provided on the master tkinter window
        :param string:
        :return:
        """
        print(str(string))
        self._widgets['CommandLine'].insert(INSERT, f'\n{string}')
        return


it = Interface()
