import  os, glob, errno
import pandas as pd
import numpy as np
import tkinter as tk
from scipy import stats
import plotly.express as px
from tkinter import ttk
from tkinter.messagebox import Message
from tkinter.filedialog import askopenfile, askdirectory

class Scoring_Function:
    '''
    tab 1 : Input
    '''
    t1_r1, t1_r2, t1_r3, t1_r4, t1_r5, t1_r6, t1_r7 = 0,1,2,3,4,5,6
    t1_c1,t1_c2, t1_c3, t1_c4 = 1,2,3,4

    def __init__(self, master):
        self.master = master
        self.tabControl = ttk.Notebook(master)
        self.tab1 = ttk.Frame(self.tabControl)
        self.tab2 = ttk.Frame(self.tabControl)
        self.tab3 = ttk.Frame(self.tabControl)
        self.tab4 = ttk.Frame(self.tabControl)
        # self.tab5 = ttk.Frame(self.tabControl)
        self.tabControl.add(self.tab1, text='Load')
        self.tabControl.add(self.tab2, text='RMSE', state=tk.DISABLED)
        self.tabControl.add(self.tab3, text='SF all mol', state=tk.DISABLED) # look at function called each mol
        self.tabControl.add(self.tab4, text='SF each mol', state=tk.DISABLED) # look at function called z
        # self.tabControl.add(self.tab5, text='Other', state=tk.DISABLED)
        self.tabControl.pack(expand=1, fill='both')

        self.solution_calc = None #pd.read_csv(r'/Users/mkr97/Furosemide/CCDC_COOH_useThisone/CLT(N=1000).csv', index_col=0)
        # solid_calc_tmp = None #pd.read_csv(r'/Users/mkr97/Furosemide/Solid_calc.csv', index_col=0)
        # df_calc_calc_tmp = Scoring_Function.calc_shifts(self.solution_calc,solid_calc_tmp)
        # df_calc_calc_tmp = Scoring_Function.average_H(df_calc_calc_tmp, ['H9', 'H10'], 'H9*')
        self.solid_calc =  None #pd.read_csv(r'/Users/mkr97/Furosemide/Solid_calc.csv', index_col=0)
        self.solution_exp =  None #pd.read_csv(r'/Users/mkr97/Furosemide/solution_calc_shift.csv', index_col=0)
        self.solid_exp =  None #pd.read_csv(r'/Users/mkr97/Furosemide/solid_exp.csv', index_col=0)
        # self.df_exp_exp = Scoring_Function.exp_shifts(self.solution_exp, self.solid_exp)
        self.path =  None #r'/Users/mkr97/Furosemide'
        self.in_file = ttk.Label(self.tab1)
        self.in_file.grid(row=Scoring_Function.t1_r6, columnspan=20, sticky=tk.S)
        self.view_input = ttk.Label(self.tab1)
        self.view_input.grid(row=Scoring_Function.t1_r7, columnspan=20, sticky=tk.S)
        self.label_input = ttk.Label(self.tab1)
        self.label_input.grid(row=Scoring_Function.t1_r1, column = Scoring_Function.t1_c4, sticky=tk.S)

        self.Tab1()
        self.Tab2()
        self.Tab3()
        self.Tab4()
        # self.df_calc_calc = Scoring_Function.calc_shifts(self.solution_calc, self.solid_calc)
        # self.df_exp_exp = Scoring_Function.exp_shifts(self.solution_exp, self.solid_exp)


        new_df_l = [] #new_df.index.to_list()  # Exp labels as list
        atom_label_tab1 = tk.StringVar()
        self.new_label_dropdown = ttk.OptionMenu(self.tab1, atom_label_tab1, 'Atom', *new_df_l,
                                            command=self.callback_new_label)  # Dropdown list, needs a title or would use the first input as a title
        self.new_label_dropdown.configure(state=tk.DISABLED)
        self.new_label_dropdown.grid(row=Scoring_Function.t1_r4, column=Scoring_Function.t1_c3, sticky=tk.W, padx=20)
        self.avgButton = ttk.Button(self.tab1, text='Average', command=lambda: self.get_average(), state=tk.DISABLED)
        self.avgButton.grid(row=Scoring_Function.t1_r5,column=Scoring_Function.t1_c3,sticky=tk.W, padx=16)

        self.update_buttons()

    def update_buttons(self):
        if self.solid_calc is not None and self.solution_calc is not None and self.solid_exp is not None and self.solution_exp is not None:
            new_df = pd.concat([self.solid_exp, self.solution_exp], axis=1)  # Exp labels
            self.new_label_dropdown.configure(state=tk.NORMAL)
            atom_list = new_df.index.to_list()
            atom_list.insert(0, 'Atom')
            self.new_label_dropdown.set_menu(*atom_list)
            self.df_calc_calc = Scoring_Function.calc_shifts(self.solution_calc, self.solid_calc)
            self.df_exp_exp = Scoring_Function.exp_shifts(self.solution_exp, self.solid_exp)
            Scoring_Function.switch(self.avgButton)
            self.tabControl.tab(2, state=tk.NORMAL)  # tab3
            self.tabControl.tab(3, state=tk.NORMAL)  # tab4
        if self.solid_calc is not None and self.solution_calc is not None:
            self.tabControl.tab(1, state=tk.NORMAL)  # tab2

    def Tab1(self):
        '''
        Tab 1: Inputs
        '''
        ttk.Label(self.tab1, background='red', width = 2, anchor = tk.W).grid(row=Scoring_Function.t1_r1, column=Scoring_Function.t1_c2, padx = 10)
        ttk.Label(self.tab1, background='red', width = 2, anchor = tk.W).grid(row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c2, padx = 10)
        ttk.Label(self.tab1, background='red', width = 2, anchor = tk.W).grid(row=Scoring_Function.t1_r3, column=Scoring_Function.t1_c2, padx = 10)
        ttk.Label(self.tab1, background='red', width = 2, anchor = tk.W).grid(row=Scoring_Function.t1_r4, column=Scoring_Function.t1_c2, padx = 10)
        ttk.Label(self.tab1, background='red', width = 2, anchor = tk.W).grid(row=Scoring_Function.t1_r5, column=Scoring_Function.t1_c2, padx = 10)

        ttk.Label(self.tab1, text='Calc Solution Shift').grid(row=Scoring_Function.t1_r1, sticky=tk.W)
        ttk.Label(self.tab1, text='Calc Solid Shift').grid(row=Scoring_Function.t1_r2, sticky=tk.W)
        ttk.Label(self.tab1, text='Exp Solution Shift').grid(row=Scoring_Function.t1_r3, sticky=tk.W)
        ttk.Label(self.tab1, text='Exp Solid Shift').grid(row=Scoring_Function.t1_r4, sticky=tk.W)
        ttk.Label(self.tab1, text='Output').grid(row=Scoring_Function.t1_r5, sticky=tk.W)
        ttk.Button(self.tab1, text='Open', command=lambda:self.open_sol_calc_file('Open Calculated Solution files')).grid(row=Scoring_Function.t1_r1, column=Scoring_Function.t1_c1, sticky=tk.W)
        ttk.Button(self.tab1, text='Open', command=lambda:self.open_s_calc_file('Open Calculated Solid files')).grid(row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c1, sticky=tk.W)
        ttk.Button(self.tab1, text='Open', command=lambda:self.open_sol_exp_file('Open Experimental Solution files')).grid(row=Scoring_Function.t1_r3, column=Scoring_Function.t1_c1, sticky=tk.W)
        ttk.Button(self.tab1, text='Open', command=lambda:self.open_s_exp_file('Open Experimental Solution files')).grid(row=Scoring_Function.t1_r4, column=Scoring_Function.t1_c1, sticky=tk.W)
        ttk.Button(self.tab1, text='Open',command=lambda: self.open_output()).grid(row=Scoring_Function.t1_r5, column=Scoring_Function.t1_c1, sticky=tk.W)

        ttk.Label(self.tab1, text='Select atoms to be averaged').grid(row=Scoring_Function.t1_r1, column = Scoring_Function.t1_c3,  sticky=tk.W, padx=20)
        ttk.Label(self.tab1, text='Name of new label').grid(row=Scoring_Function.t1_r3, column = Scoring_Function.t1_c3,  sticky=tk.W, padx=20)
        self.Old_labels = ttk.Entry(self.tab1, width=20)
        self.Old_labels.grid(row=Scoring_Function.t1_r2, column = Scoring_Function.t1_c3,  sticky=tk.W, padx=20)



    def Tab2(self):
        '''
        Tab 2: RMSE
        '''
        ttk.Label(self.tab2, text = 'Calculate RMSE (solid)').grid(row=Scoring_Function.t1_r1, sticky=tk.W)
        # ttk.Label(self.tab2, text='Calculate RMSE (solution)').grid(row=Scoring_Function.t1_r2, sticky=tk.W)
        ttk.Button(self.tab2, text='Proton', command=lambda: self.solid_solid(self.solid_exp, self.solid_calc, 'H')).grid(row=Scoring_Function.t1_r1,column=Scoring_Function.t1_c1,sticky=tk.W)
        ttk.Button(self.tab2, text='Carbon', command=lambda: self.solid_solid(self.solid_exp, self.solid_calc, 'C')).grid(row=Scoring_Function.t1_r1,column=Scoring_Function.t1_c2,sticky=tk.W)
        self.tab2_tree = ttk.Treeview(self.tab2, height=10)
        self.tab2_error_label = ttk.Label(self.tab2)
        self.tab2_error_label.grid(row=Scoring_Function.t1_r2, sticky=tk.W, columnspan=10)

    def Tab3(self):
        '''
        Tab 3: Scoring function for each molecule in unit cell
        '''
        ttk.Label(self.tab3, text = 'Scoring function for each mol: ').grid(row=Scoring_Function.t1_r1, sticky=tk.W)
        ttk.Button(self.tab3, text='Load Proton', command=lambda: self.scoring_function_each_mol(self.df_exp_exp, self.df_calc_calc,  'H')).grid(row=Scoring_Function.t1_r1,column=Scoring_Function.t1_c1,sticky=tk.W)
        ttk.Button(self.tab3, text='Load Carbon', command=lambda: self.scoring_function_each_mol(self.df_exp_exp, self.df_calc_calc, 'C')).grid(row=Scoring_Function.t1_r1,column=Scoring_Function.t1_c2,sticky=tk.W)
        ttk.Label(self.tab3, text = 'Delete Atoms').grid(row=Scoring_Function.t1_r2, sticky=tk.W)
        self.tab3_tree = ttk.Treeview(self.tab3, height=10)
        self.tab3_list_box_proton = tk.Listbox(self.tab3, height=5, width=10,  highlightthickness=0, activestyle='none')
        self.tab3_list_box_proton.grid(row=Scoring_Function.t1_r2,column=Scoring_Function.t1_c1, sticky=tk.S)
        self.tab3_list_box_carbon = tk.Listbox(self.tab3, height=5, width=10,  highlightthickness=0)
        self.tab3_list_box_carbon.grid(row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c2, sticky=tk.S)
        ttk.Label(self.tab3, text = 'Run').grid(row=Scoring_Function.t1_r3, sticky=tk.W)
        self.view_graph_inp_t3 = tk.StringVar()
        ttk.Checkbutton(self.tab3, text='View graph', onvalue='Yes', offvalue='No', variable=self.view_graph_inp_t3).grid(row=Scoring_Function.t1_r3,column=Scoring_Function.t1_c2,sticky=tk.W)
        self.save_graph_inp_t3 = tk.StringVar()
        ttk.Checkbutton(self.tab3, text='Save graph', onvalue='Yes', offvalue='No', variable=self.save_graph_inp_t3).grid(row=Scoring_Function.t1_r3,column=Scoring_Function.t1_c3,sticky=tk.W)
        self.save_table_inp_t3 = tk.StringVar()
        ttk.Checkbutton(self.tab3, text='Save Table', onvalue='Yes', offvalue='No', variable=self.save_table_inp_t3).grid(row=Scoring_Function.t1_r3,column=Scoring_Function.t1_c4,sticky=tk.W)
        ttk.Button(self.tab3, text='Run', command=lambda: self.run_scoring_function_each_mol(self.df_exp_exp_mol, self.df_calc_calc_mol, show_fig=self.view_graph_inp_t3.get(), save_fig=self.save_graph_inp_t3.get(), save_table=self.save_table_inp_t3.get())).grid(row=Scoring_Function.t1_r3,column=Scoring_Function.t1_c1,sticky=tk.W)


    def Tab4(self):
        '''
        Tab 4: Scoring function for asymmetric unit cell
        '''
        ttk.Label(self.tab4, text = 'Scoring function for Z\': ').grid(row=Scoring_Function.t1_r1, sticky=tk.W)
        ttk.Button(self.tab4, text='Load Proton', command=lambda: self.scoring_function_Z(self.df_exp_exp, self.df_calc_calc,  'H')).grid(row=Scoring_Function.t1_r1,column=Scoring_Function.t1_c1,sticky=tk.W)
        ttk.Button(self.tab4, text='Load Carbon', command=lambda: self.scoring_function_Z(self.df_exp_exp, self.df_calc_calc, 'C')).grid(row=Scoring_Function.t1_r1,column=Scoring_Function.t1_c2,sticky=tk.W)
        ttk.Label(self.tab4, text = 'Delete Atoms').grid(row=Scoring_Function.t1_r2, sticky=tk.W)
        self.tab4_tree = ttk.Treeview(self.tab4, height=10)
        self.tab4_list_box_proton = tk.Listbox(self.tab4, height=5, width=10,  highlightthickness=0, activestyle='none')
        self.tab4_list_box_proton.grid(row=Scoring_Function.t1_r2,column=Scoring_Function.t1_c1, sticky=tk.S)
        self.tab4_list_box_carbon = tk.Listbox(self.tab4, height=5, width=10,  highlightthickness=0)
        self.tab4_list_box_carbon.grid(row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c2, sticky=tk.S)
        ttk.Label(self.tab4, text = 'Run').grid(row=Scoring_Function.t1_r3, sticky=tk.W)
        self.view_graph_inp = tk.StringVar()
        ttk.Checkbutton(self.tab4, text='View graph', onvalue='Yes', offvalue='No', variable=self.view_graph_inp).grid(row=Scoring_Function.t1_r3,column=Scoring_Function.t1_c2,sticky=tk.W)
        self.save_graph_inp = tk.StringVar()
        ttk.Checkbutton(self.tab4, text='Save graph', onvalue='Yes', offvalue='No', variable=self.save_graph_inp).grid(row=Scoring_Function.t1_r3,column=Scoring_Function.t1_c3,sticky=tk.W)
        self.save_table_inp = tk.StringVar()
        ttk.Checkbutton(self.tab4, text='Save Table', onvalue='Yes', offvalue='No', variable=self.save_table_inp).grid(row=Scoring_Function.t1_r3,column=Scoring_Function.t1_c4,sticky=tk.W)
        ttk.Button(self.tab4, text='Run', command=lambda: self.run_scoring_function_z(self.df_exp_exp2, self.df_calc_calc2, show_fig=self.view_graph_inp.get(), save_fig=self.save_graph_inp.get(), save_table=self.save_table_inp.get())).grid(row=Scoring_Function.t1_r3,column=Scoring_Function.t1_c1,sticky=tk.W)


    def get_average(self):
        old_labels_str = str(self.Old_labels.get())
        old_labels = old_labels_str.split(',')
        print(self.New_labels)
        if self.solid_calc is not None:
            try:
                Scoring_Function.average_H(self.df_calc_calc, old_labels, self.New_labels)
                Scoring_Function.average_H(self.solid_calc, old_labels, self.New_labels)
                self.view_input.config(text=self.solid_calc.loc[self.New_labels, :])
                self.label_input.config(text='{} = ({})'.format(self.New_labels, ', '.join(old_labels)))
            except:
                Message(title='Error!', message='Atom {} not found'.format(old_labels_str), master=self.master).show()
                print('Labels not in table')
                self.label_input.config(text='Atoms ({}) not found'.format(old_labels_str))


        print(old_labels)

    def open_sol_calc_file(self, title):
        sol_calc_file = askopenfile(title=title, mode='r', filetypes=(('CSV Files', '*.csv'),))
        self.solution_calc = self.open_files(sol_calc_file)
        if self.solution_calc is not None:
            ttk.Label(self.tab1, background='#000fff000', width=2, anchor=tk.W).grid(row=Scoring_Function.t1_r1,
                                                                              column=Scoring_Function.t1_c2)
        self.update_buttons()
    def open_s_calc_file(self, title):
        s_calc_file = askopenfile(title=title, mode='r', filetypes=(('CSV Files', '*.csv'),))
        self.solid_calc = self.open_files(s_calc_file)
        if self.solid_calc is not None:
            ttk.Label(self.tab1, background='#000fff000', width=2, anchor=tk.W).grid(row=Scoring_Function.t1_r2,
                                                                              column=Scoring_Function.t1_c2)
        self.update_buttons()
    def open_sol_exp_file(self, title):
        sol_exp_file = askopenfile(title=title, mode='r', filetypes=(('CSV Files', '*.csv'),))
        self.solution_exp = self.open_files(sol_exp_file)
        if self.solution_exp is not None:
            ttk.Label(self.tab1, background='#000fff000', width=2, anchor=tk.W).grid(row=Scoring_Function.t1_r3,
                                                                              column=Scoring_Function.t1_c2)
        self.update_buttons()
    def open_s_exp_file(self, title):
        s_exp_file = askopenfile(title=title, mode='r', filetypes=(('CSV Files', '*.csv'),))
        self.solid_exp = self.open_files(s_exp_file)
        if self.solid_exp is not None:
            ttk.Label(self.tab1, background='#000fff000', width=2, anchor=tk.W).grid(row=Scoring_Function.t1_r4,
                                                                              column=Scoring_Function.t1_c2)
        self.update_buttons()
    def open_files(self, file):
        if file != None:
            self.in_file.config(text='Opened: %s'%file.name)
            df = pd.read_csv(file, index_col=0)
            self.view_input.config(text = df.head(5))
            return df
        else:
            Message(title='Error!', message='Error no file selected', master=self.master).show()
            self.in_file.config(text='No file selected')
            self.view_input.config(text='')
            return None
        self.update_buttons()

    def open_output(self):
        dir = askdirectory()
        self.path = glob.glob(dir)[0]
        print(self.path)
        if self.path != '':
            print(self.path)
            ttk.Label(self.tab1, background='#000fff000', width=2, anchor=tk.W).grid(row=Scoring_Function.t1_r5,
                                                                                     column=Scoring_Function.t1_c2)
            self.in_file.config(text='Output directory: %s' % self.path)
        else:
            Message(title='Error!', message='Error no directory selected', master=self.master).show()
            self.in_file.config(text='No directory selected')
            self.view_input.config(text='')

    def solid_solid(self,df_exp, df_calc, Atom):
        '''
        Compare solid exp vs Solid Calc cia RMSE
        Try to see if you can differentiate structuers via solid only
        '''
        df_calc = df_calc[df_calc.index.str.contains(Atom)]
        df_exp = df_exp[df_exp.index.str.contains(Atom)]
        df_1 = pd.DataFrame()
        for col in df_exp.columns:
            conformer = []
            exp_conf = []
            r_sq = []
            m = []
            c = []
            rmse_l = []
            print(col)
            for i in df_calc.columns:
                solid_df = pd.concat([df_calc[i], df_exp], axis=1)
                x = solid_df.loc[:, col]
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, solid_df.iloc[:, 0])
                r_sqrd = r_value ** 2
                m.append(slope.round(2))
                c.append(intercept.round(2))
                r_sq.append(r_sqrd.round(2))
                conformer.append(i)
                exp_conf.append(col)
                pred_y = slope * x + intercept
                rmse = Scoring_Function.get_rmse(pred_y, solid_df.iloc[:, 0])
                rmse_l.append(rmse.round(2))
            df = pd.DataFrame(data={'Conformer': conformer, 'Exp': exp_conf,
                                    'R_sq': r_sq, 'Gradient': m, 'Intercept': c,
                                    'RMSE': rmse_l})
            self.tab2_error_label.config(text='Atoms used: {}'.format(solid_df.index.to_list()))
            if df.isnull().values.any():
                print('Mismatch between calculated labels and experimental labels')
                continue
            else:
                # outfile = r'{}_Solid_exp_Vs_Solid_calc_{}.csv'.format(col, Atom)
                df_1 = df_1.append(df)
                # if self.path !=None:
                #     df_1.to_csv(self.path + outfile, index=None)
        if df_1.empty:
            print('Mismatch between calculated labels and experimental labels')
            self.tab2_error_label.config(text='Mismatch between calculated labels and experimental labels see below')
            calc_labels = df_calc.index.to_list()
            exp_labels = df_exp.index.to_list()
            df_error = pd.concat([pd.Series(calc_labels),pd.Series(exp_labels)], ignore_index=True, axis=1)
            df_error.columns = ['Calculated', 'Experimental']
            Scoring_Function.tree_print(df_error, self.tab2_tree)
        else:
            df_1.sort_values(['Exp','RMSE'], inplace=True)  # Sort by RMSE values
            Scoring_Function.tree_print(df_1, self.tab2_tree)
            if self.path != None:
                df_1.to_csv(self.path + r'/Solid_exp_Vs_Solid_calc_{}.csv'.format(Atom), index=None)
                print('saving')
                print(self.path + r'/Solid_exp_Vs_Solid_calc_{}.csv')
    def scoring_function_Z(self,df_exp, df_calc, Atom):
        '''
        :param df1: Exp solution
        :param df2: Exp solid
        :param df3: Calc solution
        :param df4: Calc solid
        :param Atom: H or C
        '''
        df_exp_sub = df_exp.copy()
        df_calc_sub = df_calc.copy()
        self.df_calc_calc2 = df_calc[df_calc.index.str.contains(Atom)] #Filter by atom, need to create new df
        self.df_exp_exp2 = df_exp[df_exp.index.str.contains(Atom)]
        df_inputs = pd.concat([self.df_exp_exp2, self.df_calc_calc2], axis=1, join="inner")
        df_inputs.dropna(inplace=True)
        l_inputs = df_inputs.index.to_list()

        if Atom == 'H':
            self.tab4_list_box_proton.delete('0', 'end')
            self.tab4_list_box_carbon.delete('0', 'end')
            for atom in l_inputs:
                self.tab4_list_box_proton.insert(tk.END,atom)
                ttk.Button(self.tab4, text='Delete', command=lambda: self.delete_item(self.tab4_list_box_proton)).grid(
                    row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c3, sticky=tk.W)
                ttk.Button(self.tab4, text='Reset', command=lambda: self.reset_list(self.tab4_list_box_proton, l_inputs)).grid(
                    row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c4, sticky=tk.W)

        else:
            self.tab4_list_box_proton.delete('0', 'end')
            self.tab4_list_box_carbon.delete('0', 'end')
            for atom in l_inputs:
                self.tab4_list_box_carbon.insert(tk.END,atom)
                ttk.Button(self.tab4, text='Delete', command=lambda: self.delete_item(self.tab4_list_box_carbon)).grid(
                    row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c3, sticky=tk.W)
                ttk.Button(self.tab4, text='Reset', command=lambda: self.reset_list(self.tab4_list_box_carbon, l_inputs)).grid(
                    row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c4, sticky=tk.W)

    def update_listbox(self):
        for item in self.tab4_list_box_proton.get(0, tk.END):
            self.atoms_selec.append(item)
        for item in self.tab3_list_box_proton.get(0, tk.END):
            self.atoms_selec.append(item)
        for item in self.tab4_list_box_carbon.get(0, tk.END):
            self.atoms_selec.append(item)
        for item in self.tab3_list_box_carbon.get(0, tk.END):
            self.atoms_selec.append(item)

    def run_scoring_function_z(self, df_exp, df_calc,show_fig = 'No', save_fig = 'No', save_table= 'No'):
        '''
        Run scoring function for each molecule in assymetric unit
        '''
        self.atoms_selec = []
        name = []
        rsq = []
        ml = []
        cl = []
        pl = []
        exp_poly = []
        self.update_listbox()
        print(df_calc)
        print(df_exp)
        if 'H' in str(self.atoms_selec):
            Atom ='H'
        if 'C' in str(self.atoms_selec):
            Atom = 'C'
        print(Atom)
        df_calc = df_calc[df_calc.index.isin(self.atoms_selec)]
        df_print = pd.DataFrame()
        for mol in range(0,len(df_exp.columns)):
            for i in range(len(df_calc.columns)):
                y = df_calc.iloc[:, i]
                name.append(y.name)
                y_title = y.name
                y.name = 'y'
                x = df_exp.iloc[:, mol]
                x_title = x.name
                x.name = 'x'
                df_calc_exp = pd.concat([x, y], axis=1)
                df_calc_exp.dropna(axis=0, inplace=True)
                df_calc_exp['Atom'] = df_calc_exp.index
                slope, intercept, r_value, p_value, std_err = stats.linregress(df_calc_exp['x'], df_calc_exp['y'])
                r_sqr = r_value ** 2
                p_value = p_value / 2  # for one-tailed
                n = len(x)  # number of samples
                ttest = (r_value * np.sqrt(n - 2) / np.sqrt(1 - r_sqr))  # ttest for two unkowns
                if ttest < 0:
                    p_value = 1 - p_value  # 'greater'

                title = "Scoring function plot of calc %s vs exp %s" % (y_title, x_title)
                fig = px.scatter(df_calc_exp, x=x.name, y=y.name, opacity=0.65, text=df_calc_exp.index,
                                 trendline='ols', hover_data={'x': False, 'y': False, 'Atom': True})
                fig.update_layout(
                    title=title,
                    xaxis_title="Experimental solid-solution",
                    yaxis_title="Calculated solid-solution",
                )
                if show_fig == 'Yes':
                    fig.show()
                if save_fig == 'Yes':
                    Scoring_Function.save_figures(fig, self.path + r'/figures/%s/%s/' % (x_title, Atom) + title + '_' + Atom + '.html')
                print('r_squared = {:.2f} \n m = {:.2f} \n c = {:.2f} \n p = {:.3f}'.format(r_value, slope, intercept,
                                                                                            p_value))
                rsq.append(round(r_sqr, 2))
                ml.append(round(slope, 2))
                cl.append(round(intercept, 2))
                pl.append(round(p_value, 5))
                exp_poly.append(x_title)
            df = pd.DataFrame(data={'Polymorph': name, 'Experimental': exp_poly, 'r_sq': rsq, 'm': ml, 'c': cl, 'p': pl})
            df.sort_values(by='p', inplace=True)
            print(df)
            df_print = df_print.append(df)
            if save_table=='Yes':
                df.to_csv(self.path + r'/scoring_function_1000_%s_%s.csv' % (x_title, Atom), index=None)
        df_print.sort_values('p', inplace=True) # Sort by p values
        df_out = df_print.groupby(['Experimental', 'p'], as_index=False).head() # Group the exp & p values and get the first values
        col = df_out.pop('p') # Remove the p
        col_length= len(df_out.columns) # get last col
        df_out.insert(col_length,col.name,col) # insert p into last col
        Scoring_Function.tree_print(df_out, self.tab4_tree) # show in gui

    def callback_new_label(self,selection):
        '''
        Selects value from dropdown
        '''
        self.New_labels = selection

    def scoring_function_each_mol(self,df_exp, df_calc, Atom):
        '''
        :param df1: Exp solution
        :param df2: Exp solid
        :param df3: Calc solution
        :param df4: Calc solid
        :param Atom: H or C
        '''
        df_exp_sub = df_exp.copy()
        df_calc_sub = df_calc.copy()
        self.df_calc_calc_mol = df_calc[df_calc.index.str.contains(Atom)]
        self.df_exp_exp_mol = df_exp[df_exp.index.str.contains(Atom)]
        df_inputs = pd.concat([self.df_exp_exp_mol, self.df_calc_calc_mol], axis=1, join="inner")
        df_inputs.dropna(inplace=True)
        l_inputs = df_inputs.index.to_list()
        if Atom == 'H':
            self.tab3_list_box_proton.delete('0', 'end')
            self.tab3_list_box_carbon.delete('0', 'end')
            for atom in l_inputs:
                self.tab3_list_box_proton.insert(tk.END,atom)
                ttk.Button(self.tab3, text='Delete', command=lambda: self.delete_item(self.tab3_list_box_proton)).grid(
                    row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c3, sticky=tk.W)
                ttk.Button(self.tab3, text='Reset', command=lambda: self.reset_list(self.tab3_list_box_proton, l_inputs)).grid(
                    row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c4, sticky=tk.W)
        else:
            self.tab3_list_box_proton.delete('0', 'end')
            self.tab3_list_box_carbon.delete('0', 'end')
            for atom in l_inputs:
                self.tab3_list_box_carbon.insert(tk.END,atom)
                ttk.Button(self.tab3, text='Delete', command=lambda: self.delete_item(self.tab3_list_box_carbon)).grid(
                    row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c3, sticky=tk.W)
                ttk.Button(self.tab3, text='Reset', command=lambda: self.reset_list(self.tab3_list_box_carbon, l_inputs)).grid(
                    row=Scoring_Function.t1_r2, column=Scoring_Function.t1_c4, sticky=tk.W)


    def run_scoring_function_each_mol(self,df_exp, df_calc, show_fig = 'No', save_fig = 'No', save_table= 'No'):
        '''
        This one combines the two molecules into 1 for exp and calc, because just combining exp didnt really work
        '''
        self.atoms_selec = []
        self.update_listbox()
        print(self.atoms_selec)
        if 'H' in str(self.atoms_selec):
            Atom ='H'
        if 'C' in str(self.atoms_selec):
            Atom = 'C'
        df_calc = df_calc[df_calc.index.isin(self.atoms_selec)]
        # df_exp = df_exp[df_exp.index.str.contains(Atom)]
        # print(df_calc)
        # print(df_exp)
        if len(df_exp.columns) > 1:
            select_col = df_exp[df_exp.filter(like='mol').columns]
            exp_name = [s for s in select_col.columns if '_' in s]
            exp_name = [i.split('_', 1)[0] for i in exp_name]
            exp_name = exp_name[0]
            x = Scoring_Function.join_mols(df_exp, len(df_exp.columns))
            x_title = 'Exp_({})'.format(exp_name)
            print(x_title)
            x.name = 'x'
            name = []
            rsq = []
            ml = []
            cl = []
            pl = []
            # Combines the two molecules into 1
            df_z_2 = df_calc[df_calc.filter(like='mol').columns]
            mul_forms = [s for s in df_z_2.columns if '_' in s]
            mul_forms = [i.split('_', 1)[0] for i in mul_forms]
            each_form = list(dict.fromkeys(mul_forms))
            for i in range(len(each_form)):
                df = df_z_2.filter(regex=each_form[i])
                y = Scoring_Function.join_mols(df, len(df_exp.columns))
                name.append(each_form[i])
                y_title = each_form[i]
                y.name = 'y'
                df_calc_exp = pd.concat([x, y], axis=1)
                df_calc_exp.dropna(axis=0, inplace=True)
                df_calc_exp['Atom'] = df_calc_exp.index
                slope, intercept, r_value, p_value, std_err = stats.linregress(df_calc_exp['x'], df_calc_exp['y'])
                r_sqr = r_value ** 2
                p_value = p_value / 2  # for one-tailed
                n = len(x)  # number of samples
                ttest = (r_value * np.sqrt(n - 2) / np.sqrt(1 - r_sqr))  # ttest for two unkowns
                if ttest < 0:
                    p_value = 1 - p_value  # 'greater'

                title = "Scoring function plot of calc %s vs exp %s" % (y_title, x_title)
                fig = px.scatter(df_calc_exp, x=x.name, y=y.name, opacity=0.65, text=df_calc_exp.index,
                                 trendline='ols', hover_data={'x': False, 'y': False, 'Atom': True})
                fig.update_layout(
                    title=title,
                    xaxis_title="Experimental solid-solution",
                    yaxis_title="Calculated solid-solution",
                )
                if show_fig == 'Yes':
                    fig.show()
                if save_fig == 'Yes':
                    Scoring_Function.save_figures(fig, self.path + r'/figures/%s/%s/' % (x_title, Atom) + title + '_' + Atom + '.html')
                # print('r_squared = {:.2f} \n m = {:.2f} \n c = {:.2f} \n p = {:.3f}'.format(r_value, slope, intercept, p_value))
                rsq.append(round(r_sqr, 2))
                ml.append(round(slope, 2))
                cl.append(round(intercept, 2))
                pl.append(round(p_value, 5))
            df = pd.DataFrame(data={'Polymorph': name, 'r_sq': rsq, 'm': ml, 'c': cl, 'p': pl})
            print(df)
            if save_table == 'Yes':
                df.to_csv(self.path + r'/scoring_function_1000_%s_%s.csv' % (x_title, Atom), index=None)
            Scoring_Function.tree_print(df, self.tab3_tree)  # show in gui
    @staticmethod
    def drop_first_col(df):
        return df.drop(df.columns[[0]], axis=1, inplace=True)

    @staticmethod
    def average_H(df, l, x):
        '''
        :param df: dataframe to be averaged
        :param l: list of protons to be averaged
        :param x: new assignemnt from exp
        Need to average protons that are avereaged out experimentally
        :return:
        '''
        df.loc[x] = df.loc[l].mean()
        df.drop(index=l, inplace=True)
        return df

    @staticmethod
    def get_rmse(prediction, targets):
        '''
        Calculates RMSE by passing predicted and actual values
        '''
        difference = prediction - targets
        difference_sq = difference ** 2
        difference_sq_mean = difference_sq.mean()
        rmse_val = np.sqrt(difference_sq_mean)
        return rmse_val

    @staticmethod
    def save_figures(fig, filename):
        if not os.path.exists(os.path.dirname(filename)):
            try:
                os.makedirs(os.path.dirname(filename))
                fig.write_html(filename)
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        else:
            fig.write_html(filename)

    @staticmethod
    def calc_shifts(df1, df2):
        '''
        DFT calculated
        :df1: solution
        :df2: solid
        :return: Difference df of solution and solid
        '''
        df_calc = pd.concat([df1['Average'], df2], axis=1)
        df_calc.dropna(axis=0, inplace=True)  # drops rows that dont exits in solid or solution
        df_calc = df_calc.subtract(df_calc.iloc[:, 0], axis=0)
        Scoring_Function.drop_first_col(df_calc)
        return df_calc

    @staticmethod
    def exp_shifts(df1,df2):
        '''
        Experimental
        :df1: solution
        :df2: solid
        :return: Difference df of solution and solid
        '''
        df_exp = pd.concat([df1, df2], axis=1)
        df_exp = df_exp.subtract(df_exp.iloc[:, 0], axis=0)
        Scoring_Function.drop_first_col(df_exp) # Drops extra column in soliution
        return df_exp

    @staticmethod
    def drop_first_col(df):
        return df.drop(df.columns[[0]], axis=1, inplace=True)

    @staticmethod
    def join_mols(df, nmols):
        '''
        :param df: df exp or df calc
        :param nmols: number of molecules in df exp
        '''
        df_z_2 = df[df.filter(like='mol').columns]
        mul_forms = [s for s in df_z_2.columns if '_' in s]
        mul_forms = [i.split('_', 1)[1] for i in mul_forms]
        each_form = list(dict.fromkeys(mul_forms))
        if len(each_form) > 0 and len(each_form) <= nmols:  # If there are more than one molecule
            dfn = []  # Store all the molecules
            for i in range(0, nmols):
                nme = each_form[i].split('mol')[1]  # split them by molecule
                dftmp = df.iloc[:, i]  # select each molecule
                dftmp = dftmp.add_prefix('%s_' % nme)  # add a prefix to index to identify each molecule
                dfn.append(dftmp)  # combine the molecules
        df2 = pd.concat(dfn, axis=0)  # combine the whole molecule to output
        return df2

    @staticmethod
    def delete_item(my_list):
        my_list.delete(tk.ANCHOR)

    @staticmethod
    def reset_list(my_list, l_inputs):
        my_list.delete('0', 'end')
        for atom in l_inputs:
            my_list.insert(tk.END, atom)


    @staticmethod
    def clear_tree(my_tree):
        '''
        Clears old tree
        '''
        my_tree.delete(*my_tree.get_children())


    @staticmethod
    def tree_print(df, my_tree):
        Scoring_Function.clear_tree(my_tree)
        my_tree['column'] = df.columns.to_list()
        my_tree['show'] = 'headings'
        for column in my_tree['column']:
            my_tree.column(str(column), width=100, minwidth=25, stretch =False, anchor='center')
            my_tree.heading(column, text=column)
        df_rows = df.to_numpy().tolist()
        count = 0
        for row in df_rows:
            my_tree.insert('', 'end', iid = count, values=row)
            count +=1
        my_tree.grid(row=Scoring_Function.t1_r7, columnspan=40)


    @staticmethod
    def switch(b1):
        if b1["state"] == tk.NORMAL:
            b1["state"] = tk.DISABLED
        else:
            b1["state"] = tk.NORMAL


root = tk.Tk()
style = ttk.Style()
style.theme_use('classic')
root.geometry('600x600')
root.title('Scoring Function')
scoring_function = Scoring_Function(root)
root.mainloop()



