import os, glob
import pandas as pd
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfile, askdirectory
from tkinter.messagebox import Message

'''
Get solution files
'''

def open_dir():
    global in_dir, pass_label, mag_files
    dir = askdirectory()
    if dir is not None:
        mag_files = glob.glob(os.path.join(dir, '*.magres'))
        in_dir.config(text=dir)
        in_dir.grid(row=7,columnspan=5, sticky=tk.S)
        if not mag_files:
            print('Error no Magres file found in %s' % dir)
            Message(title='Error!', message='Error no Magres file found in %s' % dir, master=root).show()
            pass_label.config(text ='Error no Magres file found in %s' % dir)
            pass_label.grid(row = 8,columnspan=5, sticky = tk.S)
        else:
            print(mag_files)
            pass_label.config(text ='A total of %s magres files found in %s' %(len(mag_files),dir))
            pass_label.grid(row = 8,columnspan=5, sticky = tk.S)


def open_file():
    global in_file, mag_labels
    file = askopenfile(title = 'Select file containing magres labels', mode ='r', filetypes = (('CSV Files','*.csv'),))
    in_file.config(text=file.name)
    in_file.grid(row=7,columnspan=5, sticky=tk.S)
    if file is not None:
        mag_labels = pd.read_csv(file, index_col=0)
        pass_label.config(text=mag_labels.head())
        pass_label.grid(row=7,columnspan=5, sticky=tk.S)
        print(mag_labels)


def getInput():
    global Proton_Get, proton_label, Carbon_Get, carbon_label, Proton_ref, Carbon_ref
    Proton_ref = float(Proton_Get.get())
    proton_label.config(text=Proton_ref)
    proton_label.grid(row=0, column=2, sticky=tk.S)
    Carbon_ref = float(Carbon_Get.get())
    carbon_label.config(text=Carbon_ref)
    carbon_label.grid(row=1,column=2, sticky=tk.S)


def solution_magres_output(labels,select_files):
    solution_labels = labels.copy().iloc[:, :1]
    solution_labels['atom'] = solution_labels.iloc[:, 0].str[:1]
    array = []
    for filename in select_files:  # itinerate over every file in the s
        fname = os.path.basename(filename)  # Get file name from each .magres file
        with open(filename, "r") as f:
            name_list = []
            shifts_list = []
            df = pd.DataFrame()
            for line in f:
                if ("Isotropic" not in line):
                    continue
                ki = line.split()
                name = "%s%s" % (ki[0], ki[1])
                shift = float(ki[3])
                if (name[0] == 'H'):
                    shift = Proton_ref - shift
                elif (name[0] == 'C'):
                    shift = Carbon_ref - shift
                name_list.append(name)
                shifts_list.append(shift)
            df['Solution'] = name_list
            try:
                df['{}'.format(fname)] = shifts_list
            except:
                print('Error! Atom %s not found in %s' % (name, fname))
            df.set_index('Solution')
            array.append(df)
    df_combined = pd.concat(array, axis=1)
    df_combined = df_combined.loc[:, ~df_combined.columns.duplicated()]
    dt = df_combined.set_index('Solution').T.to_dict('dict')
    d = []
    solution_labels.iloc[:, 0].apply(lambda x: d.append(dt[x]))
    df = pd.DataFrame.from_dict(d)
    df.index = solution_labels.index
    # print(df)
    return df

def out_folder():
    global out_file, path
    dir = askdirectory()
    path = str(glob.glob(dir))[1:-1]
    out_file.config(text=path)
    out_file.grid(row = 9,columnspan=5, sticky = tk.S)

def getOutput():
    global Output_name, mag_labels, mag_files
    getInput()
    out_name = Output_name.get()
    print('Output file name %s'%out_name)
    path1 = os.path.join(path[1:-1], out_name+'.csv')
    print(path1)
    df_out = solution_magres_output(mag_labels, mag_files)
    df_out.to_csv(path1 )
    pass_label.config(text ='File saved in {}\n{}\n{}'.format(path1, df_out.head(3), df_out.tail(3)))
    pass_label.grid(row = 10, columnspan=5, sticky = tk.S)

'''
Solid files
'''

def open_solid():
    global solid_files
    dir = askdirectory()
    if dir is not None:
        solid_files = glob.glob(os.path.join(dir, '*.magres'))
        in_dir.config(text=dir)
        in_dir.grid(row=7,columnspan=5, sticky=tk.S)
        if not solid_files:
            print('Error no Magres file found in %s' % dir)
            Message(title='Error!', message='Error no Magres file found in %s' % dir, master=root).show()
            pass_label_2.config(text ='Error no Magres file found in %s' % dir)
            pass_label_2.grid(row = 9,columnspan=5, sticky = tk.S)
        else:
            print(solid_files)
            pass_label_2.config(text ='A total of %s magres files found in %s' %(len(solid_files),dir))
            pass_label_2.grid(row = 9,columnspan=5, sticky = tk.S)

def open_solid_file():
    global solid_labels
    file = askopenfile(title = 'Select file containing magres labels', mode ='r', filetypes = (('CSV Files','*.csv'),))
    in_file.config(text=file.name)
    in_file.grid(row=8,columnspan=5, sticky=tk.S)
    if file is not None:
        solid_labels = pd.read_csv(file, index_col=0)
        pass_label_2.config(text=solid_labels.head())
        pass_label_2.grid(row=7,columnspan=5, sticky=tk.S)
        print(solid_labels)


def getInput_solid():
    global Proton_Get_solid, proton_label_solid, \
        Carbon_Get_solid, carbon_label_solid, Proton_ref_solid, Carbon_ref_solid
    Proton_ref_solid = float(Proton_Get_solid.get())
    proton_label_solid.config(text=Proton_ref_solid)
    proton_label_solid.grid(row=0, column=2, sticky=tk.S)
    Carbon_ref_solid = float(Carbon_Get_solid.get())
    carbon_label_solid.config(text=Carbon_ref_solid)
    carbon_label_solid.grid(row=1,column=2, sticky=tk.S)

def solid_magres_output(labels,select_files):
    '''
    :param labels: Solid labels that match a reference list
    :param select_files: Solid magres lables
    :return: Chemical shift for each file
    Need the column of solid labels to match the files exactly without the extenstion
    Do not use periods in file names except for the .magres at the end, there is a '.'seperator
    '''
    solid_labels =labels.copy()
    # original_label = solid_labels.columns
    print(solid_labels.columns)
    df = pd.DataFrame()
    for filename in select_files:  # itinerate over every file in the s
        fname = os.path.basename(filename)  # Get file name from each .magres file
        fname_noext = fname.split('.')[0]
        print(fname_noext)
        tmpdf = solid_labels.filter(regex=fname_noext) # Take subset of df for each file matching the labels
        with open(filename, "r") as f:
            d = {}
            for line in f:
                if ("Isotropic" not in line):
                    continue
                ki = line.split()
                name = "%s%s" % (ki[0], ki[1])
                shift = float(ki[3])
                if (name[0] == 'H'):
                    shift = Proton_ref_solid - shift
                elif (name[0] == 'C'):
                    shift = Carbon_ref_solid - shift
                d[name] = shift
            tmpdf2 = tmpdf.replace(d) #Match labels and replace with float
            df = pd.concat([df,tmpdf2], axis=1) #Combine into one df
    print(df)
    return df

def get_Solid_Output():
    global solid_labels, solid_files
    getInput_solid()
    out_name = Output_name_solid.get()
    print(out_name)
    path1 = os.path.join(path[1:-1], out_name+'.csv')
    df_out = solid_magres_output(solid_labels, solid_files)
    df_out.to_csv(path1)
    pass_label_2.config(text ='File saved in {}\n{}\n{}'.format(path1, df_out.head(3), df_out.tail(3)))
    pass_label_2.grid(row = 14, columnspan=5, sticky = tk.S)

'''
Tkinter
'''

def main():
    global pass_label, proton_label, carbon_label, \
        in_dir, in_file, out_file, Proton_Get, Carbon_Get, \
        Output_name, root, Output_name_solid,\
        Proton_Get_solid, Carbon_Get_solid, proton_label_solid, carbon_label_solid, pass_label_2
    root = tk.Tk()
    root.title('Chemical Shifts')
    root.geometry('600x400')
    style = ttk.Style()
    style.theme_use('classic')
    tabControl = ttk.Notebook(root)
    tab1 = ttk.Frame(tabControl)
    tab2 = ttk.Frame(tabControl)
    tabControl.add(tab1, text='Solution')
    tabControl.add(tab2, text='Solid')
    tabControl.pack(expand=1, fill='both')

    '''
    tab 1 : solution
    '''
    ttk.Label(tab1, text='Proton reference').grid(row=0, sticky=tk.W)
    ttk.Label(tab1, text='Carbon reference').grid(row=1, sticky=tk.W)
    ttk.Label(tab1, text='Solution labels').grid(row=2, sticky=tk.W)
    ttk.Label(tab1, text = 'Solution file directory').grid(row = 3, sticky = tk.W)
    ttk.Label(tab1, text='Output').grid(row=4, sticky=tk.W)

    Proton_Get = ttk.Entry(tab1, width =5)
    Proton_Get.insert(0, 30)
    Carbon_Get = ttk.Entry(tab1, width =5)
    Carbon_Get.insert(0, 170)
    Output_name = ttk.Entry(tab1, width =10)

    pass_label = ttk.Label(tab1)
    proton_label = ttk.Label(tab1)
    carbon_label = ttk.Label(tab1)
    in_dir = ttk.Label(tab1)
    in_file = ttk.Label(tab1)
    out_file = ttk.Label(tab1)

    Proton_Get.grid(row=0, column=1)
    Carbon_Get.grid(row=1, column=1)
    Proton_Get.bind('<Return>', (lambda event: getInput()))
    Carbon_Get.bind('<Return>', (lambda event: getInput()))
    Output_name.grid(row=4, column=1)

    ttk.Button(tab1, text='Open', command=lambda: open_file()).grid(row=2, column=1, sticky=tk.W)  #Solu labels
    ttk.Button(tab1, text ='Open', command = lambda:open_dir()).grid(row = 3, column = 1, sticky = tk.W) # Solu magres
    ttk.Button(tab1, text="Run", command=getOutput).grid(row=4, column=3, sticky=tk.W) # Solution Run
    ttk.Button(tab1, text='Folder', command=lambda: out_folder()).grid(row=4, column=2, sticky=tk.W)
    '''
    tab 2: solid
    '''
    pass_label_2 = ttk.Label(tab2)
    proton_label_solid = ttk.Label(tab2)
    carbon_label_solid = ttk.Label(tab2)

    ttk.Label(tab2, text='Proton reference').grid(row=0, sticky=tk.W)
    ttk.Label(tab2, text='Carbon reference').grid(row=1, sticky=tk.W)
    ttk.Label(tab2, text='Solid labels').grid(row=2, sticky=tk.W)
    ttk.Label(tab2, text = 'Solid files directory').grid(row=3, sticky=tk.W)
    ttk.Label(tab2, text='Output').grid(row=4, sticky=tk.W)

    Proton_Get_solid = ttk.Entry(tab2, width =5)
    Proton_Get_solid.insert(0, 30)
    Carbon_Get_solid = ttk.Entry(tab2, width =5)
    Carbon_Get_solid.insert(0, 170)
    Output_name_solid = ttk.Entry(tab2, width =10)

    Proton_Get_solid.grid(row=0, column=1)
    Carbon_Get_solid.grid(row=1, column=1)
    Proton_Get_solid.bind('<Return>', (lambda event: getInput()))
    Carbon_Get_solid.bind('<Return>', (lambda event: getInput()))
    Output_name_solid.grid(row=4, column=1)

    ttk.Button(tab2, text ='Open', command = lambda:open_solid_file()).grid(row = 2, column = 1, sticky = tk.W) # Solid labels
    ttk.Button(tab2, text ='Open', command = lambda:open_solid()).grid(row = 3, column = 1, sticky = tk.W) # Solid magres
    ttk.Button(tab2, text="Run", command=get_Solid_Output).grid(row=4, column=3, sticky=tk.W) # Solid Run
    ttk.Button(tab2, text='Folder', command=lambda: out_folder()).grid(row=4, column=2, sticky=tk.W)

    root.mainloop()
main()


