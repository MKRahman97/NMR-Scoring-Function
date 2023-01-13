import os, glob
import pandas as pd
import numpy as np
import regex as re
import plotly.express as px
import plotly.graph_objects as go
import scipy.stats as stats
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfile, askdirectory
from tkinter.messagebox import Message

def read_df():
    '''
    Load the chemical shift data
    '''
    file = askopenfile(title='Select chemical shift data', mode='r', filetypes=(('CSV Files', '*.csv'),))
    if file is not None:
        df = pd.read_csv(file.name, index_col=0)
        try:
            global df_out
            df_out = df.astype(float)
            chem_shift_path.config(text=file.name)
            chem_shift_path.grid(row=2, sticky=tk.W,columnspan=8)
            print_lab.config(text=df.head())
            print_lab.grid(row=10, sticky=tk.W, columnspan=8)

            #Atoms for CLT and Bootstrapping
            get_atom_list_CLT(df)
            get_atom_list_boots(df)
            #Atoms for DOFvs DOF plot
            DOF_Atom(df_out)
            Atom_table(df_out)


        except:
            print('Error! please select the correct data frame')
            Message(title='Error!', message='Error! %s does not contain chemical shift data, please select the correct dataframe' %file.name,
                    master=root).show()
            chem_shift_path.config(text='Error! no data selected from %s'%file.name)
            chem_shift_path.grid(row=2, sticky=tk.W,columnspan=8)
    else:
        Message(title='Error!', message='Error no chemical shift data selected', master=root).show()

def get_DOF_data():
    '''
    Load the torsion data
    '''
    file = askopenfile(title='Select Torsion data', mode='r', filetypes=(('Text Files', '*.txt'),)) #'All Files', '*.*'
    if file is not None:
        try:
            torsion_path.config(text=file.name)
            torsion_path.grid(row=4, sticky=tk.W, columnspan=8)
            with open(file.name,'r') as f:
                reader = f.read()
                df = pd.DataFrame()
                for i, part in enumerate(reader.split('}')[1:]):
                    lines = part.split('\n')
                    data = [x for x in lines if any(kw in x for kw in 'SDF')]
                    try:
                        df[str(i+1)] = data
                    except:
                        df2 = pd.DataFrame()
                        df2[str(i+1)] = data
                        df = pd.concat([df,df2], axis=1)
                    df['DOF_' +str(i+1)] = df[str(i+1)].str.extract('(\d+\.\d+)')
                df.dropna(inplace=True)
                df['index_col'] = df['1'].str.extract('(\d+)').astype(int)
                df = df.set_index('index_col')
                df = df.loc[:, df.columns.str.startswith('DOF')]
                df = df.astype(float)
                # Centre points around zero
                df[df > int(180)] = df[df > int(180)] - int(360)
            global DOF_df
            DOF_df = df
            print_lab.config(text=df.head())
            print_lab.grid(row=10, sticky=tk.W, columnspan=8)

            #Create DOF list for Dof vs chem shift plot
            Dof_list_energy_plot(DOF_df)
            Dof_X(DOF_df)
            Dof_Y(DOF_df)

            # Dof list for selecting angles
            Dof_angle(DOF_df)

            # Select colour for 2D contour plot
            get_colour()
        except:
            print('Error! problem opening file %s' %file.name)
            Message(title='Error!',
                    message='Error! problem opening file %s' % file.name,
                    master=root).show()
            torsion_path.config(text='Error! no data selected from %s' % file.name)
            torsion_path.grid(row=4, sticky=tk.W, columnspan=8)
    else:
        Message(title='Error!', message='Error no torsion data selected', master=root).show()

def output_dir():
    '''
    Get output directory to save all the generated tables
    '''
    dir = askdirectory()
    global path_out
    path_out = str(glob.glob(dir))[2:-2] #Remove the list brackets and speechmarks
    out_path.config(text=path_out)
    out_path.grid(row = 6, sticky = tk.W, columnspan=8)

def re_order(df):
    '''
    Re-orders the dataframe and returns the columns with the corresponding number
    '''
    # Extract the number from filename, replace the col names with the numbers and reorder it
    conf_name = [int(re.findall('\d+', i)[0]) for i in df.columns.to_list()]
    df.columns  = conf_name
    df = df.reindex(sorted(df.columns), axis=1)
    return df

def reset_index(df):
    df = df.reset_index()
    return df

def format_df(df):
    dfT=df.T
    dfT.reset_index(inplace=True)
    dfT.pop('index')
    return dfT

def clear_old_labels():
    plot_save.grid_forget()  # Clears the old labal
    run_avr_out_tab3.grid_forget()
    save_xyz.grid_forget()

def fig_save_plot(fig, Title):
    '''
    Save the figures as .html in the output directory
    '''
    clear_old_labels()
    try:
        out_name = '{}.html'.format(Title)
        output_dir = os.path.join(path_out, out_name)
        fig.write_html(output_dir)
        run_avr_out.config(text='File saved in {}'.format(output_dir))
        run_avr_out.grid(row=9, sticky=tk.W, columnspan = 5)

        plot_save.config(text='File saved in {}'.format(output_dir))
        plot_save.grid(row=8, sticky=tk.W, column=2, columnspan=5)
    except:
        Message(title='Error!', message='Can not save file, check your output directory', master=root).show()
        run_avr_out.config(text='Error! No file was saved')
        run_avr_out.grid(row=9, sticky=tk.W, columnspan = 5)

        plot_save.config(text='File saved in {}'.format(output_dir))
        plot_save.grid(row=8, sticky=tk.W, column=2, columnspan=5)

def callback_CLT(selection):
    '''
    Selects the value from the list for CLT
    '''
    global atom_val_clt
    atom_val_clt = selection
    print(selection)

def get_atom_list_CLT(df):
    '''
    Create a list of atoms to be used for selecting which atom to view in histogram plots for CLT
    '''
    # Create dropdown list for selecting atom
    atom_list_clt = df.index.to_list() #List of atoms
    atom_sel_clt = tk.StringVar(tab2) #Tkinter select atom
    atom_drop_list_clt = ttk.OptionMenu(tab2, atom_sel_clt, 'Atom',*atom_list_clt, command= callback_CLT) #Dropdown list, needs a title or would use the first input as a title
    atom_drop_list_clt.grid(row=6, column=1)

def callback_boots(selection):
    '''
    Selects the value from the list for CLT
    '''
    global atom_val_boots
    atom_val_boots = selection
    print(selection)

def get_atom_list_boots(df):
    '''
    Create a list of atoms to be used for selecting which atom to view in histogram plots for CLT
    '''
    # Create dropdown list for selecting atom
    atom_list_boots = df.index.to_list() #List of atoms
    atom_sel_boots = tk.StringVar(tab2) #Tkinter select atom
    atom_drop_list_boots = ttk.OptionMenu(tab2, atom_sel_boots, 'Atom',*atom_list_boots, command= callback_boots) #Dropdown list, needs a title or would use the first input as a title
    atom_drop_list_boots.grid(row=8, column=1)

def running_avr_std_plot(df, Title, dtype):
    '''
    Plot the chemical shift for each conformational ensemble
    '''
    df['x'] = df.index + 1
    df_melt = df.melt(id_vars="x", value_vars=df.columns[:-1])
    df_melt.rename(columns={'x': 'Conformational Ensemble', 'value': dtype, df_melt.columns[1]: 'Legend'}, inplace=True)
    fig = px.line(df_melt, x='Conformational Ensemble', y=dtype, color='Legend', title=Title)
    show_plot_c = save_plot.get()
    save_plot_c = show_plot.get()
    if show_plot_c == 'No' and save_plot_c =='No':
        run_avr_out.config(text='Average and standard deviation calculations completed for {} files'.format(len(df)+1))
        run_avr_out.grid(row=10, sticky=tk.W, columnspan = 5)
    if show_plot_c == 'Yes':
        fig.show()
    if save_plot_c == 'Yes':
        fig_save_plot(fig, Title)
    if show_plot_c =='Yes' and save_plot_c =='Yes':
        fig.show()
        fig_save_plot(fig, Title)

def running_average(df):
    '''
    Calculate the running average and standard deviation, then plot it.
    '''
    df = re_order(df.copy())
    global df_avr, df_std_clt
    df_avr = df.expanding(axis=1).mean()
    df_std = df.expanding(axis=1).std()
    # Generate a list of the number of columns
    Num_of_Sample = df.T.reset_index()
    Num_of_Sample = Num_of_Sample.index.to_list()
    # First element is 1
    Num_of_Sample = [x + 1 for x in Num_of_Sample]
    # Sqrt as CLT takes sqrt N
    Num_of_Sample = [np.sqrt(x) for x in Num_of_Sample]
    #Std/sqrt(N)
    df_std_clt = df_std.div(Num_of_Sample)
    df_std_clt.dropna(axis=1, inplace=True)
    #Plot the averages and standard deviation
    df_avr_T = format_df(df_avr)
    df_std_clt_T = format_df(df_std_clt)
    running_avr_std_plot(df_avr_T, 'Running Average', 'Chemcial shift')
    running_avr_std_plot(df_std_clt_T, 'Running deviation', 'Std.dev/sqrt(N)')
    return df_std_clt

def select_tol(df, tolerance):
    '''
    View the number of samples within a set tolerance
    '''
    print('running convergence test')

    # Create Tkinter list to store output
    global Tol_list, Tol_scroll
    Tol_list = tk.Listbox(tab2, height=5)
    Tol_list.grid(row=13,  sticky='w', columnspan=5)
    Tol_scroll = tk.Scrollbar(tab2)
    Tol_scroll.grid(row=13, column=5, sticky='w')

    df = df * 1.96
    atom_list = df.index.to_list()
    for atom in atom_list:
        each_atom = df.loc[atom, :].round(2)
        within_tol = each_atom[each_atom < tolerance]
        if within_tol.size == 0:
            conv_status = 'Atom {} did not converge in the defined tolerance of {} '.format(atom, tolerance)
            Tol_list.insert(tk.END, conv_status)
        else:
            # Find the first convergance valie
            conv_val = within_tol.iloc[0]
            # Get the index value, which is the value of N
            Num_of_samp = within_tol[within_tol == conv_val]
            # Drop the duplicates corresponding to the convergance value
            Num_of_samp.drop_duplicates(inplace=True)
            conv_status = 'Atom {} converged after {} samples with a 95% CI of {}'.format(atom, Num_of_samp.index.to_list()[0],conv_val)
            Tol_list.insert(tk.END, conv_status)
    Tol_list.config(yscrollcommand=Tol_scroll.set, width=50)
    Tol_scroll.config(command=Tol_list.yview)


def get_tol_value_H(df):
    '''
    Get the tolerance for convergence testing
    '''
    df = running_average(df)
    df = df[df.index.str.contains('H')]
    print(df)
    Tol_shiftH = None
    try:
        Tol_shiftH = float(Tol_get_H.get())
    except:
        print('Error! Enter a valid number')
        Message(title='Error!', message='Enter a valid number', master=root).show()

    try:
        select_tol(df, Tol_shiftH)
    except:
        Message(title='Error!', message='Please run the running averages function and enter a valid tolerance', master=root).show()
        print('Please run the running averages and enter a valid tolerance')

def get_tol_value_C(df):
    '''
    Get the tolerance for convergence testing
    '''
    df = running_average(df)
    df = df[df.index.str.contains('C')]
    print(df)
    Tol_shiftC = None
    try:
        Tol_shiftC = float(Tol_get_C.get())
    except:
        print('Error! Enter a valid number')
        Message(title='Error!', message='Enter a valid number', master=root).show()
    try:
        select_tol(df, Tol_shiftC)
    except:
        Message(title='Error!', message='Please run the running averages function and enter a valid tolerance', master=root).show()
        print('Please run the running averages and enter a valid tolerance')


def plot_hist(df, atom, Mu,LCI, UCI, Title):
    fig = px.histogram(df, x=atom, labels={atom:'Chemical shift (ppm)'},
                       title='Atom {} from {}'.format(atom, Title),
                       color_discrete_sequence=['rgb(27,158,119)'],
                       template= 'simple_white')

    fig.add_vline(x=Mu, line_width=2, line_dash='dash', line_color='black', opacity=1)
    fig.add_vrect(x0=LCI, x1=UCI, fillcolor='red', opacity=0.15, line_width=0)
    fig.show()
    if save_plot_clt_hist.get() == 'Yes' or save_plot_boots_hist.get() =='Yes':
        fig_save_plot(fig, Title)

def plot_qq(df, atom, numSample):
    title = 'Normal Q-Q plot for atom {} N = {}'.format(atom, numSample)
    fig = go.Figure()
    qq = stats.probplot(df[atom], dist='norm')
    x = np.array([qq[0][0][0], qq[0][0][-1]])
    fig.add_scatter(x=qq[0][0], y=qq[0][1], mode='markers')
    fig.add_scatter(x=x, y=qq[1][1] + qq[1][0] * x, mode='lines')
    fig.layout.update(showlegend=False, title=title,
                      xaxis_title='Theoretical Quantiles', yaxis_title='Sample Quantiles')
    fig.show()
    if save_plot_clt_hist.get() == 'Yes' or save_plot_boots_hist.get() =='Yes':
        fig_save_plot(fig, title)


def get_N(df):
    '''
    Get the number of samples for CLT then run the CLT
    '''
    Num_sampl = int(Num_samp.get())
    max_val = int(len(df.columns))
    if Num_sampl > max_val:
        Num_sampl = max_val
        Message(title='Error!', message='Value greater than number of samples, setting N as {}'.format(max_val),
                master=root).show()

    try:
        CLT(df, Num_sampl, atom_val_clt)
    except:
        CLT(df, Num_sampl, None)

def get_M(df):
    '''
    Get the number of samples for Bootstrapping then run the Bootstrap function
    '''
    Num_sampl = int(M_samp.get())
    max_val = int(len(df.columns))
    if Num_sampl > max_val:
        Num_sampl = max_val
        Message(title='Error!', message='Value greater than number of samples, setting N as {}'.format(max_val), master=root).show()
    try:
        bootstrap(df, Num_sampl, atom_val_boots)
    except:
        bootstrap(df, Num_sampl, None)

def save_df(df, Title,tab):
    clear_old_labels()
    try:
        out_name = '{}.csv'.format(Title)
        output_dir = os.path.join(path_out, out_name)
        df.to_csv(output_dir)
        if tab == 'tab2':
            run_avr_out.config(text='File saved in {}'.format(output_dir))
            run_avr_out.grid(row=11, sticky=tk.W, columnspan = 5)
        if tab == 'tab3':
            run_avr_out_tab3.config(text='File saved in {}'.format(output_dir))
            run_avr_out_tab3.grid(row=8, sticky=tk.W, column=2, columnspan=5)
    except:
        Message(title='Error!', message='Can not save file, check your output directory', master=root).show()
        if tab == 'tab2':
            run_avr_out.config(text='Error! No file was saved')
            run_avr_out.grid(row=11, sticky=tk.W, columnspan = 5)
        if tab == 'tab3':
            run_avr_out_tab3.config(text='Error! No file was saved')
            run_avr_out_tab3.grid(row=8, sticky=tk.W, column=2, columnspan=5)


def CLT(df, numSample, atom):
    '''
    Central limit theorum
    '''
    df = re_order(df.copy())
    df = df.iloc[:, :numSample]
    s_mean = df.mean(axis=1)
    s_std = df.std(axis=1)
    df_clt = pd.DataFrame()
    df_clt['Average'] = s_mean
    df_clt['St.dev'] = s_std
    df_clt['95% CI'] = (s_std/np.sqrt(numSample))*1.96
    df_clt['mu - CI'] = s_mean - df_clt['95% CI']
    df_clt['mu + CI'] = s_mean + df_clt['95% CI']
    save_df(df_clt, 'CLT(N={})'.format(numSample), 'tab2')
    if atom !=None and plot_hist_c.get()=='Yes':
        LCI = df_clt.loc[atom, 'mu - CI']
        UCI = df_clt.loc[atom, 'mu + CI']
        Mu = df_clt.loc[atom, 'Average']
        plot_hist(df.T, atom, Mu, LCI ,UCI, 'CLT')
    print(df_clt)

def bootstrap(df, numSample, atom):
    sample_list = []
    #Take a subset of data of size N used for bootsrapping
    df = df.iloc[:, :numSample]
    # Generate random numbers without repetition at a range of 0- 1mil with a total length of numSample(M)
    randomstate = np.random.default_rng(0)
    ran_num = randomstate.choice(10**6, numSample, replace=False, shuffle=False)
    for sample, num in zip(range(1,numSample+1), ran_num):
        # Take M sampels from the data frameß
        sample = df.sample(n=numSample, replace=True, random_state=num, axis=1)
        meas_sample = sample.mean(axis=1)
        sample_list.append(meas_sample)
    means_df = pd.concat(sample_list, axis=1)
    # Relabel the column so that it starts from 1
    means_df.columns = [x+1 for x in means_df.columns]
    means_of_means = means_df.mean(axis=1)
    st_dev = means_df.std(axis=1)
    CI = 1.96 * st_dev
    LCI = means_of_means - CI
    UCI = means_of_means + CI
    boostrap_df = pd.DataFrame(data={'Means of means': means_of_means, 'St.dev': st_dev,
                                     '95% CI': CI, 'mu - CI': LCI, 'mu + CI': UCI})
    if save_plot_boots_hist.get()=='Yes':
        save_df(boostrap_df, 'Bootstrap(N={})'.format(numSample), 'tab2')

    if atom !=None and plot_boots_hist.get()=='Yes':
        LCI = boostrap_df.loc[atom, 'mu - CI']
        UCI = boostrap_df.loc[atom, 'mu + CI']
        Mu = boostrap_df.loc[atom, 'Means of means']
        plot_hist(means_df.T, atom, Mu, LCI ,UCI, 'Bootstrap means')
    if atom != None and plot_qq_f.get()=='Yes':
        '''
        Plot each mean in a QQ plot
        '''
        means_T = means_df.T
        plot_qq(means_T, atom, numSample)

def _Transpose_df(df):
    '''
    Transpose the df so that the file names are in the columns
    Extract the numbers from the filename and reorder the rows
    '''
    df_T = df.T
    df_T.reset_index(inplace=True)
    df_T['file_num'] = df_T['index'].str.extract('(\d+)').astype(int) + 1 # Need to add one to each file name as it starts from 0 instead of 1
    df_T.set_index('file_num', inplace=True)
    df_T.rename(columns={'index':'File'}, inplace=True)
    df_T.sort_index(inplace=True)
    return df_T

def getNumBins():
    global NumBins, bins, labels
    NumBins = int(NumBins_val.get())
    '''
    Enter the number of bins and print the size of the bin
    '''
    # Create a list of number of bins for the angles
    bins = []
    labels = []
    if NumBins !=0:
        width = int(360 / NumBins)
        print('Width of the bins are: {} degrees'.format(width))
        NumBins_label.config(text='Width of bins: {} degrees'.format(width))
        NumBins_label.grid(row=0, column=3, columnspan =4, sticky=tk.W)
        for i in range(NumBins+1):
            ub = int(-180) + width/2 +width*i
            lb = int(-180) - width/2 +width*i
            bins.append(lb)
            bins.append(ub)
            label = int(-180) + width*i
            labels.append(label)
        # Drop duplicates
        bins = list(dict.fromkeys(bins))

    else:
        print('Number of bins used: {}'.format(NumBins))
        NumBins_label.config(text='Number of bins used: {}'.format(NumBins))
        NumBins_label.grid(row=0, column=3, columnspan =4 , sticky=tk.W)

def bin_df(df, bins,labels, func):
    '''
    Groups the df into bins by the angles
    '''
    tmpdf =  df.groupby(pd.cut(df.index, bins=bins, labels=labels)).mean()
    if 180 or -180 in tmpdf.index:
        tmpdf2 = df.groupby(pd.cut(df.index, bins=bins, labels=labels)).count()  # Count number of rows
        tmpdf2_tot = tmpdf2.loc[-180] + tmpdf2.loc[180] #Combined count
        tmpdf.loc[180] = (tmpdf.loc[180] * tmpdf2.loc[180] + tmpdf.loc[-180] * tmpdf2.loc[
            -180]) / tmpdf2_tot  # Mean of -180 and 180, symmetrical about zero
        tmpdf.loc[-180] = tmpdf.loc[180]  # Set -180 to equal to combined mean
    if func == 'mean':
        df1 = tmpdf #tmpdf is mean calculation
    elif func == 'stdev':
        df1 = df.groupby(pd.cut(df.index, bins=bins, labels=labels)).std()
        print(df1.loc[[-180, 180]])
        if 180 or -180 in df1.index:
            mean_df = df.groupby(pd.cut(df.index, bins=bins, labels=labels)).mean() #calc mean again
            tmpdf = tmpdf.subtract(mean_df)
            tmpdf = tmpdf**2 # Square to remove negatives
            tmpdf = tmpdf.add(df1**2) #Add original variance NOT stdev
            tmpdf = tmpdf.mul(tmpdf2) # multiply by the count
            tot = tmpdf.loc[180] + tmpdf.loc[-180] # sume Ex2
            s2 = tot/tmpdf2_tot
            s = np.sqrt(s2) #Calculate combined sdev
            df1.loc[[180,-180]] = [s,s] #Set the sdev equal to each other
    elif func == 'min':
        df1 = df.groupby(pd.cut(df.index, bins=bins, labels=labels)).min()
        if 180 or -180 in df1.index:
            mindf = df1.loc[[180,-180]].min()
            df1.loc[[180, -180]] = [mindf, mindf]
    elif func == 'max':
        df1 = df.groupby(pd.cut(df.index, bins=bins, labels=labels)).max()
        if 180 or -180 in df1.index:
            maxdf = df1.loc[[180,-180]].max()
            df1.loc[[180, -180]] = [maxdf, maxdf]
    elif func =='range':
        maxdf = df.groupby(pd.cut(df.index, bins=bins, labels=labels)).max()
        mindf = df.groupby(pd.cut(df.index, bins=bins, labels=labels)).min()
        df1 = maxdf - mindf
        if 180 or -180 in mindf.index:
            mindf180 = mindf.loc[[180,-180]].min()
            maxdf180 = maxdf.loc[[180,-180]].max()
            range180 = maxdf180-mindf180
            df1.loc[[180, -180]] = [range180, range180]
    return df1

def Dof_list_energy_plot(df):
    '''
    DOF list for plot vs Chemical shift
    '''
    dof_list = df.columns.to_list()
    which_dof = tk.StringVar(tab3)
    which_dof_drop_list_cs = ttk.OptionMenu(tab3, which_dof, 'DOF',*dof_list, command= callback_DOF_cs) #Dropdown list, needs a title or would use the first input as a title
    which_dof_drop_list_cs.grid(row=1, column=1)

def callback_DOF_cs(selection):
    '''
    Selects the value from the list for CLT
    '''
    global DOF_cs
    DOF_cs = selection
    print(selection)

def Atom_table(df):
    '''
    Select atom for exporting table of stats from DOF energy plot
    '''
    atom_list = df.index.to_list()
    which_atom = tk.StringVar(tab3)
    which_atom_drop_list_cs = ttk.OptionMenu(tab3, which_atom, 'DOF',*atom_list, command= callback_atom_cs) #Dropdown list, needs a title or would use the first input as a title
    which_atom_drop_list_cs.grid(row=1, column=6)

def callback_atom_cs(selection):
    '''
    Selects the value from the list for DOF
    '''
    global Atom_cs
    Atom_cs = selection
    print(selection)

def DOF_energy_plot(DOF_df, df):
    try:
        NumBins
    except:
        Message(title='Error!', message='Please select a bin number', master=root).show()

    try:
        DOF_cs
    except:
        Message(title='Error!', message='Please select a DOF', master=root).show()
    df = _Transpose_df(df.copy())
    print(len(df))
    Each_Dof = DOF_df.loc[:,DOF_cs]
    print(len(Each_Dof))
    #Combine the DOF data with chemical shift, matches based on indices
    data = pd.concat([df,Each_Dof], join='inner', axis=1)
    data.sort_values(by=Each_Dof.name, ascending=True, inplace=True)
    data = data.set_index(Each_Dof.name)
    # Create a list of angles of each file for the DOF
    File_to_Ang_s = data.loc[:,'File']
    # Remove the file names to create plot with data
    data.pop('File')
    if plot_dof_cs.get() == 'Yes':
        '''
        Plot Dof vs shift
        '''
        if NumBins <= 1:
            data['x'] = data.index
            df_melt = data.melt(id_vars='x', value_vars=data.columns[:-1]).round(2)
            df_melt.rename(columns={'x': 'Angle', 'value': 'Chemical Shift', 'variable': 'Legend'}, inplace=True)
            min_data = df_melt.groupby(['Legend']).min()
            max_data = df_melt.groupby(['Legend']).max()
            range_data = round(max_data['Chemical Shift'] - min_data['Chemical Shift'], 2)
            df_melt['Range'] = df_melt['Legend'].map(range_data)
            fig = px.scatter(df_melt, x='Angle', y='Chemical Shift', color='Legend', title=Each_Dof.name,
                             hover_data=['Range'])
            fig.show()
        else:
            if plot_sdev_cs.get() == 'Yes':
                stat = 'stdev'
                title = 'Standard deviation'
            else:
                stat = 'range'
                title = 'Range'
            # Group angles into the bins and find the mean and std. dev
            data_binned_mean = bin_df(data, bins, labels, 'mean').round(2)
            data_binned_error = bin_df(data, bins, labels, stat).round(2)
            # Create columns to be used as the x axis in the plot
            data_binned_mean['x'] = data_binned_mean.index
            data_binned_error['x'] = data_binned_mean.index
            # Create dfs to be used by plotly
            bin_mean_melt = data_binned_mean.melt(id_vars='x', value_vars=data_binned_mean.columns[:-1])
            bin_mean_melt.rename(columns={'x': 'Angle','value': 'Mean Chemical Shift', 'variable': 'Legend'}, inplace=True)
            bin_error_melt = data_binned_error.melt(id_vars='x', value_vars=data_binned_error.columns[:-1])
            bin_error_melt.rename(columns={'value': title}, inplace=True)
            bin_plot = pd.concat([bin_mean_melt, bin_error_melt[title]], axis=1)
            print(bin_plot)
            min_data = bin_plot.groupby(['Legend']).min()
            max_data = bin_plot.groupby(['Legend']).max()
            range_data = round(max_data['Mean Chemical Shift'] - min_data['Mean Chemical Shift'], 2)
            bin_plot['Overall range'] = bin_plot['Legend'].map(range_data)
            bin_plot['yError'] = title
            fig = px.line(bin_plot, x='Angle', y='Mean Chemical Shift', color='Legend', error_y=title,
                          title=Each_Dof.name, hover_data=['Overall range', 'yError'])
            fig.show()
        if save_ang_cs.get() == 'Yes':
            fig_save_plot(fig, 'Chemical shift {} with {} bins'.format(DOF_cs, NumBins))

    elif export_csv_cs.get() == 'Yes':
        '''
        Export table as csv
        '''
        try:
            Atom_cs
        except:
            Message(title='Error!', message='Please select an Atom', master=root).show()
        if NumBins <= 1:
            Message(title='Error!', message='Please select more than 1 bin', master=root).show()
        elif NumBins > 1:
            df1 = pd.DataFrame()
            df1['Mean'] = bin_df(data.loc[:,Atom_cs], bins, labels, 'mean').round(2)
            df1['St Dev'] = bin_df(data.loc[:, Atom_cs], bins, labels, 'stdev').round(2)
            df1['Min'] = bin_df(data.loc[:, Atom_cs], bins, labels, 'min').round(2)
            df1['Max'] = bin_df(data.loc[:, Atom_cs], bins, labels, 'max').round(2)
            df1['Range'] = bin_df(data.loc[:, Atom_cs], bins, labels, 'range').round(2)
            df1.dropna(axis=0, inplace=True) #Drop blank rows
            print(df1)
            save_df(df1, '{} in {} Data'.format(Atom_cs, DOF_cs), 'tab3')
    else:
        Message(title='Error!', message='Please select an action', master=root).show()
    return File_to_Ang_s


def get_File(Angle, s):
    try:
        file = s.loc[Angle]
        print(file)
    except:
        print('Error! Angle %s not found in %s' %(Angle, s.index.name))

def Dof_X(df):
    '''
    DOF_X for contour plot
    '''
    dof_list = df.columns.to_list()
    which_dof = tk.StringVar(tab2)
    which_dof_drop_list_cs = ttk.OptionMenu(tab3, which_dof, 'DOF_X',*dof_list, command= callback_DOF_X) #Dropdown list, needs a title or would use the first input as a title
    which_dof_drop_list_cs.grid(row=2, column=1)

def callback_DOF_X(selection):
    '''
    DOF_X for contour plot
    '''
    global DOF_X
    DOF_X = selection
    print(selection)


def Dof_Y(df):
    '''
    DOF_Y for contour plot
    '''
    dof_list = df.columns.to_list()
    which_dof = tk.StringVar(tab2)
    which_dof_drop_list_cs = ttk.OptionMenu(tab3, which_dof, 'DOF_Y',*dof_list, command= callback_DOF_Y) #Dropdown list, needs a title or would use the first input as a title
    which_dof_drop_list_cs.grid(row=2, column=2)

def callback_DOF_Y(selection):
    '''
    DOF_Y for contour plot
    '''
    global DOF_Y
    DOF_Y = selection
    print(selection)

def DOF_Atom(df):
    '''
    Select atom for the pair of DOFs
    '''
    atom_list = df.index.to_list() #List of atoms
    atom_sel = tk.StringVar(tab3) #Tkinter select atom
    atom_drop_list = ttk.OptionMenu(tab3, atom_sel, 'Atom',*atom_list, command= callback_DOF_Atom)
    atom_drop_list.grid(row=2, column=3)

def callback_DOF_Atom(selection):
    '''
    Select atom for the pair of DOFs
    '''
    global Dof_atom
    Dof_atom = selection
    print(selection)

def get_colour():
    '''
    Choose a colour for the contour plots
    '''
    # Get a list of all colour schemes
    named_colorscales = sorted(px.colors.named_colorscales())
    global colour_2D
    col_list_2D = tk.StringVar(tab3) #Tkinter select atom
    col_drop_list = ttk.OptionMenu(tab3, col_list_2D, 'Colour',*named_colorscales, command=colour_select_2D)
    # Set as default
    col_list_2D.set('deep')
    col_drop_list.grid(row=2, column=4)

def colour_select_2D(selection):
    global colour_2D
    colour_2D = selection
    print(selection)

def Contour_plot(Dof_df, df):
    # Raise error if variables are not selected
    try: colour_2D
    except:
        Message(title='Error!', message='Please choose a colour', master=root).show()
    try: DOF_X
    except:
        Message(title='Error!', message='Please select DOF_X', master=root).show()
    try: DOF_Y
    except:
        Message(title='Error!', message='Please select DOF_Y', master=root).show()
    try: Dof_atom
    except:
        Message(title='Error!', message='Please select an atom', master=root).show()

    Bins_DOF_X = int(Dof_x_bin.get())
    Bins_DOF_Y = int(Dof_y_bin.get())

    Contour_bins_lab.config(text='DOF_X = {} DOF_Y={}'.format(Bins_DOF_X,Bins_DOF_Y))
    Contour_bins_lab.grid(row=3, column=5, columnspan=3, sticky=tk.W)

    print(Bins_DOF_X)
    print(Bins_DOF_Y)

    try: Bins_DOF_X
    except:
        Message(title='Error!', message='Please select a bin number for DOF_X', master=root).show()

    try: Bins_DOF_Y
    except:
        Message(title='Error!', message='Please select a bin number for DOF_Y', master=root).show()

    df = _Transpose_df(df.copy())
    Dof1 = Dof_df[DOF_X].astype(int)
    Dof2 = Dof_df[DOF_Y].astype(int)
    atom = df[Dof_atom]
    df1 = pd.concat([atom,Dof1, Dof2], join='inner', axis=1)
    tmpdf = pd.DataFrame(data={DOF_X: [180, -180], DOF_Y: [180, -180]})
    df2 = pd.concat([tmpdf, df1], axis=0)

    fig = px.density_contour(df2, x=DOF_X, y=DOF_Y, z=Dof_atom, histfunc='avg',
                             nbinsx=Bins_DOF_X, nbinsy=Bins_DOF_Y, title = Dof_atom + '(%s vs %s)' %(DOF_X,DOF_Y)
                             )
    fig.update_traces(contours_coloring='fill', colorscale = colour_2D, contours_showlines=False)
    fig.show()
    if save_ang_DofDof.get() =='Yes':
        fig_save_plot(fig, 'Contour Plot of Atom {}({} vs {})'.format(Dof_atom, DOF_X, DOF_Y))

def plot_3D(Dof_df, df):
    '''
    Plot a 3D scatter plot for DOF_X vs DOF_Y
    '''
    try:
        DOF_X
    except:
        Message(title='Error!', message='Please select DOF_X', master=root).show()
    try:
        DOF_Y
    except:
        Message(title='Error!', message='Please select DOF_Y', master=root).show()


    Dof1 = Dof_df[DOF_X].round(2)
    Dof2 = Dof_df[DOF_Y].round(2)
    df = _Transpose_df(df.copy())
    df.pop('File')
    df_data = pd.concat([df, Dof1, Dof2], join='inner', axis=1)

    tmpdf = pd.DataFrame(data={DOF_X:[180,-180], DOF_Y:[180,-180]})
    df2 = pd.concat([tmpdf,df_data], axis=0)

    # df_data = df_data.groupby(np.arange(len(df_data)) // NumBins).mean()
    df_melt = df2.melt(id_vars=[DOF_X,DOF_Y], value_vars=df_data.columns[:-2] )
    fig = px.scatter_3d(df_melt, x = DOF_X, y= DOF_Y, z = 'value', color='variable', labels={'value':'Chemical shift (ppm)' } )
    fig.show()

    if save_3D_Scatter.get() =='Yes':
        fig_save_plot(fig, '3D Scatter plot of {} vs {}'.format(DOF_X, DOF_Y))


def plot_Mesh(Dof_df, df):
    try: colour_2D
    except:
        Message(title='Error!', message='Please choose a colour', master=root).show()
    try: DOF_X
    except:
        Message(title='Error!', message='Please select DOF_X', master=root).show()
    try: DOF_Y
    except:
        Message(title='Error!', message='Please select DOF_Y', master=root).show()
    try: Dof_atom
    except:
        Message(title='Error!', message='Please select an atom', master=root).show()
    df = _Transpose_df(df.copy())
    Dof1 = Dof_df[DOF_X]
    Dof2 = Dof_df[DOF_Y]
    atom = df[Dof_atom]
    df.pop('File')
    df_data = pd.concat([atom, Dof1, Dof2], join='inner', axis=1).round(2)
    # df_data.sort_values(by=[DOF2, DOF1], ascending=[True, True], inplace=True)
    fig = go.Figure(data=[go.Mesh3d(x=df_data[DOF_X], y=df_data[DOF_Y], z=df_data[Dof_atom],
                                    opacity=1, intensity=df_data[Dof_atom], colorscale=colour_2D,
                                    hovertemplate=DOF_X + ': %{x}' + \
                                                  '<br>' + DOF_Y + ': %{y}' + \
                                                  '<br>' + Dof_atom + ': %{z}<extra></extra>'
                                    )],
                    layout=go.Layout(scene=dict(
                        xaxis_title=DOF_X,
                        yaxis_title=DOF_Y,
                        zaxis_title='Chemical shift (ppm)', ),
                        title=go.layout.Title(
                            text='Mesh plot of {} vs {} for atom {}'.format(DOF_X, DOF_Y, Dof_atom))
                    ))
    fig.show()
    if save_Mesh.get() == 'Yes':
        fig_save_plot(fig, 'Mesh plot of {} vs {} for atom {}'.format(DOF_X, DOF_Y, Dof_atom))

def get_limits():
    global upper_val, lower_val
    lower_val = int(lower_val_label.get())
    upper_val = int(upper_val_lab.get())

def plot_DOFvsDOF_static(Dof_df, df):
    try:
        get_limits()
    except:
        Message(title='Error!', message='Please enter an upper or lower angle', master=root).show()
    df = _Transpose_df(df.copy())
    Dof1 = Dof_df[DOF_X]
    Dof2 = Dof_df[DOF_Y]
    df1 = pd.concat([df, Dof1, Dof2], join='inner', axis=1)
    df1.pop('File')
    if upper_val <= lower_val:
        print('Error! No value exists from %s to %s'%(upper_val, lower_val))
        Message(title='Error!', message='Invalid entry of from %s to %s'%(upper_val, lower_val), master=root).show()
    else:
        df_sel = df1[(df1[DOF_X] < upper_val) & (df1[DOF_X]>lower_val)]
        print('{} data points in this range'.format(len(df_sel)))
        df_sel.pop(DOF_X)
        data = df_sel.set_index(DOF_Y)
        data.sort_index(inplace=True)
        data['x'] = data.index
        df_melt = data.melt(id_vars='x', value_vars=data.columns[:-1]).round(2)
        df_melt.rename(columns={'x': 'Angle', 'value': 'Chemical Shift', 'variable': 'Legend'}, inplace=True)
        min_data = df_melt.groupby(['Legend']).min()
        max_data = df_melt.groupby(['Legend']).max()
        range_data = round(max_data['Chemical Shift'] - min_data['Chemical Shift'],2)
        mean_data = df_melt.groupby(['Legend']).mean().round(2)
        std_data = df_melt.groupby(['Legend']).std().round(2)
        df_melt['Range'] = df_melt['Legend'].map(range_data)
        df_melt['Mean'] = df_melt['Legend'].map(mean_data['Chemical Shift'])
        df_melt['Std'] = df_melt['Legend'].map(std_data['Chemical Shift'])
        df_melt['Mean'] = df_melt['Mean'].astype(str) + ' ' + u'\u00B1' + ' ' +df_melt['Std'].astype(str)
        print(df_melt)
        fig = px.line(df_melt, x='Angle', y='Chemical Shift', color='Legend',
                      title='{} for angles {}<{}<{}'.format(DOF_Y,lower_val,DOF_X,upper_val), hover_data=['Range']+['Mean'])
        fig.show()
        if save_static.get() == 'Yes':
            fig_save_plot(fig, '{} for angles {}<{}<{}'.format(DOF_Y,lower_val,DOF_X,upper_val))

def Dof_angle(df):
    '''
    Select DOF for the file search when a range of angles are entered
    '''
    dof_list = df.columns.to_list()
    which_dof = tk.StringVar(tab3)
    which_dof_drop_list_cs = ttk.OptionMenu(tab3, which_dof, 'DOF',*dof_list, command= callback_Dof_angle)
    which_dof_drop_list_cs.grid(row=7, column=1)

def callback_Dof_angle(selection):
    '''
    DOF_Y for contour plot
    '''
    global DOF_angle
    DOF_angle = selection
    print(selection)

def get_angle(df, df2):
    '''
    Get the number of samples for CLT then run the CLT
    '''

    global angle_x, angle_y
    angle_x = None
    angle_y = None
    try:
        angle_x = float(Angle_x.get())
    except:
        print('Angle {} is not a valid number please enter a valid number'.format(Angle_x.get()))
        Message(title='Error!', message='Angle {} is not a valid number please enter a valid number'.format(Angle_x.get()), master=root).show()

    try:
        angle_y = float(Angle_y.get())
    except:
        print('Angle {} is not a valid number please enter a valid number'.format(Angle_y.get()))
        Message(title='Error!', message='Angle {} is not a valid number please enter a valid number'.format(Angle_x.get()), master=root).show()

    if angle_x > angle_y:
        print('Lower angle of {} can not be greater than the upper limit of {}'.format(angle_x,angle_y))
        Message(title='Error!', message='Lower angle of {} can not be greater than the upper limit of {}'.format(angle_x,angle_y), master=root).show()

    elif angle_x != None and angle_y != None:
        print(angle_x)
        print(angle_y)
        # try:
        get_angles(df, df2)
        # except:
        #     print('Please load chemical shift data and torsion data')
    else:
        print('No angles entered')
        Message(title='Error!', message='No angles entered', master=root).show()


def get_angles(df, df2):
    df2 = _Transpose_df(df2.copy())
    s = df[DOF_angle]
    s = s[s.between(angle_x, angle_y, inclusive=True)]
    df = pd.concat([s, df2], axis=1, join='inner')
    if len(df)<1:
        print('No angles found in this range')
    else:
        print(df)

    # Clear old tree view
    clear_tree()
    # Sort by angles
    df = df.sort_values(by=[DOF_angle])
    # make tree
    my_tree['column'] = df[[DOF_angle, 'File']].columns.to_list()
    my_tree['show'] = 'headings'

    # loop through column list
    for column in my_tree['column']:
        my_tree.heading(column, text=column)

    # Put data in treeview
    # convert to np array then to list
    df_rows = df.to_numpy().tolist()
    for row in df_rows:
        my_tree.insert('', 'end', values=row)
    my_tree.grid(row=9, columnspan=8)

    if save_DOF_cs.get() =='Yes':
        save_df(df, 'Angles {} to {} for'.format(angle_x,angle_y,DOF_angle), 'tab3')


def clear_tree():
    '''
    Clears old tree
    '''
    my_tree.delete(*my_tree.get_children())

def magres_to_xyz():
    clear_old_labels()
    file = askopenfile(initialfile= r'/Users/mkr97/', title='Select magres file', mode='r', filetypes=(('Magres files', '*.magres'),))
    file = file.name
    print(file)
    try:
        with open(file, 'r') as f:
            fname = os.path.basename(file)
            fname_noext = os.path.splitext(fname)[0]
            atoms = []
            for line in f:
                if ('atom' not in line):
                    continue
                ki = line.split()
                if (len(ki)<4):
                    continue
                atoms.append([ki[2],ki[3],ki[4],ki[5],ki[6]])
            with open(os.path.join(path_out, fname_noext + '.xyz'), 'w') as outfile:
                outfile.write('%s\n'%len(atoms))
                outfile.write('%s\n'%fname_noext)
                for i in atoms:
                    outfile.write('%s%s \t%f \t%f \t%f\n'%(i[0],i[1],float(i[2]),float(i[3]),float(i[4])))
                print('File %s.xyz saved to %s'%(fname_noext, path_out))
        save_xyz.config(text='File {}.xyz was saved in {}'.format(fname_noext, path_out))
        save_xyz.grid(row=8, sticky=tk.W, column=2,columnspan=5)
    except:
        Message(title='Error!', message='Can not save file, check your output directory', master=root).show()
        save_xyz.config(text='Error! No file was saved')
        save_xyz.grid(row=8, sticky=tk.W, column=2, columnspan=5)

def reset():
    print('Clear list')
    Tol_list.delete(0, tk.END)

def main():
    global root, tab2, tab3
    root = tk.Tk()
    style = ttk.Style()
    style.theme_use('classic')
    root.title('Solution Ensemble')
    root.geometry('800x400')
    tabControl = ttk.Notebook(root)

    tab1 = ttk.Frame(tabControl)
    tab2 = ttk.Frame(tabControl)
    tab3 = ttk.Frame(tabControl)


    tabControl.add(tab1, text='Load data')
    tabControl.add(tab2, text='Average')
    tabControl.add(tab3, text='Chem_vs_angles')
    tabControl.pack(expand=1, fill='both')


    '''
    Read input and get output directory
    '''

    ttk.Label(tab1, text='Chemical shift data').grid(row=1, sticky=tk.W)
    ttk.Label(tab1, text='Torsion data').grid(row=3, sticky=tk.W)
    ttk.Label(tab1, text='Output directory').grid(row=5, sticky=tk.W)

    global chem_shift_path, torsion_path, out_path, print_lab
    chem_shift_path = ttk.Label(tab1)
    torsion_path = ttk.Label(tab1)
    out_path = ttk.Label(tab1)
    print_lab = ttk.Label(tab1)

    ttk.Button(tab1, text='Open', command=lambda: read_df()).grid(row=1, column=1, sticky=tk.W)
    ttk.Button(tab1, text='Open', command=lambda: get_DOF_data()).grid(row=3, column=1, sticky=tk.W)
    ttk.Button(tab1, text='Open', command=lambda: output_dir()).grid(row=5, column=1, sticky=tk.W)

    '''
    Analysis:
    Running average and standard deviaiton
    '''
    global run_avr_out, save_plot, show_plot
    run_avr_out = ttk.Label(tab2)
    save_plot = tk.StringVar()
    show_plot = tk.StringVar()

    ttk.Label(tab2, text='Running average:').grid(row=1, sticky=tk.W)

    ttk.Button(tab2, text='Run', command=lambda: running_average(df_out)).grid(row=1, column=3, sticky=tk.W)
    # Show plots
    c1 = tk.Checkbutton(tab2, text='Plot graphs', variable=save_plot, onvalue='Yes', offvalue='No')
    c1.deselect()
    c1.grid(row=1, column=1, sticky=tk.W)
    # Save plots
    c2 = tk.Checkbutton(tab2, text='Save graphs', variable=show_plot, onvalue='Yes', offvalue='No')
    c2.deselect()
    c2.grid(row=1, column=2, sticky=tk.W)

    #Enter chemical shift and do CLT test for convergence
    #Tolerance for H
    global Tol_get_H, Tol_get_C, Tol_list,Tol_scroll
    ttk.Label(tab2, text='Tolerance (H):').grid(row=3, sticky=tk.W)
    Tol_get_H = ttk.Entry(tab2)
    ttk.Button(tab2, text='Run', command=lambda: get_tol_value_H(df_out)).grid(row=3, column=3, sticky=tk.W)
    Tol_get_H.grid(row=3, column =1, columnspan=2)
    # Tolerance for C
    ttk.Label(tab2, text='Tolerance (C):').grid(row=4, sticky=tk.W)
    Tol_get_C = ttk.Entry(tab2)
    ttk.Button(tab2, text='Run', command=lambda: get_tol_value_C(df_out)).grid(row=4, column=3, sticky=tk.W)
    Tol_get_C.grid(row=4, column =1, columnspan=2)

    # Tol_list = tk.Listbox(tab2)
    # Tol_scroll = tk.Scrollbar(tab2)

    # Perform CLT and save df
    global Num_samp, plot_hist_c, save_plot_clt_hist
    ttk.Label(tab2, text='Num of sample (N):').grid(row=5, sticky=tk.W)
    Num_samp = ttk.Entry(tab2)
    ttk.Button(tab2, text='Run', command=lambda: get_N(df_out)).grid(row=6, column=3, sticky=tk.W)
    Num_samp.grid(row=5, column =1, columnspan=2)

    ttk.Label(tab2, text='Histogram(CLT):').grid(row=6, sticky=tk.W)

    #Checkbox to view plot
    plot_hist_c = tk.StringVar()
    plot_hist_cb = tk.Checkbutton(tab2, text='Plot hist', variable=plot_hist_c, onvalue='Yes', offvalue='No')
    plot_hist_cb.deselect()
    plot_hist_cb.grid(row=5, column=3, sticky=tk.W)

    #Checkbox to save plot
    save_plot_clt_hist = tk.StringVar()
    save_plot_clt_hist_c = tk.Checkbutton(tab2, text='Save', variable=save_plot_clt_hist, onvalue='Yes', offvalue='No')
    save_plot_clt_hist_c.deselect()
    save_plot_clt_hist_c.grid(row=6, column=2, sticky=tk.W)

    # Bootstrap
    global M_samp, plot_boots_hist, save_plot_boots_hist, plot_qq_f
    ttk.Label(tab2, text='Num of sample (N):').grid(row=7, sticky=tk.W)
    M_samp = ttk.Entry(tab2)
    ttk.Button(tab2, text='Run', command=lambda: get_M(df_out)).grid(row=8, column=3, sticky=tk.W)
    M_samp.grid(row=7, column =1, columnspan=2)
    ttk.Label(tab2, text='Histogram(Bootstrap):').grid(row=8, sticky=tk.W)

    #Checkbox to view plot
    plot_boots_hist = tk.StringVar()
    plot_boots_hist_c = tk.Checkbutton(tab2, text='Plot hist', variable=plot_boots_hist, onvalue='Yes', offvalue='No')
    plot_boots_hist_c.deselect()
    plot_boots_hist_c.grid(row=7, column=3, sticky=tk.W)

    plot_qq_f = tk.StringVar()
    plot_qq_c = tk.Checkbutton(tab2, text='Plot QQ', variable=plot_qq_f, onvalue='Yes', offvalue='No')
    plot_qq_c.deselect()
    plot_qq_c.grid(row=7, column=4, sticky=tk.W)

    #Checkbox to save plot
    save_plot_boots_hist = tk.StringVar()
    save_plot_boots_hist_c = tk.Checkbutton(tab2, text='Save', variable=save_plot_boots_hist, onvalue='Yes', offvalue='No')
    save_plot_boots_hist_c.deselect()
    save_plot_boots_hist_c.grid(row=8, column=2, sticky=tk.W)


    #Reset
    ttk.Button(tab2, text='Reset', command=lambda: reset()).grid(row=9, sticky=tk.W)
    '''
    Plotting
    Angle vs chemical shift
    2D DOF vs DOF
    '''

    global NumBins_val, NumBins_label
    #Number of bins create entry box
    ttk.Label(tab3, text='Select number of bins:').grid(row=0, sticky=tk.W, columnspan=1)

    #Shows the entered number of bins
    NumBins_label = ttk.Label(tab3)
    NumBins_val = ttk.Entry(tab3, width=5)
    NumBins_val.insert(0, 0)
    ttk.Button(tab3, text='Run', command=lambda: getNumBins()).grid(row=0, column=2, sticky=tk.W)
    NumBins_val.bind('<Return>', (lambda event: getNumBins()))
    NumBins_val.grid(row=0, column=1)


    #Chem shift vs angles
    ttk.Label(tab3, text='Plot torsion with chemical shift:').grid(row=1, sticky=tk.W, columnspan=1)
    #Select DOF

    # Save plot
    global save_ang_cs, plot_save
    plot_save = ttk.Label(tab3)

    save_ang_cs = tk.StringVar()
    save_ang_cs_c = tk.Checkbutton(tab3, text='Save', variable=save_ang_cs, onvalue='Yes', offvalue='No')
    save_ang_cs_c.deselect()
    save_ang_cs_c.grid(row=1, column=2, sticky=tk.W)
    ttk.Button(tab3, text='Run', command=lambda: DOF_energy_plot(DOF_df, df_out)).grid(row=1, column=7, sticky=tk.W)

    # Plot range or stdev checkbox
    # Export table checkbox
    global plot_sdev_cs, export_csv_cs, plot_dof_cs

    plot_sdev_cs = tk.StringVar()
    export_csv_cs = tk.StringVar()
    plot_dof_cs = tk.StringVar()

    plot_sdev_cs_c= tk.Checkbutton(tab3, text='stdev', variable=plot_sdev_cs, onvalue='Yes', offvalue='No')
    export_csv_cs_c = tk.Checkbutton(tab3, text='Table', variable=export_csv_cs, onvalue='Yes', offvalue='No')
    plot_dof_cs_c = tk.Checkbutton(tab3, text='Plot', variable=plot_dof_cs, onvalue='Yes', offvalue='No')

    plot_sdev_cs_c.deselect()
    export_csv_cs_c.deselect()
    plot_dof_cs_c.deselect()

    plot_sdev_cs_c.grid(row=1, column=3, sticky=tk.W, padx=(0, 10))
    export_csv_cs_c.grid(row=1, column=5, sticky=tk.W)
    plot_dof_cs_c.grid(row=1, column=4, sticky=tk.W)

    #Creat dropdown list for table


    # 2D DOF Plots
    ttk.Label(tab3, text='DOF_X vs DOF_Y:').grid(row=2, sticky=tk.W)
    global save_ang_DofDof
    save_ang_DofDof = tk.StringVar()
    save_ang_cs_c = tk.Checkbutton(tab3, text='Save', variable=save_ang_DofDof, onvalue='Yes', offvalue='No')
    save_ang_cs_c.deselect()
    save_ang_cs_c.grid(row=3, column=3, sticky=tk.W)

    #Select bins for each DOF
    global Dof_x_bin, Dof_y_bin, Contour_bins_lab
    ttk.Label(tab3, text='Bins for each DOF:').grid(row=3, sticky=tk.W)
    Dof_x_bin = ttk.Entry(tab3, width=5)
    Dof_x_bin.insert(0, 0)
    Dof_y_bin = ttk.Entry(tab3, width=5)
    Dof_y_bin.insert(0, 0)

    ttk.Button(tab3, text='Run', command=lambda: Contour_plot(DOF_df, df_out)).grid(row=3, column=4, sticky=tk.W)

    #Prints the number of bins used for contour plot
    Contour_bins_lab = ttk.Label(tab3)

    Dof_x_bin.grid(row=3, column=1)
    Dof_y_bin.grid(row=3, column=2)

    # 3D Scatter plot
    global save_3D_Scatter
    ttk.Label(tab3, text='3D Scatter:').grid(row=4, sticky=tk.W)
    save_3D_Scatter = tk.StringVar()
    save_3D_Scatter_c = tk.Checkbutton(tab3, text='Save', variable=save_3D_Scatter, onvalue='Yes', offvalue='No')
    save_3D_Scatter_c.deselect()
    save_3D_Scatter_c.grid(row=4, column=1, sticky=tk.W, padx=(10,0))
    ttk.Button(tab3, text='Run', command=lambda: plot_3D(DOF_df, df_out)).grid(row=4, column=2, sticky=tk.W)

    # Mesh plot
    global save_Mesh
    ttk.Label(tab3, text='Mesh plot:').grid(row=5, sticky=tk.W)
    save_Mesh = tk.StringVar()
    save_Mesh_c = tk.Checkbutton(tab3, text='Save', variable=save_Mesh, onvalue='Yes', offvalue='No')
    save_Mesh_c.deselect()
    save_Mesh_c.grid(row=5, column=1, sticky=tk.W, padx=(10,0))
    ttk.Button(tab3, text='Run', command=lambda: plot_Mesh(DOF_df, df_out)).grid(row=5, column=2, sticky=tk.W)

    # Static plot
    global upper_val_lab, lower_val_label, save_static
    ttk.Label(tab3, text='DOF_X range:').grid(row=6, sticky=tk.W)

    #Select lower limit
    lower_val_label = ttk.Entry(tab3, width=5)
    lower_val_label.insert(0, 0)
    lower_val_label.grid(row=6, column=1)

    #Select upper limit
    upper_val_lab = ttk.Entry(tab3, width=5)
    upper_val_lab.insert(0, 0)
    upper_val_lab.grid(row=6, column=2)

    #Save
    save_static = tk.StringVar()
    save_static_c = tk.Checkbutton(tab3, text='Save', variable=save_static, onvalue='Yes', offvalue='No')
    save_static_c.deselect()
    save_static_c.grid(row=6, column=3, sticky=tk.W)
    ttk.Button(tab3, text='Run', command=lambda: plot_DOFvsDOF_static(DOF_df, df_out)).grid(row=6, column=4, sticky=tk.W)

    # Find the files for a range of angles
    # Angle X is smaller than Angle Y
    global Angle_x, Angle_y
    ttk.Label(tab3, text='Search files:').grid(row=7, sticky=tk.W)
    Angle_x = ttk.Entry(tab3, width=5)
    Angle_y = ttk.Entry(tab3, width=5)
    ttk.Label(tab3, text='Angle_LA:').grid(row=7, column=2, sticky=tk.W)
    ttk.Label(tab3, text='Angle_UA:').grid(row=7, column=4, sticky=tk.W)

    ttk.Button(tab3, text='Run', command=lambda: get_angle(DOF_df, df_out)).grid(row=7, column=7, sticky=tk.W)

    Angle_x.grid(row=7, column=3)
    Angle_y.grid(row=7, column=5)

    #Create treeview for df display
    global my_tree, save_DOF_cs, run_avr_out_tab3
    my_tree = ttk.Treeview(tab3, height=5)
    #Save df for tab3
    run_avr_out_tab3 = ttk.Label(tab3)
    # Save the df for chemcial shift of each angle
    save_DOF_cs = tk.StringVar()
    save_DOF_cs_c = tk.Checkbutton(tab3, text='Save', variable=save_DOF_cs, onvalue='Yes', offvalue='No')
    save_DOF_cs_c.deselect()
    save_DOF_cs_c.grid(row=7, column=6, sticky=tk.W)

    #Conver magres to xyz?
    global save_xyz
    save_xyz = ttk.Label(tab3)

    ttk.Label(tab3, text='Convert magres to xyz:').grid(row=8, columnspan=1, sticky=tk.W)
    ttk.Button(tab3, text='Open', command=lambda:magres_to_xyz()).grid(row=8, column=1, sticky=tk.W)

    root.mainloop()

main()

