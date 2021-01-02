import numpy as np
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
from pert import PERT
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

        
root = tk.Tk()

root.title('enerSYS')
root.iconbitmap('Logo_apk2.ico')

def no_zero(array):
    rbaru = []
    for i in array:
        if i < 0:
            rbaru.append(0)
        else:
            rbaru.append(i)
    return np.array(rbaru)

my_notebook = ttk.Notebook(root)
my_notebook.pack()

my_frame1 = tk.Frame(my_notebook)
my_frame2 = tk.Frame(my_notebook)

my_frame1.pack(fill="both", expand=1)
my_frame2.pack(fill="both", expand=1)

my_notebook.add(my_frame1, text="enerVOLTION")
my_notebook.add(my_frame2, text="enerECOTION")

clicked_tab = tk.StringVar()
clicked_tab.set("Choose Parameters")

###Label & Entry#####

minimum = tk.Label(my_frame1, text="Low")
minimum.grid(row=0, column=1)

most_likely = tk.Label(my_frame1, text="Most-Likely")
most_likely.grid(row=0, column=2)

maximum = tk.Label(my_frame1, text="High")
maximum.grid(row=0, column=3)

maximum = tk.Label(my_frame1, text="Distribution Type")
maximum.grid(row=0, column=4)

#--Reservoir Area--#
def simulation_reservoir_area():
    dist = clicked_reservoir_area.get()
    low = float(reservoir_area_e_min.get())
    most = float(reservoir_area_e_most.get())
    high = float(reservoir_area_e_max.get())
    trial = 200000

    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
        
    elif(dist == "Triangular"):
        r = np.random.triangular(low, most, high, trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')        
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')
    
    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

    return {
        "res": r,
        "val_min": round(min(r), 2),
        "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }

reservoir_area_l = tk.Label(my_frame1, text="Reservoir Area [km^2]")
reservoir_area_l.grid(row=1, column=0)

reservoir_area_e_min = tk.Entry(my_frame1, width=5)
reservoir_area_e_most = tk.Entry(my_frame1, width=5)
reservoir_area_e_max = tk.Entry(my_frame1, width=5)

reservoir_area_e_min.grid(row=1, column=1)
reservoir_area_e_most.grid(row=1, column=2)
reservoir_area_e_max.grid(row=1, column=3)

clicked_reservoir_area = tk.StringVar()
clicked_reservoir_area.set("Constant")

drop_reservoir_area = tk.OptionMenu(my_frame1, clicked_reservoir_area, "Constant", "Normal", "Triangular", "Beta-PERT", "Log Normal")
drop_reservoir_area.grid(row=1, column=4)

button_reservoir_area = tk.Button(my_frame1, text="Generate", command=simulation_reservoir_area)
button_reservoir_area.grid(row=1, column=5)

#--Reservoir Thickness--#
def simulation_reservoir_thickness():
    dist = clicked_reservoir_thickness.get()
    low = float(reservoir_thickness_e_min.get())
    most = float(reservoir_thickness_e_most.get())
    high = float(reservoir_thickness_e_max.get())
    trial = 200000

    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
        
    elif(dist == "Triangular"):

        r = np.random.triangular(low, most, high, trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')        
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')

    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

    return {
        "res": r,
        "val_min": round(min(r), 2),
        "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }

reservoir_thickness_l = tk.Label(my_frame1, text="Reservoir Thickness [m]")
reservoir_thickness_l.grid(row=2, column=0)

reservoir_thickness_e_min = tk.Entry(my_frame1, width=5)
reservoir_thickness_e_most = tk.Entry(my_frame1, width=5)
reservoir_thickness_e_max = tk.Entry(my_frame1, width=5)

reservoir_thickness_e_min.grid(row=2, column=1)
reservoir_thickness_e_most.grid(row=2, column=2)
reservoir_thickness_e_max.grid(row=2, column=3)

clicked_reservoir_thickness = tk.StringVar()
clicked_reservoir_thickness.set("Constant")

drop_reservoir_thickness = tk.OptionMenu(my_frame1, clicked_reservoir_thickness, "Constant", "Normal", "Triangular", "Beta-PERT", "Log Normal")
drop_reservoir_thickness.grid(row=2, column=4)

button_reservoir_thickness = tk.Button(my_frame1, text="Generate", command=simulation_reservoir_thickness)
button_reservoir_thickness.grid(row=2, column=5)

#--Rock Density--#
def simulation_rock_density():
    dist = clicked_rock_density.get()
    low = float(rock_density_e_min.get())
    most = float(rock_density_e_most.get())
    high = float(rock_density_e_max.get())
    trial = 200000

    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
        
    elif(dist == "Triangular"):

        r = np.random.triangular(low, most, high, trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')        
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')

    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

    return {
        "res": r,
        "val_min": round(min(r), 2),
        "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }

rock_density_l = tk.Label(my_frame1, text="Rock Density [kg/m^3]")
rock_density_l.grid(row=3, column=0)

rock_density_e_min = tk.Entry(my_frame1, width=5)
rock_density_e_most = tk.Entry(my_frame1, width=5)
rock_density_e_max = tk.Entry(my_frame1, width=5)

rock_density_e_min.grid(row=3, column=1)
rock_density_e_most.grid(row=3, column=2)
rock_density_e_max.grid(row=3, column=3)

clicked_rock_density = tk.StringVar()
clicked_rock_density.set("Constant")

drop_rock_density = tk.OptionMenu(my_frame1, clicked_rock_density, "Constant", "Normal", "Triangular", "BETA-Pert", "Log Normal")
drop_rock_density.grid(row=3, column=4)

button_reservoir_thickness = tk.Button(my_frame1, text="Generate", command=simulation_rock_density)
button_reservoir_thickness.grid(row=3, column=5)

#--Porosity--#
def simulation_porosity():
    dist = clicked_porosity.get()
    low = float(porosity_e_min.get())
    most = float(porosity_e_most.get())
    high = float(porosity_e_max.get())
    trial = 200000

    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
        
    elif(dist == "Triangular"):

        r = np.random.triangular(low, most, high, trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')        
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')

    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

    return {
        "res": r,
        "val_min": round(min(r), 2),
        "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }
porosity_l = tk.Label(my_frame1, text="Porosity [%]")
porosity_l.grid(row=4, column=0)

porosity_e_min = tk.Entry(my_frame1, width=5)
porosity_e_most = tk.Entry(my_frame1, width=5)
porosity_e_max = tk.Entry(my_frame1, width=5)

porosity_e_min.grid(row=4, column=1)
porosity_e_most.grid(row=4, column=2)
porosity_e_max.grid(row=4, column=3)

clicked_porosity = tk.StringVar()
clicked_porosity.set("Constant")

drop_porosity = tk.OptionMenu(my_frame1, clicked_porosity, "Constant", "Normal", "Triangular", "BETA-Pert", "Log Normal")
drop_porosity.grid(row=4, column=4)

button_porosity = tk.Button(my_frame1, text="Generate", command=simulation_porosity)
button_porosity.grid(row=4, column=5)

#--Recovery Factor--#
def simulation_recovery_factor():
    dist = clicked_recovery_factor.get()
    low = float(recovery_factor_e_min.get())
    most = float(recovery_factor_e_most.get())
    high = float(recovery_factor_e_max.get())
    trial = 200000

    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
        
    elif(dist == "Triangular"):

        r = np.random.triangular(low, most, high, trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')        
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')

    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

    return {
        "res": r,
        "val_min": round(min(r), 2),
        "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }
recovery_factor_l = tk.Label(my_frame1, text="Recovery Factor [%]")
recovery_factor_l.grid(row=5, column=0)

recovery_factor_e_min = tk.Entry(my_frame1, width=5)
recovery_factor_e_most = tk.Entry(my_frame1, width=5)
recovery_factor_e_max = tk.Entry(my_frame1, width=5)

recovery_factor_e_min.grid(row=5, column=1)
recovery_factor_e_most.grid(row=5, column=2)
recovery_factor_e_max.grid(row=5, column=3)

clicked_recovery_factor = tk.StringVar()
clicked_recovery_factor.set("Constant")

drop_recovery_factor = tk.OptionMenu(my_frame1, clicked_recovery_factor, "Constant", "Normal", "Triangular", "BETA-Pert", "Log Normal")
drop_recovery_factor.grid(row=5, column=4)

button_recovery_factor = tk.Button(my_frame1, text="Generate", command=simulation_recovery_factor)
button_recovery_factor.grid(row=5, column=5)

#--Rock Specific_Heat--#
def simulation_rock_spesific_heat():
    dist = clicked_rock_spesific_heat.get()
    low = float(rock_specific_heat_e_min.get())
    most = float(rock_specific_heat_e_most.get())
    high = float(rock_specific_heat_e_max.get())
    trial = 200000

    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
        
    elif(dist == "Triangular"):

        r = np.random.triangular(low, most, high, trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')        
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')

    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

    return {
        "res": r,
        "val_min": round(min(r), 2),
        "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }
rock_specific_heat_l = tk.Label(my_frame1, text="Rock Spesific Heat [kJ/kgC]")
rock_specific_heat_l.grid(row=6, column=0)

rock_specific_heat_e_min = tk.Entry(my_frame1, width=5)
rock_specific_heat_e_most = tk.Entry(my_frame1, width=5)
rock_specific_heat_e_max = tk.Entry(my_frame1, width=5)

rock_specific_heat_e_min.grid(row=6, column=1)
rock_specific_heat_e_most.grid(row=6, column=2)
rock_specific_heat_e_max.grid(row=6, column=3)

clicked_rock_spesific_heat = tk.StringVar()
clicked_rock_spesific_heat.set("Constant")

drop_rock_spesific_heat = tk.OptionMenu(my_frame1, clicked_rock_spesific_heat, "Constant", "Normal", "Triangular", "Beta-PERT", "Log Normal")
drop_rock_spesific_heat.grid(row=6, column=4)

button_rock_spesific_heat = tk.Button(my_frame1, text="Generate", command=simulation_rock_spesific_heat)
button_rock_spesific_heat.grid(row=6, column=5)

#--Reservoir Average Temperatur--#
def simulation_reservoir_temp():
    dist = clicked_reservoir_temp.get()
    low = float(reservoir_temp_e_min.get())
    most = float(reservoir_temp_e_most.get())
    high = float(reservoir_temp_e_max.get())
    trial = 200000

    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
        
    elif(dist == "Triangular"):

        r = np.random.triangular(low, most, high, trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')        
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')

    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

    return {
        "res": r,
        "val_min": round(min(r), 2),
        "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }
reservoir_temp_l = tk.Label(my_frame1, text="Reservoir Average Temperature [C]")
reservoir_temp_l.grid(row=7, column=0)

reservoir_temp_e_min = tk.Entry(my_frame1, width=5)
reservoir_temp_e_most = tk.Entry(my_frame1, width=5)
reservoir_temp_e_max = tk.Entry(my_frame1, width=5)

reservoir_temp_e_min.grid(row=7, column=1)
reservoir_temp_e_most.grid(row=7, column=2)
reservoir_temp_e_max.grid(row=7, column=3)

clicked_reservoir_temp = tk.StringVar()
clicked_reservoir_temp.set("Constant")

drop_reservoir_temp = tk.OptionMenu(my_frame1, clicked_reservoir_temp, "Constant", "Normal", "Triangular", "Beta-PERT", "Log Normal")
drop_reservoir_temp.grid(row=7, column=4)

button_reservoir_temp = tk.Button(my_frame1, text="Generate", command=simulation_reservoir_temp)
button_reservoir_temp.grid(row=7, column=5)

#--Heat-Electricity Conversion Efficiency--#
def simulation_heat_conversion():
    dist = clicked_heat_conversion.get()
    low = float(heat_conversion_e_min.get())
    most = float(heat_conversion_e_most.get())
    high = float(heat_conversion_e_max.get())
    trial = 200000

    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
        
    elif(dist == "Triangular"):

        r = np.random.triangular(low, most, high, trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')        
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')

    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

    return {
        "res": r,
        "val_min": round(min(r), 2),
        "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }
heat_conversion_l = tk.Label(my_frame1, text="Heat-Electricity Conversion Efficiency [C]")
heat_conversion_l.grid(row=8, column=0)

heat_conversion_e_min = tk.Entry(my_frame1, width=5)
heat_conversion_e_most = tk.Entry(my_frame1, width=5)
heat_conversion_e_max = tk.Entry(my_frame1, width=5)

heat_conversion_e_min.grid(row=8, column=1)
heat_conversion_e_most.grid(row=8, column=2)
heat_conversion_e_max.grid(row=8, column=3)

clicked_heat_conversion = tk.StringVar()
clicked_heat_conversion.set("Constant")

drop_heat_conversion = tk.OptionMenu(my_frame1, clicked_heat_conversion, "Constant", "Normal", "Triangular", "Beta-PERT", "Log Normal")
drop_heat_conversion.grid(row=8, column=4)

button_heat_conversion = tk.Button(my_frame1, text="Generate", command=simulation_heat_conversion)
button_heat_conversion.grid(row=8, column=5)

#--Plant Life--#
def simulation_plant_life():
    dist = clicked_plant_life.get()
    low = float(plant_life_e_min.get())
    most = float(plant_life_e_most.get())
    high = float(plant_life_e_max.get())
    trial = 200000

    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
        
    elif(dist == "Triangular"):

        r = np.random.triangular(low, most, high, trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')        
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')

    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

    return {
        "res": r,
        "val_min": round(min(r), 2),
        "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }
plant_life_l = tk.Label(my_frame1, text="Plant Life [Year]")
plant_life_l.grid(row=9, column=0)

plant_life_e_min = tk.Entry(my_frame1, width=5)
plant_life_e_most = tk.Entry(my_frame1, width=5)
plant_life_e_max = tk.Entry(my_frame1, width=5)

plant_life_e_min.grid(row=9, column=1)
plant_life_e_most.grid(row=9, column=2)
plant_life_e_max.grid(row=9, column=3)

clicked_plant_life = tk.StringVar()
clicked_plant_life.set("Constant")

drop_plant_life = tk.OptionMenu(my_frame1, clicked_plant_life, "Constant", "Normal", "Triangular", "Beta-PERT", "Log Normal")
drop_plant_life.grid(row=9, column=4)

button_plant_life = tk.Button(my_frame1, text="Generate", command=simulation_plant_life)
button_plant_life.grid(row=9, column=5)

#--Load Factor--#
def simulation_load_factor():
    dist = clicked_load_factor.get()
    low = float(load_factor_e_min.get())
    most = float(load_factor_e_most.get())
    high = float(load_factor_e_max.get())
    trial = 200000

    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
        
    elif(dist == "Triangular"):

        r = np.random.triangular(low, most, high, trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')        
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')

    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

    return {
        "res": r,
        "val_min": round(min(r), 2),
        "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }
load_factor_l = tk.Label(my_frame1, text="Load Factor [%]")
load_factor_l.grid(row=10, column=0)

load_factor_e_min = tk.Entry(my_frame1, width=5)
load_factor_e_most = tk.Entry(my_frame1, width=5)
load_factor_e_max = tk.Entry(my_frame1, width=5)

load_factor_e_min.grid(row=10, column=1)
load_factor_e_most.grid(row=10, column=2)
load_factor_e_max.grid(row=10, column=3)

clicked_load_factor = tk.StringVar()
clicked_load_factor.set("Constant")

drop_load_factor = tk.OptionMenu(my_frame1, clicked_load_factor, "Constant", "Normal", "Triangular", "Beta-PERT", "Log Normal")
drop_load_factor.grid(row=10, column=4)

button_load_factor = tk.Button(my_frame1, text="Generate", command=simulation_load_factor)
button_load_factor.grid(row=10, column=5)

#--Abandonment Temperature--#
def simulation_abandonment_temp():
    dist = clicked_abandonment_temp.get()
    low = float(abandonment_temp_e_min.get())
    most = float(abandonment_temp_e_most.get())
    high = float(abandonment_temp_e_max.get())
    trial = 200000

    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
        
    elif(dist == "Triangular"):

        r = np.random.triangular(low, most, high, trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        fig, ax1 = plt.subplots()
        count, bins, ignored = ax1.hist(r, 50, density=True, rwidth=0.9, color='teal')        
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        data_sorted = np.sort(r)
        p = 1. * np.arange(len(r)) / (len(r) - 1)
        ax2.plot(data_sorted, p, color='orangered', linewidth=5)
        ax2.tick_params(axis='y')

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')

    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

    return {
        "res": r,
        "val_min": round(min(r), 2),
        "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }
abandonment_temp_l = tk.Label(my_frame1, text="Abandonment Temperature [C]")
abandonment_temp_l.grid(row=11, column=0)

abandonment_temp_e_min = tk.Entry(my_frame1, width=5)
abandonment_temp_e_most = tk.Entry(my_frame1, width=5)
abandonment_temp_e_max = tk.Entry(my_frame1, width=5)

abandonment_temp_e_min.grid(row=11, column=1)
abandonment_temp_e_most.grid(row=11, column=2)
abandonment_temp_e_max.grid(row=11, column=3)

clicked_abandonment_temp = tk.StringVar()
clicked_abandonment_temp.set("Constant")

drop_abandonment_temp = tk.OptionMenu(my_frame1, clicked_abandonment_temp, "Constant", "Normal", "Triangular", "Beta-PERT", "Log Normal")
drop_abandonment_temp.grid(row=11, column=4)

button_abandonment_temp = tk.Button(my_frame1, text="Generate", command=simulation_abandonment_temp)
button_abandonment_temp.grid(row=11, column=5)

#--Volumetric Calculation--#

def simulation(low, most, high, dist = ""):
    trial = 200000
    if(dist == "Normal"):
        mu = (low+most+high)/3
        sigma = (high-low)/6
        r = np.random.normal(mu, sigma, trial)
        if min(r) < 0:
                r = no_zero(r)
        count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')
        
    elif(dist == "Triangular"):
        r = np.random.triangular(low, most, high, trial)
        count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')
    
    elif(dist == "Beta-PERT"):
        mu = (low + 4*most + high)/6
        sigma = (high-low)/6
        pert = PERT(low, most, high)
        r = pert.rvs(trial)
        count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')        

    elif(dist == "Constant"):
        r = most
        # count, bins, ignored = plt.hist(r, 50, density=True, rwidth=0.9, color='teal')

    return {
        "res": r,
        # "val_min": round(min(r), 2),
        # "val_max": round(max(r), 2),
        "mean": round(np.mean(r), 2),
        # "mode": round(mode, 2),
        "p50": round(np.percentile(r, 50), 2),
        "p1": round(np.percentile(r, 1), 2),
        "p10": round(np.percentile(r, 10), 2),
        "p15": round(np.percentile(r, 15), 2),
        "p85": round(np.percentile(r, 85), 2),
        "p90": round(np.percentile(r, 90), 2),
        "p99": round(np.percentile(r, 99), 2),
        # "pdf": [pdf, bins],
        }

def simulation_volumetric() :
    res_area = simulation(float(reservoir_area_e_min.get()), float(reservoir_area_e_most.get()), float(reservoir_area_e_max.get()),
                        clicked_reservoir_area.get())
    res_thick = simulation(float(reservoir_thickness_e_min.get()), float(reservoir_thickness_e_most.get()), float(reservoir_thickness_e_max.get()),
                        clicked_reservoir_thickness.get())
    rock_dens = simulation(float(rock_density_e_min.get()), float(rock_density_e_most.get()), float(rock_density_e_max.get()),
                        clicked_rock_density.get())
    por = simulation(float(porosity_e_min.get()), float(porosity_e_most.get()), float(porosity_e_max.get()),
                        clicked_porosity.get())
    rf = simulation(float(recovery_factor_e_min.get()), float(recovery_factor_e_most.get()), float(recovery_factor_e_max.get()),
                        clicked_recovery_factor.get())
    rock_heat = simulation(float(rock_specific_heat_e_min.get()), float(rock_specific_heat_e_most.get()), float(rock_specific_heat_e_max.get()),
                        clicked_rock_spesific_heat.get())
    res_temp = simulation(float(reservoir_temp_e_min.get()), float(reservoir_temp_e_most.get()), float(reservoir_temp_e_max.get()),
                        clicked_reservoir_temp.get())
    heat_conversion = simulation(float(heat_conversion_e_min.get()), float(heat_conversion_e_most.get()), float(heat_conversion_e_max.get()),
                        clicked_heat_conversion.get())
    plant_life = simulation(float(plant_life_e_min.get()), float(plant_life_e_most.get()), float(plant_life_e_max.get()),
                        clicked_plant_life.get())
    load_factor = simulation(float(load_factor_e_min.get()), float(load_factor_e_most.get()), float(load_factor_e_max.get()),
                        clicked_load_factor.get())
    abandonment_temp = simulation(float(abandonment_temp_e_min.get()), float(abandonment_temp_e_most.get()), float(abandonment_temp_e_max.get()),
                        clicked_abandonment_temp.get())
    ae_initial_rock = (1-por["res"])*(rock_dens["res"]*rock_heat["res"]*res_temp["res"])
    ae_initial_fluid = por["res"]*((783.63*1182.78*1)+(23.71*280.08*0))
    ae_final_rock = (1-por["res"])*(rock_dens["res"]*rock_heat["res"]*abandonment_temp["res"])
    ae_final_fluid = por["res"]*((783.63*1031.97*0.9)+(23.71*232.87*0.1))
    ae_rock = ae_initial_rock-ae_final_rock
    ae_liquid = ae_initial_fluid-ae_final_fluid
    res_areaa = res_area["res"]*(10**6)
    hei = (ae_initial_rock + ae_initial_fluid)*res_thick["res"]*res_areaa
    hef = (ae_final_rock + ae_final_fluid)*res_thick["res"]*res_areaa
    hth = hei-hef
    hde = rf["res"]*hth
    hre = hde/(plant_life["res"]*365*24*3600*1000)
    hel = hde*0.1/(plant_life["res"]*365*24*3600*1000)

    fig, ax1 = plt.subplots()
    count, bins, ignored = ax1.hist(hel, 60, density=False, rwidth=0.9, color='teal')
    data_sorted = np.sort(hel)
    p = 1. * np.arange(len(hel)) / (len(hel) - 1)
    ax1.set_xlabel('Resource Estimation (MWe)', fontsize=15, fontweight='bold')
    ax1.set_ylabel('Frequency', fontsize=15, fontweight='bold')
    ax1.set_title('Volumetric Stored Heat', fontsize=17, fontweight='bold')
    ax1.tick_params(axis='y')

    ax2 = ax1.twinx()
    ax2.plot(data_sorted, p, color='orangered', linewidth=5)
    ax2.tick_params(axis='y')
    ax2.set_ylim(0,1)
    ax2.set_ylabel("Cummulatice Distribution Function", fontsize=15, fontweight='bold')

    canvas = FigureCanvasTkAgg(fig,my_frame1)
    canvas.get_tk_widget().grid(row = 1, rowspan=11, column=11)
    toolbarFrame = tk.Frame(master=my_frame1)
    toolbarFrame.grid(row=12,column=11)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

button_abandonment_temp = tk.Button(my_frame1, text="Calculate Volumetric", command=simulation_volumetric)
button_abandonment_temp.grid(row=12, columnspan=6, sticky='nsew')


root.mainloop()