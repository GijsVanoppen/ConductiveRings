import matplotlib.pyplot as plt
from itertools import islice
import numpy as np

def save(filepath, fig=None):
    if not fig:
        fig = plt.gcf()
    plt.subplots_adjust(0,0,1,1,0,0)
    fig.savefig(filepath,  bbox_inches='tight')


def nth_index(iterable, value, n):
    matches = (idx for idx, val in enumerate(iterable) if val == value)
    return next(islice(matches, n-1, n), None)

def read_plot_pars() : 
    """ 
    Parameters are (in this order):
    Minimum box width,
    Maximum box width,
    Box width iterations,
    Minimum box length,
    Maximum box length,
    Box length iterations,
    Voltage difference
    """
    def extract_parameter_from_string(string):
        #returns the part of the string after the ':' sign
        parameter = ""
        start_index = string.find(':')
        for i in range(start_index+1, len(string)-1):
            parameter += string[i]
        return parameter

    f = open("input.txt", "r")
    pars = []

    line_counter = 0
    for line in f:
        if ((line_counter > 0) and (line_counter < 8)): 
            pars.append(extract_parameter_from_string(line)) 
        line_counter += 1



    return pars
    
def draw_rings(file_name, show_ring_index):
    f = open(file_name, "r")
    for line in f:
        #first, get positions of spaces in the line so we can extract the parameters
        space_index = nth_index(line, ' ', 1)
        space_index2 = nth_index(line, ' ', 2)
        space_index3 = nth_index(line, ' ', 3)
        space_index4 = nth_index(line, ' ', 4)
        
        #extract parameters 
        ring_x = float(line[space_index+1:space_index2])
        ring_y = float(line[space_index2+1:space_index3])
        r = float(line[space_index3+1:space_index4])
        R = float(line[space_index4+1:-1])
        
        #insert ring_index if requested
        if show_ring_index:
            ring_index = int(line[0:space_index])
            plt.text(ring_x,ring_y, str(ring_index)) 

        #draw ring
        c = plt.Circle((ring_x, ring_y), radius = r, fill= False)
        plt.gca().add_artist(c)
    f.close()
    
    #plotting
    pars = read_plot_pars()
    plt.axis([0, float(pars[1]), 0, float(pars[4])])
    plt.gca().set_aspect("equal", adjustable = "box")
    plt.xlabel("x (a.u.)")
    plt.ylabel("y (a.u.)")
    save("images\\rings.png")


def plotV_of_x(junctions_file_name, results_file_name):
    junctions = []
    results = []
    #read x-values of junctions
    f = open(junctions_file_name)
    for line in f:
        space_index = nth_index(line, ' ', 1)
        junctions.append(float(line[0:space_index]))
    f.close()

    #read voltages (results)
    f = open(results_file_name)
    for line in f:
        results.append(float(line))
    f.close()

    #linear fit
    coef = np.polyfit(junctions, results,1)
    poly1d_fn = np.poly1d(coef)

    #plotting and text box
    fig, ax = plt.subplots()
    ax.plot(junctions, results, 'r+', junctions, poly1d_fn(junctions), 'k')
    textstr = "Linear fit (V=mx+b):\n" + "m = " + str(round(coef[0], 10)) + "\nb = " + str(round(coef[1], 11))
    props = dict(boxstyle='round', facecolor='white', alpha=1)
    ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize = 14, verticalalignment = 'top', bbox = props)
    plt.xlabel("x (a.u.)")
    plt.ylabel("Voltage (a.u.)")
    pars = read_plot_pars()
    ax.set_ylim([0, float(pars[6])])
    save("images\V_of_x.png")

def plotV_of_xy(junctions_file_name, results_file_name):
    junctions_x = []
    junctions_y = []
    results = []
    #read x- and y-values of junctions
    f = open(junctions_file_name)
    for line in f:
        space_index = nth_index(line, ' ', 1)
        junctions_x.append(float(line[0:space_index]))
        junctions_y.append(float(line[space_index+1:-1]))
    f.close()

    #read voltages (results)
    f = open(results_file_name)
    for line in f:
        results.append(float(line))
    f.close()
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(junctions_x, junctions_y, results, c = 'r', marker = '+')
    ax.set_xlabel('x (a.u.)')
    ax.set_ylabel('y (a.u.)')
    ax.set_zlabel('Voltage (a.u.)')
    save("images/V_of_xy.png")

def plotR_of_box_dimensions(resistances_file_name, width_or_length):
    
    #read resistances
    resistances = []
    f = open(resistances_file_name)
    for line in f:
        resistances.append(float(line))
    f.close()


    if width_or_length == "width":
        #get width values
        width_list = []
        pars = read_plot_pars()
        box_width_min = float(pars[0])
        box_width_max = float(pars[1])
        box_width_iterations = int(pars[2])
        
        #put width values in list to plot them
        box_width = box_width_min
        if (box_width_iterations!=1):    
            for _ in range(box_width_iterations):
                width_list.append(box_width)
                box_width += (box_width_max-box_width_min)/(box_width_iterations-1)
        else:
            width_list.append(box_width)
        #plotting
        plt.plot(width_list, resistances, 'r+')

    elif width_or_length == "length":
        #get length values
        length_list = []
        pars = read_plot_pars()
        box_length_min = float(pars[3])
        box_length_max = float(pars[4])
        box_length_iterations = int(pars[5])
        
        #put length values in list to plot them
        box_length = box_length_min
        for _ in range(box_length_iterations):
            length_list.append(box_length)
            box_length += (box_length_max-box_length_min)/box_length_iterations

        #plotting
        plt.plot(length_list, resistances, 'r+')

    else:
        print("Failed to plot R. Second parameter should be \"width\" or \"length\"")
    
    plt.xlabel("Box " + width_or_length + " (a.u.)")
    plt.ylabel("Resistance (a.u.)")
    save("images/R_of_dims.png")
    






draw_rings("rings_all.txt", False)
# plotV_of_x("junctions.txt", "results.txt")
# plotV_of_xy("junctions.txt", "results.txt")
# plotR_of_box_dimensions("resistances.txt", "width")