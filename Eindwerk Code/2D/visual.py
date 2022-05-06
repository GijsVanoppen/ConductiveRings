import matplotlib.pyplot as plt
from itertools import islice

def nth_index(iterable, value, n):
    matches = (idx for idx, val in enumerate(iterable) if val == value)
    return next(islice(matches, n-1, n), None)

def read_plot_pars(file_name) : #parameters are: box_width, box_length
    f = open(file_name, "r")
    pars = []
    line_counter = 0
    for line in f:
        pars.append(float(line))
        line_counter = line_counter+1
    return pars
    

def draw_circles(file_name):
    f = open(file_name, "r")
    ring_index = 0
    for line in f:
        space_index = nth_index(line, ' ', 1)
        space_index2 = nth_index(line, ' ', 2)
        space_index3 = nth_index(line, ' ', 3)
        circle_x = float(line[0:space_index])
        circle_y = float(line[space_index+1:space_index2])
        r = float(line[space_index2+1:space_index3])
        R = float(line[space_index3+1:-1])
        c = plt.Circle((circle_x, circle_y), radius = r, fill= False)
        plt.text(circle_x,circle_y, str(ring_index)) 
        plt.gca().add_artist(c)
        ring_index = ring_index + 1    
    pars = read_plot_pars("plot_pars.txt")
    plt.axis([0, pars[0], 0, pars[1]])
    plt.gca().set_aspect("equal", adjustable = "box")
    f.close()
    plt.show()

draw_circles("rings_main.txt")