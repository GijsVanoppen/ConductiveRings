import matplotlib.pyplot as plt
from itertools import islice

def nth_index(iterable, value, n):
    matches = (idx for idx, val in enumerate(iterable) if val == value)
    return next(islice(matches, n-1, n), None)


def draw_circles(file_name, r):
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
        #print(circle_x, circle_y, r)
        plt.gca().add_artist(c)
        ring_index = ring_index + 1      
    plt.axis([0, 10, 0, 10])
    plt.gca().set_aspect("equal", adjustable = "box")
    f.close()
    plt.show()

draw_circles("rings.txt", 0.75)