import matplotlib.pyplot as plt
def draw_circles(file_name, r):
    f = open(file_name, "r")
    for line in f:      #lines are built like this:  circle_x circle_y r
        space_index = line.index(" ")
        space_index2 = line.rindex(" ")
        circle_x = float(line[0:space_index])
        circle_y = float(line[space_index+1:space_index2])
        r = float(line[space_index2+1:-1])
        c = plt.Circle((circle_x, circle_y), radius = r, fill= False)
        print(circle_x, circle_y, r)
        plt.gca().add_artist(c)
    plt.axis([0, circle_x, -3*r, 3*r])
    plt.gca().set_aspect("equal", adjustable = "box")
    f.close()
    plt.show()

draw_circles("circles.txt", 0.75)