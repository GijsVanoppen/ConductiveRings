from turtle import distance
import cairo
import math


WIDTH, HEIGHT = 1024, 1024

surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
ctx = cairo.Context(surface)
ctx.scale(WIDTH, HEIGHT)


circles = []

def get_circles():
    f = open("circles.txt", "r")
    for line in f:                          #lines are built like this:  x, y, r
        space_index = line.index(" ")       #index of the 'space' character
        space_index2 = line.rindex(" ")
        circle_x = float(line[0:space_index])
        circle_y = float(line[space_index+1:space_index2])
        r = float(line[space_index2+1:-1])
        circles.append([circle_x, circle_y, r])
    f.close()
    return circles

circles = get_circles()

def get_V():
    V = []
    f = open("V.txt")
    for line in f:
        V.append(float(line))
    return V

V = get_V()

def dist(c_1, c_2):
    return math.sqrt((c_1[0]-c_2[0])**2 + (c_1[1]-c_2[1])**2)


def calc_junctions(circles, i):
    r = circles[i][-1]
    x_i = circles[i][0]    
    y_i = circles[i][1]        
    x_prev = circles[i-1][0]
    y_prev = circles[i-1][1]
    x_next = circles[i+1][0]
    y_next = circles[i+1][1]

    #between first 2 circles
    d = dist(circles[i-1], circles[i])
    alpha = math.asin((y_i - y_prev)/d)
    M_x = x_prev + math.cos(alpha)*d*0.5 
    M_y = y_prev + math.sin(alpha)*d*0.5

    junction_1_x = M_x - math.sin(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_1_y = M_y + math.cos(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_2_x = M_x + math.sin(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_2_y = M_y - math.cos(alpha)*math.sqrt(r**2 - 0.25*d**2)
    
    #between last 2 circles
    d = dist(circles[i], circles[i+1])
    alpha = math.asin((y_next - y_i)/d)
    M_x = x_i + math.cos(alpha)*d*0.5
    M_y = y_i + math.sin(alpha)*d*0.5

    junction_3_x = M_x - math.sin(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_3_y = M_y + math.cos(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_4_x = M_x + math.sin(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_4_y = M_y - math.cos(alpha)*math.sqrt(r**2 - 0.25*d**2)

    #put the coörds in arrays
    junction_1 = [junction_1_x, junction_1_y]
    junction_2 = [junction_2_x, junction_2_y]
    junction_3 = [junction_3_x, junction_3_y]
    junction_4 = [junction_4_x, junction_4_y]

    #put the arrays in an array, to be returned
    junctions = [junction_1, junction_2, junction_3, junction_4]
    return junctions



def voltage_to_rgb(voltage):
    r = 1 - voltage
    b = voltage
    return [r, 0, b]


"""
def draw_arc():
    ctx.arc(0.5, 0.5, 0.05, 3*math.pi/2, 0)
    ctx.set_source_rgba(0, 0, 0, 0)

    ctx.set_line_width(0.005)
    pattern = cairo.LinearGradient(0.5, 0.45, 0.5, 0.55)
    pattern.add_color_stop_rgb(1, 1, 0.5, 0.5)

    pattern.add_color_stop_rgb(0.1, 0.5, 0.5, 1)
    ctx.set_source(pattern)
    ctx.stroke()

draw_arc()
"""


def draw_circle(circles, i, V):
    #x_1 = circles[i][0]
    #x_2 = circles[i][1]
    #y_1 = circles[i][0]
    #y_2 = circles[i][1]
    
    c_x = circles[i][0]
    c_y = circles[i][1]
    max_x = circles[-1][0]
    r = circles[0][2]

    junctions = calc_junctions(circles, i)
    
    #UPPER HORIZONTAL WIRE
    x_1 = junctions[0][0]       #upper left node (1)
    y_1 = junctions[0][1]       #upper left node (1)
    x_2 = junctions[2][0]       #upper right node (3)
    y_2 = junctions[2][1]       #upper right node (3)
    beta = math.asin(y_2-y_1)/(x_2- x_1)          #angle between horizontal and the line connection the two points
    #Apply rotation
    ctx.save()
    ctx.translate(c_x/max_x, c_y/max_x)
    ctx.rotate(math.pi/2 - beta)
    ctx.translate(-c_x/max_x, -c_y/max_x)
    # Define gradient
    rgb_1 = voltage_to_rgb(V[i*4-2])        #rgb of upper left node
    rgb_2 = voltage_to_rgb(V[i*4])          #rgb of upper right node 
    pattern = cairo.LinearGradient(x_1, y_1, x_2, y_2)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Restore user space
    ctx.restore()
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x, math.pi+math.atan((y_1-c_y)/(c_x-x_1)), 2*math.pi-math.atan((y_2-c_y)/(x_2-c_x)))
    #ctx.arc(0.5, 0.5, 0.05, 3*math.pi/2, 0)
    ctx.set_line_width(0.005)
    ctx.set_source(pattern)
    ctx.stroke()


    #BOTTOM HORIZONTAL WIRE
    x_1 = junctions[1][0]       #upper left node (2)
    y_1 = junctions[1][1]       #upper left node (2)
    x_2 = junctions[3][0]       #upper right node (4)
    y_2 = junctions[3][1]       #upper right node (4)
    beta = math.asin(y_2-y_1)/(x_2- x_1)          #angle between horizontal and the line connection the two points
    #Apply rotation
    ctx.save()
    ctx.translate(c_x/max_x, c_y/max_x)
    ctx.rotate(math.pi/2 - beta)
    ctx.translate(-c_x/max_x, -c_y/max_x)
    # Define gradient
    rgb_1 = voltage_to_rgb(V[i*4-1])        #rgb of bottom left node
    rgb_2 = voltage_to_rgb(V[i*4+1])          #rgb of bottom right node 
    pattern = cairo.LinearGradient(x_1, y_1, x_2, y_2)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Restore user space
    ctx.restore()
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x, math.atan((c_y-y_2)/(x_2-c_x)), 0.5*math.pi+math.atan((c_x - x_1)/(c_y-y_1)))
    #ctx.arc(0.5, 0.5, 0.05, 3*math.pi/2, 0)
    ctx.set_line_width(0.005)
    ctx.set_source(pattern)
    ctx.stroke()






#draw_arc(circles, 1, [0, 0, 0.9, 0, 0.8, 0])

for i in range(1,9):
    draw_circle(circles, i, V)









surface.write_to_png("1D.png")