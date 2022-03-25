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


def calc_junctions_first(circles):
    r = circles[0][-1]
    x_i = circles[0][0]    
    y_i = circles[0][1]        
    x_next = circles[1][0]
    y_next = circles[1][1]

    d = dist(circles[0], circles[1])
    alpha = math.asin((y_next - y_i)/d)
    M_x = x_i + math.cos(alpha)*d*0.5
    M_y = y_i + math.sin(alpha)*d*0.5

    junction_3_x = M_x - math.sin(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_3_y = M_y + math.cos(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_4_x = M_x + math.sin(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_4_y = M_y - math.cos(alpha)*math.sqrt(r**2 - 0.25*d**2)

    #put the coörds in arrays
    junction_1 = [0, r]
    junction_2 = [0, -r]
    junction_3 = [junction_3_x, junction_3_y]
    junction_4 = [junction_4_x, junction_4_y]

    #put the arrays in an array, to be returned
    junctions = [junction_1, junction_2, junction_3, junction_4]
    return junctions

def calc_junctions_last(circles):
    N = len(circles)

    r = circles[N-1][-1]
    x_i = circles[N-1][0]    
    y_i = circles[N-1][1]        
    x_prev = circles[N-2][0]
    y_prev = circles[N-2][1]


    #between first 2 circles
    d = dist(circles[N-2], circles[N-1])
    alpha = math.asin((y_i - y_prev)/d)
    M_x = x_prev + math.cos(alpha)*d*0.5 
    M_y = y_prev + math.sin(alpha)*d*0.5

    junction_1_x = M_x - math.sin(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_1_y = M_y + math.cos(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_2_x = M_x + math.sin(alpha)*math.sqrt(r**2 - 0.25*d**2)
    junction_2_y = M_y - math.cos(alpha)*math.sqrt(r**2 - 0.25*d**2)

    #put the coörds in arrays
    junction_1 = [junction_1_x, junction_1_y]
    junction_2 = [junction_2_x, junction_2_y]
    junction_3 = [0, r]
    junction_4 = [0, -r]

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

def draw_upper_arc(circles, i, V):
    c_x = circles[i][0]
    c_y = circles[i][1]
    max_x = circles[-1][0]
    r = circles[0][2]

    junctions = calc_junctions(circles, i)
    
    x_1 = junctions[0][0]       #upper left node (1)
    y_1 = junctions[0][1]       #upper left node (1)
    x_2 = junctions[2][0]       #upper right node (3)
    y_2 = junctions[2][1]       #upper right node (3)
    # Define gradient
    rgb_1 = voltage_to_rgb(V[i*4-2])        #rgb of upper left node
    rgb_2 = voltage_to_rgb(V[i*4])          #rgb of upper right node 
    pattern = cairo.LinearGradient(x_1/max_x, y_1/max_x, x_2/max_x, y_2/max_x)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x, math.pi+math.atan((y_1-c_y)/(c_x-x_1)), 2*math.pi-math.atan((y_2-c_y)/(x_2-c_x)))

    ctx.set_line_width(0.05/len(circles))
    ctx.set_source(pattern)
    ctx.stroke()


def draw_bottom_arc(circles, i, V):    
    c_x = circles[i][0]
    c_y = circles[i][1]
    max_x = circles[-1][0]
    r = circles[0][2]

    junctions = calc_junctions(circles, i)
    
    x_1 = junctions[1][0]       #bottom left node (2)
    y_1 = junctions[1][1]       #bottom left node (2)
    x_2 = junctions[3][0]       #bottom right node (4)
    y_2 = junctions[3][1]       #bottom right node (4)
    # Define gradient
    rgb_1 = voltage_to_rgb(V[i*4-1])        #rgb of bottom left node
    rgb_2 = voltage_to_rgb(V[i*4+1])          #rgb of bottom right node 
    pattern = cairo.LinearGradient(x_1/max_x, y_1/max_x, x_2/max_x, y_2/max_x)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x, math.atan((c_y-y_2)/(x_2-c_x)), 0.5*math.pi+math.atan((c_x - x_1)/(c_y-y_1)))
    ctx.set_line_width(0.05/len(circles))
    ctx.set_source(pattern)
    ctx.stroke()    
    


def draw_left_arc(circles, i, V):
    c_x = circles[i][0]
    c_y = circles[i][1]
    max_x = circles[-1][0]
    r = circles[0][2]

    junctions = calc_junctions(circles, i)
        
    x_1 = junctions[0][0]       #upper left node (1)
    y_1 = junctions[0][1]       #upper left node (1)
    x_2 = junctions[1][0]       #bottom left node (2)
    y_2 = junctions[1][1]       #bottom left node (2)
    # Define gradient
    rgb_1 = voltage_to_rgb(V[i*4-2])        #rgb of upper left node
    rgb_2 = voltage_to_rgb(V[i*4-1])          #rgb of bottom left node 
    pattern = cairo.LinearGradient(x_1/max_x, y_1/max_x, x_2/max_x, y_2/max_x)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x,0.5*math.pi+math.atan((c_x-x_2)/(c_y-y_2)), math.pi+math.atan((y_1-c_y)/(c_x-x_1)))
    ctx.set_line_width(0.05/len(circles))
    ctx.set_source(pattern)
    ctx.stroke()    
    

def draw_right_arc(circles, i, V):
    c_x = circles[i][0]
    c_y = circles[i][1]
    max_x = circles[-1][0]
    r = circles[0][2]

    junctions = calc_junctions(circles, i)
        
    x_1 = junctions[2][0]       #upper right node (3)
    y_1 = junctions[2][1]       #upper right node (3)
    x_2 = junctions[3][0]       #bottom right node (4)
    y_2 = junctions[3][1]       #bottom right node (4)
    # Define gradient
    rgb_1 = voltage_to_rgb(V[i*4])            #rgb of upper right node
    rgb_2 = voltage_to_rgb(V[i*4+1])          #rgb of bottom right node 
    pattern = cairo.LinearGradient(x_1/max_x, y_1/max_x, x_2/max_x, y_2/max_x)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x, 2*math.pi-math.atan((y_1-c_y)/(x_1-c_x)), math.atan((c_y-y_2)/(x_2- c_x)))
    ctx.set_line_width(0.05/len(circles))
    ctx.set_source(pattern)
    ctx.stroke()   


def draw_first_circle(circles):
    c_x = circles[0][0]
    c_y = circles[0][1]
    max_x = circles[-1][0]
    r = circles[0][2]

    junctions = calc_junctions_first(circles)
    # UPPER HORIZONTAL WIRE
    x_1 = junctions[0][0]       #upper left node (1)
    y_1 = junctions[0][1]       #upper left node (1)
    x_2 = junctions[2][0]       #upper right node (3)
    y_2 = junctions[2][1]       #upper right node (3)
    # Define gradient
    rgb_1 = voltage_to_rgb(1)        #rgb of upper left node
    rgb_2 = voltage_to_rgb(V[0])          #rgb of upper right node 
    pattern = cairo.LinearGradient(x_1/max_x, y_1/max_x, x_2/max_x, y_2/max_x)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x, 1.5*math.pi, 2*math.pi-math.atan((y_2-c_y)/(x_2-c_x)))

    ctx.set_line_width(0.05/len(circles))
    ctx.set_source(pattern)
    ctx.stroke()


    # HORIZONTAL RIGHT WIRE
    x_1 = junctions[2][0]       #upper right node (3)
    y_1 = junctions[2][1]       #upper right node (3)
    x_2 = junctions[3][0]       #bottom right node (4)
    y_2 = junctions[3][1]       #bottom right node (4)
    # Define gradient
    rgb_1 = voltage_to_rgb(V[0])            #rgb of upper right node
    rgb_2 = voltage_to_rgb(V[1])          #rgb of bottom right node 
    pattern = cairo.LinearGradient(x_1/max_x, y_1/max_x, x_2/max_x, y_2/max_x)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x, 2*math.pi-math.atan((y_1-c_y)/(x_1-c_x)), math.atan((c_y-y_2)/(x_2- c_x)))
    ctx.set_line_width(0.05/len(circles))
    ctx.set_source(pattern)
    ctx.stroke()   

    # DRAW BOTTOM HORIZONTAL LINE
    x_1 = junctions[1][0]       #bottom left node (2)
    y_1 = junctions[1][1]       #bottom left node (2)
    x_2 = junctions[3][0]       #bottom right node (4)
    y_2 = junctions[3][1]       #bottom right node (4)
    # Define gradient
    rgb_1 = voltage_to_rgb(1)        #rgb of bottom left node
    rgb_2 = voltage_to_rgb(V[1])          #rgb of bottom right node 
    pattern = cairo.LinearGradient(x_1/max_x, y_1/max_x, x_2/max_x, y_2/max_x)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x, math.atan((c_y-y_2)/(x_2-c_x)), -0.5*math.pi)
    ctx.set_line_width(0.05/len(circles))
    ctx.set_source(pattern)
    ctx.stroke()    

    

def draw_last_circle(circles):
    N = len(circles)
    c_x = circles[N-1][0]
    c_y = circles[N-1][1]
    max_x = circles[-1][0]
    r = circles[0][2]

    junctions = calc_junctions_last(circles)
    
    # UPPER WIRE
    x_1 = junctions[0][0]       #upper left node (1)
    y_1 = junctions[0][1]       #upper left node (1)
    x_2 = junctions[2][0]       #upper right node (3)
    y_2 = junctions[2][1]       #upper right node (3)
    # Define gradient
    rgb_1 = voltage_to_rgb(V[-2])        #rgb of upper left node
    rgb_2 = voltage_to_rgb(0)          #rgb of upper right node 
    pattern = cairo.LinearGradient(x_1/max_x, y_1/max_x, x_2/max_x, y_2/max_x)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x, math.pi+math.atan((y_1-c_y)/(c_x-x_1)), 1.5*math.pi)
    ctx.set_line_width(0.05/len(circles))
    ctx.set_source(pattern)
    ctx.stroke()

    # LEFT WIRE
    x_1 = junctions[0][0]       #upper left node (1)
    y_1 = junctions[0][1]       #upper left node (1)
    x_2 = junctions[1][0]       #bottom left node (2)
    y_2 = junctions[1][1]       #bottom left node (2)
    # Define gradient
    rgb_1 = voltage_to_rgb(V[-2])        #rgb of upper left node
    rgb_2 = voltage_to_rgb(V[-1])          #rgb of bottom left node 
    pattern = cairo.LinearGradient(x_1/max_x, y_1/max_x, x_2/max_x, y_2/max_x)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x,0.5*math.pi+math.atan((c_x-x_2)/(c_y-y_2)), math.pi+math.atan((y_1-c_y)/(c_x-x_1)))
    ctx.set_line_width(0.05/len(circles))
    ctx.set_source(pattern)
    ctx.stroke()   

    # BOTTOM WIRE
    x_1 = junctions[1][0]       #bottom left node (2)
    y_1 = junctions[1][1]       #bottom left node (2)
    x_2 = junctions[3][0]       #bottom right node (4)
    y_2 = junctions[3][1]       #bottom right node (4)
    # Define gradient
    rgb_1 = voltage_to_rgb(V[-1])        #rgb of bottom left node
    rgb_2 = voltage_to_rgb(0)          #rgb of bottom right node 
    pattern = cairo.LinearGradient(x_1/max_x, y_1/max_x, x_2/max_x, y_2/max_x)
    pattern.add_color_stop_rgb(0, rgb_1[0], 0, rgb_1[2])
    pattern.add_color_stop_rgb(1, rgb_2[0], 0, rgb_2[2])
    # Draw the circle
    ctx.arc(c_x/max_x, -c_y/max_x+0.5, r/max_x, 0.5*math.pi, 0.5*math.pi+math.atan((c_x - x_1)/(c_y-y_1)))
    ctx.set_line_width(0.05/len(circles))
    ctx.set_source(pattern)
    ctx.stroke()    
     


for i in range(1,len(circles)-1):
    draw_upper_arc(circles, i, V)
    draw_bottom_arc(circles, i, V)
    draw_left_arc(circles, i, V)
    draw_right_arc(circles, i, V)
draw_first_circle(circles)
draw_last_circle(circles)







surface.write_to_png("1D.png")