from turtle import color
import matplotlib.pyplot as plt
import math
def distance(p1,p2):
    return (math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2))

def circle_test(circle1, circle2, r):
    x_1 = circle1[0]
    y_1 = circle1[1]
    x_2 = circle2[0]
    y_2 = circle2[1]

    d = distance(circle1, circle2)
    alpha = math.asin((y_2-y_1)/d)

    M_x = x_1 + math.cos(alpha)*0.5*d
    M_y = y_1 + math.sin(alpha)*0.5*d

    junction1 = [M_x - math.sin(alpha)*math.sqrt(r**2-0.25*d**2), M_y + math.cos(alpha)*math.sqrt(r**2-0.25*d**2)]
    junction2 = [M_x + math.sin(alpha)*math.sqrt(r**2-0.25*d**2), M_y - math.cos(alpha)*math.sqrt(r**2-0.25*d**2)]
    #plt.plot(circle1[0], circle1[1], color = 'r')
    #plt.plot(circle2[0], circle2[1], color = 'r')
    plt.scatter([circle1[0], circle2[0]], [circle1[1], circle2[1]])
    plt.scatter(junction1[0], junction1[1])
    plt.scatter(junction2[0], junction2[1])
    plt.show()
c1 = [-0.1, 0.25]
c2 = [1, -0.25]
r = 0.75
circle_test(c1, c2, r)