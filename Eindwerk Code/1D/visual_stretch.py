import matplotlib.pyplot as plt
def draw_graph(file_name, r):   
    #draw the R_tot(a)-graph
    a_lst = []
    R_tot_lst = []
    f = open(file_name, "r")
    for line in f:                          #lines are built like this:  a R_tot
        space_index = line.index(" ")       #index of the 'space' character
        a = float(line[0:space_index])
        R_tot = float(line[space_index+1:-1])
        a_lst.append(a)
        R_tot_lst.append(R_tot)
        print(a, " ", R_tot)
    plt.xlabel("Lattice constant a (a.u.)")
    plt.ylabel("Total resistance R_tot (a.u.)")
    plt.scatter(a_lst, R_tot_lst)
    f.close()
    plt.show()

draw_graph("results_stretch.txt", 0.75)