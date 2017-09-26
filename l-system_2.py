# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26  2017

L-systems physics... done without turtle graphics

F = draw a line forward
+ = turn right
- = turn left
[ = save position and heading (as a list)
] = restore last saved poisition, pop it from list

@author: Hanchak, Mike
"""
import numpy as np
import matplotlib.pyplot as plt

# dictionary containing the mapping rule
# this is a tree-like L-System (reference: The Coding Train on YouTube)
rule = {'F':'FF+[+F-F-F]-[-F+F+F]'}

# initial condition
sentence = 'F'
new_sentence = ''

# number of recursive iterations
N = 5
#################### generate the sentence of moves #####################
# generate the sentence of moves
for i in range(N):
    for char in sentence:
        if char == 'F':  #char == 'X' or char == 'F':
            new_sentence += rule[char]
        else:
            new_sentence += char
    
    sentence = new_sentence
    new_sentence = ''
    #print(sentence)

#################### pot the resulting moves as lines ####################
######################### using MATPLOTLIB graphics ######################
l = 20  # length of forward movement
ang = 20 * np.pi / 180  # turning angle
pos = [(0,0)]  # initial position (list of tuples)
heading = [90]    # initial heading (list of floats)

# set initial orientation and position of "pen"
x = 0
y = 0
th = 90 * np.pi / 180
scale = 0.66

colors = plt.get_cmap('viridis')

ax = plt.axes()

# now create a list of lists to save all of the x,y points for later plotting
points = [[] for i in range(N+1)]
level = 0  # level for later coloring

def draw(x,y,th,pos):
    '''get the coordinates of the points that make up a single line'''
    oldx = x
    oldy = y
    x = oldx + l*np.cos(th)
    y = oldy + l*np.sin(th)

    return x,y

# go through the sentence on character at a time
for char in sentence:
    
    if char == 'F':
        
        x,y = draw(x,y,th,pos)
        points[level].append([x,y])
        
    elif char == '+':
        
        th -= ang
        
    elif char == '-':
        
        th += ang
        
    elif char == '[':
        level += 1
        # save transformation state by append position to list
        pos.append((x,y))
        heading.append(th)
        l = l * scale
        points[level].append([np.nan,np.nan])  # break in the plotting
        points[level].append([x,y])

    elif char == ']':
        level -= 1
        # go to last saved position
        x,y = pos[-1]
        th = heading[-1]
        l = l / scale
        
        # delete last saved position
        pos.pop()
        heading.pop()
        
        points[level].append([np.nan,np.nan])  # break in the plotting
        points[level].append([x,y])

#### this is the actual plot generation
for i,point in enumerate(points):
    
    pt_array = np.array(point)  # make an array for each level

    plt.plot(pt_array[:,0],pt_array[:,1],linewidth=0.2,color=colors(i / N))


plt.axis('equal')
plt.axis('off')
plt.gcf().set_facecolor('black')
