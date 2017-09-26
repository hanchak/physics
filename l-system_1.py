# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 08:56:15 2017

L-systems physics

@author: HanchaMS
"""

#import numpy as np
#import matplotlib.pyplot as plt
import turtle as tl

# dictionary containing the mapping rules
#rule = {'a':'ab','b':'a'}
#rule = {'X':'F[−X][X]F[−X]+FX', 'F':'FF'}
#rule = {'F':'F-F++F-F'}
rule = {'F':'FF+[+F-F-F]-[-F+F+F]'}
#rule = {'F':'F+F-F-F+F'}

sentence = 'F'
new_sentence = ''

for i in range(3):
    for char in sentence:
        if char == 'F':  #char == 'X' or char == 'F':
            new_sentence += rule[char]
        else:
            new_sentence += char
    
    sentence = new_sentence
    new_sentence = ''
    print(sentence)

# draw the shape using turtle graphics

tl.speed(-1)

l = 15
ang = 20
pos = [(0,-200)]
heading = [90]

# set initial transformation
tl.up()
tl.setpos(pos[0])
tl.setheading(heading[0])
tl.down()


for char in sentence:
    
    if char == 'F':
        
        tl.fd(l)
        
    elif char == '+':
        
        tl.rt(ang)
        
    elif char == '-':
        
        tl.lt(ang)
        
    elif char == '[':
                
        # save transformation state by append position to list
        pos.append(tl.pos())
        heading.append(tl.heading())
        l = l * 2/3

    elif char == ']':
        
        # go to last saved position
        tl.up()
        tl.setpos(pos[-1])
        tl.setheading(heading[-1])
        tl.down()
        l = l*3/2
        
        # delete last saved position
        pos.pop()
        heading.pop()

tl.done()

