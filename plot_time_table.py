import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy
from numpy import *
import os
import shutil
from enum import Enum

class ItemType(Enum):
    TIAuto = 0
    TILinear = 1
    TIConstant = 2
    TIIncrease = 3
    TIDecrease = 4
    TILogistic = 5


def get_segment_in_arr(x, arr):
    if x < arr[0]:
        return -1;
    i = 0
    for a in arr:
        if i != 0 and a > x:
            return i - 1
        i += 1
    return i - 1

##v_in = [0, 1, 2, 3]
##v_out = [1, 1, 2, 2.5]
##v_def = [1, 3.1, 1, 1]
##v_type = [1, 4, 1, 1]

##v_in = [10, 100, 200, 300, 400]
###v_in = [1, 2, 3, 4, 5]
##v_out = [1, 0, 1, 3, 5]
##v_def = [1, 1, 8, 2.1, 1]
###v_type = [ItemType.TIConstant, ItemType.TIAuto, ItemType.TILogistic, ItemType.TIIncrease, ItemType.TILinear]
###v_type = [ItemType.TIConstant, ItemType.TIAuto, ItemType.TILogistic, ItemType.TILinear, ItemType.TILinear]

##v_in = [1, 2]
###v_in = [100, 200]
##v_out = [1, 2]
##v_def = [1, 20]
##v_type = [ItemType.TILogistic, ItemType.TILogistic]

#v_in = [1, 2, 3, 4, 5, 6]
v_in = [1, 21, 30, 34, 35, 41]
#v_in = [100, 200, 300, 400, 500, 600]
v_out = [1, 2, 1, 1.5, 0.5, 1, 0.2]
v_def = [1, 1, 1, 1, 1, 1]
v_type = [ItemType.TIIncrease, ItemType.TILogistic, ItemType.TILogistic, ItemType.TIDecrease, ItemType.TIIncrease, ItemType.TIConstant]

def calc_point(x):
    i = get_segment_in_arr(x, v_in)
    if i < 0:
        return -1.0
    last_point = i + 1 >= len(v_in)
    #print(i)
    type = v_type[i]
    if last_point:
        type = ItemType.TIConstant
    if type == ItemType.TIAuto:
        return -1.0
    elif type == ItemType.TILinear:
        #print("i = ", i)
        x_local = (x - v_in[i]) / (v_in[i + 1] - v_in[i])
        #print("x_local = ", x_local)
        return x_local * (v_out[i + 1] - v_out[i]) + v_out[i]
    elif type == ItemType.TIConstant:
        return v_out[i]
    elif type == ItemType.TIIncrease or type == ItemType.TIDecrease:
        par_sign = 1 if type == ItemType.TIIncrease else -1
        k = v_def[i]
        x0 = v_in[i]
        y0 = v_out[i]
        x1 = v_in[i + 1]
        y1 = v_out[i + 1]
        k0 = 1.
        b = 0.

        # compute k0
        #x_targ = (v_in[i + 1] - v_in[i]) / 2.0
        x_targ = par_sign * (-0.8) * (v_in[i + 1] - v_in[i])
        y_targ = math.fabs(v_out[i + 1] - v_out[i]) / 100.0
##        if y_targ <= 0:
##            print("can't use " + ("increase" if type == ItemType.TIIncrease else "decrease"))
##            return -1.
##        print("y_targ = ", y_targ)
##        print("x_targ = ", x_targ)
        k0 = y_targ ** (par_sign / x_targ)

        # set k
        k *= k0
##        print ("k = ", k)
##        print ("x1 = ", x1)

        # compute x_shift
        dy = y1 - y0
        dx = math.pow(k, par_sign * x1) - math.pow(k, par_sign * x0)
        dx = k ** (par_sign * x1) - k ** (par_sign * x0)
        x_shift = 0.
        if dx != 0:
            x_shift = -par_sign * math.log(dy / dx, k)

        def f(x):
            return k ** (par_sign * (x - x_shift)) + b
        b = y0 - f(x0)
        return f(x)
    elif ItemType.TILogistic:
        # for logistic we are computing k0, x_shift, L, b
        # for last point: L = 1, k0 = 1
        k = v_def[i]
        L = v_out[i + 1] - v_out[i]
        
        # compute k0
        x_targ = v_in[i + 1] - v_in[i]
        y_coef = 0.99999   # determines function form
        y_targ = y_coef * L
        k0 = -math.log(L / y_targ - 1) / x_targ

        # set k
        k *= k0
        
        x_shift = v_in[i]
        # compute x_shift
        x_shift = (v_in[i] + v_in[i + 1]) / 2.0
        
        x0 = x_shift # checkpoint to set b
        y0 = (v_out[i] + v_out[i + 1]) / 2.0      
        b = 0.        
        def f(x):
            return L / (1 + math.exp(-k * (x - x_shift))) + b
        b = y0 - f(x0)
        return f(x)
    return -1.0


v_x = []
v_y = []
indent = v_in[1] - v_in[0]
#indent = 0
x = v_in[0] - indent
delta_graph = (v_in[1] - v_in[0]) / 1000
while x <= v_in[-1] + indent:
    v_x.append(x)
    v_y.append(calc_point(x))
    x += delta_graph

ax = plt.gca()
#ax.set_xticks(numpy.arange(min(time), max(time), find_rounded_delta(time)))
#ax.set_yticks(numpy.arange(min(area.points), max(area.points), find_rounded_delta(area.points)))
plt.grid()

im = plt.plot(v_x, v_y, marker='.', color='b')
plt.show()
#plt.savefig('time_table.png', dpi=100)
plt.close()
