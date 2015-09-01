import matplotlib.pyplot as plt
import numpy
from numpy import *
from enum import Enum

class ItemType(Enum):
    TIAuto = 0
    TILinear = 1
    TIConstant = 2
    TILogistic = 5
    TIExponential = 6


def get_segment_in_arr(x, arr):
    if x < arr[0]:
        return -1;
    i = 0
    for a in arr:
        if i != 0 and a > x:
            return i - 1
        i += 1
    return i - 1

# tests

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
####v_in = [1, 21, 30, 34, 35, 41]
##v_in = [100, 200, 300, 400, 500, 600]
##v_out = [1, 2, 1, 1.5, 0.5, 1]
##v_def = [1.1, 2, 1, 0.97, 1, 1]
##v_type = [ItemType.TIIncrease, ItemType.TILogistic, ItemType.TILogistic, ItemType.TIDecrease, ItemType.TILinear, ItemType.TIConstant]

v_in = [0, 1e-4, 2e-4, 4e-4]
v_out = [1e-6, 2e-6, 4e-6, 2e-6]
v_def = [1.0, 0.1, 1.0, 1.0]
v_type = [ItemType.TIExponential, ItemType.TIExponential, ItemType.TIExponential, ItemType.TILinear]

##v_in = [0, 1e-4]
##v_out = [2e-6, 1e-6]
##v_def = [1.0, 1.0]
##v_type = [ItemType.TIExponential, ItemType.TIExponential]

def calc_point(x):
    i = get_segment_in_arr(x, v_in)
    if i < 0:
        return -1.0
    last_point = i + 1 >= len(v_in)
    next_auto = False if last_point else v_type[i + 1] == ItemType.TIAuto
    no_next_point = last_point or next_auto
    
    #print(i)
    type = v_type[i]
    if no_next_point:
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
    elif type == ItemType.TIExponential:
        k = v_def[i]
        x0 = v_in[i]
        y0 = v_out[i]
        x1 = v_in[i + 1]
        y1 = v_out[i + 1]
        k0 = 1.
        b = 0.
        par_sign = 1 if y1 > y0 else -1

        dx = x1 - x0;
        dy = y1 - y0;

        def_dx = 1.0
        def_dy = 1.0
##        print("par_sign = ", par_sign)
        x_targ = par_sign * (-0.8) * def_dx
        y_targ = def_dy / 100.0

        k0 = y_targ ** (par_sign / x_targ)
        k *= k0
        #print("k = ", k)

        coef_x = abs(dx / def_dx)
        coef_y = abs(dy / def_dy)
        print("coef_x = ", coef_x)
        print("coef_y = ", coef_y)

        # compute x_shift
        pow0 = k ** (par_sign * (-def_dx))
        d_pow = par_sign * (1. - pow0)
        print("d_pow = ", d_pow)
        print("def_dy = ", def_dy)
        print("def_dy / d_pow = ", def_dy / d_pow)
        x_shift = -par_sign * math.log(def_dy / d_pow) / math.log(k) + def_dx
        print("k = ", k)
        print("x_shift = ", x_shift)

        def f(x):
            x -= x0
            x /= coef_x
            ret = k ** (par_sign * (x - x_shift))
            return ret * coef_y + b
##        print("y0 / coef_y = ", y0 / coef_y)
##        print("f(x0) = ", f(x0))
##        print("f(x0) / coef_y = ", f(x0) / coef_y)
        b = y0 - f(x0)
##        print("b = ", b)
        return f(x)
    elif type == ItemType.TILogistic:
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
#x = v_in[0] - indent
x = v_in[0]
delta_graph = (v_in[1] - v_in[0]) / 10
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
