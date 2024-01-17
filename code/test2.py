import matplotlib as plt
import matplotlib.cm as pltcm
import numpy as np


z = np.array([(1,2), (3,4), (5,6)])

a = np.tile(z, (1,3))
b = np.reshape(a, (np.shape(z)[0] * 3, np.shape(z)[1]))

print(b)

"""
triple = 3
    a = np.shape(node_plt)
    repeat = np.tile(single, (1, triple))
    repeat = np.reshape(repeat, (np.shape(single)[0] * 3, np.shape(single)[1]))
    tripled = np.array(list(map(tuple, repeat)))
    #repeat = np.repeat(single, triple)
    #tripled = repeat.reshape((np.shape(node_plt)))
    """