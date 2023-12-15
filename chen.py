import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
x = []
d = 8
p = []
p = np.eye(8)
saiz = 9
for i in range(d):
    x.append(np.linspace(0, 2 * np.pi, saiz))
v = []
"""
for x1 in x[0]:
    for x2 in x[1]:
        for x3 in x[2]:
            for x4 in x[3]:
                for x5 in x[4]:
                    for x6 in x[5]:
                        for x7 in x[6]:
                            for x8 in x[7]:
                                x1
"""
for i in range(saiz):
    z = 0
    for j in range(len(x)):
        z += x[j][i] * p[j][j]
    print
    v.append(z)


def combinations_of_size(input_set, k):
    # Convert the set to a list to make it compatible with combinations
    elements = list(input_set)
    
    # Generate combinations of size k
    result = list(combinations(elements, k))
    
    return result
my_set = {1, 2, 3, 4,5,6,7,8}
k_value = 2
result_combinations = combinations_of_size(my_set, k_value)
print(result_combinations)
#for (a,b) in result_combinations:

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[0], x[1], v)

ax.set_xlabel('X1 Axis')
ax.set_ylabel('X2 Axis')
ax.set_zlabel('V Axis')
plt.show()