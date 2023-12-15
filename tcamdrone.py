import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations

def V(x,P,G):
    sum = 0
    for i in G:
        sum = sum + x.T @ P[:][:][i] @ x
    return sum

def h(psi):
    return [[(np.sin(2*psi)/4 - 1/2)*((538*np.cos(psi)**2)/747 + (1285*np.sin(psi)**2)/747 - 1285/747)],
            [-(np.sin(2*psi)/4 + 1/2)*((538*np.cos(psi)**2)/747 + (1285*np.sin(psi)**2)/747 - 1285/747)],
            [-(np.sin(2*psi)/4 - 1/2)*((538*np.cos(psi)**2)/747 + (1285*np.sin(psi)**2)/747 - 538/747)],
            [(np.sin(2*psi)/4 + 1/2)*((538*np.cos(psi)**2)/747 + (1285*np.sin(psi)**2)/747 - 538/747)]]

def dh(psi):
    return [[(np.cos(2*psi)*((538*np.cos(psi)**2)/747 + (1285*np.sin(psi)**2)/747 - 1285/747))/2 + 2*np.cos(psi)*np.sin(psi)*(np.sin(2*psi)/4 - 1/2)],
            [ - (np.cos(2*psi)*((538*np.cos(psi)**2)/747 + (1285*np.sin(psi)**2)/747 - 1285/747))/2 - 2*np.cos(psi)*np.sin(psi)*(np.sin(2*psi)/4 + 1/2)],
            [- (np.cos(2*psi)*((538*np.cos(psi)**2)/747 + (1285*np.sin(psi)**2)/747 - 538/747))/2 - 2*np.cos(psi)*np.sin(psi)*(np.sin(2*psi)/4 - 1/2)],
            [(np.cos(2*psi)*((538*np.cos(psi)**2)/747 + (1285*np.sin(psi)**2)/747 - 538/747))/2 + 2*np.cos(psi)*np.sin(psi)*(np.sin(2*psi)/4 + 1/2)]]

def combinations_of_size(input_set, k):
    # Convert the set to a list to make it compatible with combinations
    elements = list(input_set)
    
    # Generate combinations of size k
    result = list(combinations(elements, k))
    
    return result

# Example usage:
my_set = {0,1, 2, 3, 4,5,6,7}
k_value = 2
result_combinations = combinations_of_size(my_set, k_value)
print(result_combinations)

x = np.array([3, 0, 0, 0, 0, 0, 0, 0])
P = [[np.identity(8)],[np.identity(8)],[np.identity(8)]]
G =  [0] # G offset by 1 because python starts at 0

saiz = 6
V_graph = np.zeros([saiz**2, saiz**2])
x_graph = np.zeros(saiz**2)
y_graph = np.zeros(saiz**2)
fig = plt.figure()
for i,j in result_combinations:
    x=np.zeros(8)
    vix=0
    viy=0
    for xi in np.linspace(0, 2 * np.pi, saiz):
        for xj in np.linspace(0, 2 * np.pi, saiz):
            x[i] = xi
            x[j] = xj
            V_graph[vix][viy]= V(x,P,G)
            x_graph[vix]=xi
            y_graph[viy]=xj
            viy = viy + 1
        viy=0
        vix = vix + 1
            
    ax = fig.add_subplot(111, projection='3d')
    print(y_graph)
    print(x_graph)
    ax.plot(x_graph,y_graph,V_graph)
    ax.set_xlabel('X1 Axis')
    ax.set_ylabel('X2 Axis')
    ax.set_zlabel('V Axis')
    plt.show()
    
    
    



for x1 in np.linspace(0, 2 * np.pi, saiz):
    for x2 in np.linspace(0, 2 * np.pi, saiz):
        for x3 in np.linspace(0, 2 * np.pi, saiz):
            for x4 in np.linspace(0, 2 * np.pi, saiz):
                for x5 in np.linspace(0, 2 * np.pi, saiz):
                    for x6 in np.linspace(0, 2 * np.pi, saiz):
                        for x7 in np.linspace(0, 2 * np.pi, saiz):
                            for x8 in np.linspace(0, 2 * np.pi, saiz):
                                x = np.array([x1,x2,x3,x4,x5,x6,x7,x8])
                                V_graph.append(V(x,P,G))







fig, ax = plt.subplots()
ax.plot(V_graph)
plt.show()


print(V(x,P,G))
print(h(2*np.pi))
