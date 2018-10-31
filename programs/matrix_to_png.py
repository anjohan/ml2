import numpy as np
import matplotlib.pyplot as plt
import sys

method = sys.argv[1]

f = open("data/J_" + method + ".dat", "r")
lines = f.readlines()
L = int(lines[0])

n = 0

for i in range(2, len(lines), L + 1):
    n += 1
    x = np.zeros((L, L))
    for j in range(L):
        x[j, :] = list(map(float, lines[i + j].split()))
    print(x)
    plt.imsave(
        "data/J_" + method + str(n) + ".png",
        x,
        cmap=plt.get_cmap("Greys"),
        vmin=0,
        vmax=1)
