# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import recursive_counter

x, y, y_comb = [], [], []
for L in range(3, 25):
    sols = recursive_counter.solve_recursive(L)
    if L < 10:
        print(sols)
    x.append(L)
    y.append(2**(L - 2)/(len(sols)/2 - 1))
    y_comb.append(len(sols)/2 - 1)
    print(L, len(sols))
x.append(34)
y.append(2**(34 - 2)/358728827)
x.append(36)
y.append(2**(36 - 2)/1356530952)

#%%

# runs with N=6, M=8
plt.scatter([32, 34, 36], [6.5/0.612, 20.45/1.91, 105/9], label="single core, exhaustive")
plt.scatter([32, 34, 36], [1.194/0.134, 3.417/0.368, 18.34/1.668], label="multi core, exhaustive")
plt.plot(x, y, label="theoretical maximum", c="black")
plt.grid()
plt.legend()
plt.title("speedup on server (enum vs. loops impl), 128bit")
plt.grid()
plt.savefig("speedup.pdf")
plt.show()

#%%

# factor three penalty for using 128bit vector instead of 64bit
plt.scatter([32, 34, 36], [1.194, 3.417, 18.34], label="exhaustive, 128bit, enum", marker="^")
plt.scatter([32, 34, 36], [0.134, 0.368, 1.668], label="exhaustive, 128bit, loops")
plt.scatter([48], [8*60 + 55], label="N=8, M=6, 128bit, enum", marker="^")
plt.scatter([48], [5*60 + 8.4], label="N=8, M=6, 128bit, loops")
plt.scatter([48], [1*60 + 44], label="N=8, M=6, 64bit, loops")
plt.scatter([40], [0.189], label="N=7, M=5, 128bit, enum", marker="^")
plt.scatter([40], [0.514], label="N=6, M=8, 128bit, loops")
plt.yscale("log")
plt.axhline(3600*24, label="1 day")
plt.grid()
plt.xlim(31, 60)
plt.ylim(0.1, 10e5)
plt.title("runtime, multi-core")
plt.legend(loc="upper left")
plt.savefig("perf_plot.pdf")
plt.show()