#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 15:53:03 2022

@author: frank
"""

import numpy as np
import matplotlib.pyplot as plt

def check_code(c, k, h, Delta):
    L = len(c)
    cc = c*3
    r = sum((cc[L + j + Delta]*cc[L + j]*
             cc[L + j + Delta + h - k]*cc[L + j + h - k])
            for j in range(L))
    if r != 0:
        print("code error:", k, h, Delta, r)


def full_check(c, N, M):
    for k in range(N):
        for h in range(N):
            if k != h:
                for Delta in range(1, M+1):
                    check_code(c, k, h, Delta)

sols = []
f = open("brute_code/solution_summary")
L_vals = set()
for line in f:
    L = int(line.split(",")[1])
    N = int(line.split(",")[2])
    M = int(line.split(",")[3])
    c = line.split("[")[1].split("]")[0].split(",")
    code = [int(cc) for cc in c]
    full_check(code, N, M-1)
    assert len(code) == L
    if N==8 and M==6:
        print(code)
    print(N, M, sum(code))
    sols.append((N, M, code))
max_N = max(t[0] for t in sols)
max_M = max(t[1] for t in sols)
res = np.zeros((max_M, max_N-1), dtype=np.int32)
for N in range(2, max_N + 1):
    for M in range(1, max_M + 1):
        try:
            L = len(next(t for t in sols if t[0] == N and t[1] == M)[2])
        except StopIteration:
            L = -1
        res[M-1, N-2] = L
        L_vals.add(L)
res_c = np.zeros((max_M, max_N-1), dtype=np.int32)
c = list(L_vals)
c.sort()
for N in range(2, max_N + 1):
    for M in range(1, max_M + 1):
        res_c[M-1, N-2] = c.index(res[M-1, N-2])
cmap = plt.cm.get_cmap('tab10', len(L_vals))
cmap.set_under((0.9, 0.9, 0.9))
plt.imshow(res_c, cmap=cmap, origin="lower", interpolation="none")
plt.xticks(range(max_N-1), range(2, max_N + 1))
plt.yticks(range(max_M), range(1, max_M+1))
for N in range(2, max_N+1):
    for M in range(1, max_M+1):
        if res[M-1, N-2] != -1:
            plt.text(N-2, M-1, str(res[M-1, N-2]), va='center', ha='center')
        elif N >= 8:
            plt.text(N-2, M-1, ">48", va='center', ha='center')
plt.clim(1)
plt.xlabel("N")
plt.ylabel("M")
#plt.savefig("code_length_overview.pdf")
plt.show()


#%%


def doit(leading, remaining, sol, sols):
    if remaining == 0:
        #if len(sol) % 2 == 1:
        sols.append(sol)
    else:
        for digit in range(1, min(leading, remaining) + 1):
            doit(leading, remaining - digit, sol + [digit], sols)

x = []
y = []
for L in range(2, 20):
    sols = []
    for leading in range(1, L):
        doit(leading, L - leading, [leading], sols)
    x.append(L)
    if L == 10:
        f = sols
        print(sols)
    print(L, len(sols))
    y.append(len(sols)/L**2)
plt.plot(x, y)

#%%

sols = []
L = 10
alternating = sum(1 << i for i in range(0, 64, 2))
for leading in range(2, L):
    c = L - leading
    sol = [1]*(L - leading)
    b = (1 << (leading - 1)) - 1
    runs = 0
    i = 0
    while True:
        if b % 2 == 0:
            b = (b << (c + 1)) | (alternating >> (63 - c))
        else:
            b = (b << (c + 1)) | (alternating >> (62 - c))
        runs += c
        c = 0
        print(f"{leading} {b:b} {[leading] + sol}")
        sols.append([leading] + sol[:runs])
        while runs > 0 and (sol[runs - 1] == leading or c == 0):
            c += sol[runs - 1]
            b = b >> sol[runs - 1]
            sol[runs - 1] = 1
            runs -= 1
        if runs == 0:
            break
        sol[runs - 1] += 1
        c -= 1
        i += 1



