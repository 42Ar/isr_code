#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 15:26:04 2022

@author: frank
"""

import recursive_counter

L = 10
sols = recursive_counter.solve_recursive(L)
sols_check = []
for sol in sols:
    res = ""
    for i, c in enumerate(sol):
        if i%2 == 0:
            res += "1"*c
        else:
            res += "0"*c
    if res[-1] == "0":
        sols_check.append(res)
assert len(sols_check) == len(set(sols_check))

#%%

sols = []
L = 10
alternating = sum(1 << i for i in range(0, 64, 2))
for leading in list(range(2, L//2))[::-1]:
    L_algo = L - leading # the remaining leading bits are brute-forced
    sol = [None]*(L_algo - 2)
    c = L_algo - leading
    b = (1 << (leading - 1)) - 1
    for i in range(L_algo - leading):
        sol[i] = 1
    runs = 0
    b = b << (c + 1)
    c += 1
    while True:
        b = ((b >> c) << c) | (alternating >> (63 - c + (runs&1)))
        runs += c - 2
        c = sol[runs]
        if (runs&1) == 0:
            start = 1 << (c - 1)
            if start == 1:
                start += 1
            end = (1 << leading) - 2
        else:
            start = 0
            end = (1 << leading) - 1 - (1 << (c - 1))
        print(f"{b:b} {start:b} {end:b}")
        start |= (b << leading)
        end |= (b << leading)
        for lastbits in range(start, end + 1, 2):
            print(f"{lastbits:02b} {[leading] + sol[:runs+1]}")
            sols.append(tuple([leading] + sol[:runs+1] + [f"{lastbits:b}"]))
        if runs == 0:
            break
        sol[runs] = 1
        prev = sol[runs - 1]
        need_break = False
        while prev == leading:
            runs -= 1
            if runs == 0:
                need_break = True
                break
            c += prev
            sol[runs] = 1
            prev = sol[runs - 1]
        if need_break:
            break
        sol[runs - 1] = prev + 1
start = ((1 << L//2) - 1) << L//2
end = (1 << L) - 2
for brute_bits in range(start, end + 1, 2):
    sols.append(("brute", f"{brute_bits:02b}"))
sols.append(("last", f"{(alternating >> (63 - L)):b}"))
assert set(sols_check) == set(t[-1] for t in sols)

#%%

# leading | code         stack
#       5 | 5            0
#       5 | 4 1          0
#       5 | 3 2          0, 1
#       5 | 3 1 1        0
#       5 | 2 3          0, 1
#       5 | 2 2 1        0, 2
#       5 | 2 1 2        0, 2
#       5 | 2 1 1 1      0
#       5 | 1 4          1
#       5 | 1 3 1        1
#       5 | 1 2 2        1, 2
#       5 | 1 2 1 1      1
#       5 | 1 1 3        2
#       5 | 1 1 2 1      2
#       5 | 1 1 1 2      3
#       5 | 1 1 1 1 1

# 6 | 3, RUNS((1, 1, 1)), 3 | 4, RUNS(3, (1, 1, 1)), 2 | 5, RUNS(:5), 1


L = 10
sols = []
for leading in range(2, L):
    stack = [0]*20
    sol = [0]*20
    sums = [0]*20
    remaining = L - leading
    i = 0
    while remaining > leading:
        sol[i] = leading
        remaining -= leading
        sums[i] = leading*i + leading
        stack[i] = i
        i += 1
    sums[i] = leading*i + remaining
    sol[i] = remaining
    print(leading, sums, sol, stack)
    stacksize = i
    runi = i
    while True:
        print(leading, runi, sol[:runi + 1], stack[:stacksize])
        sols.append(tuple([leading] + sol[:runi + 1]))
        if sol[runi] == 1:
            # can no longer decrease => get next runi from stack
            if stacksize == 0:
                break
            stacksize -= 1
            runi = stack[stacksize]
        sol[runi] -= 1
        sums[runi] -= 1
        remaining = L - leading - sums[runi]
        runi += 1
        sol[runi] = remaining
        sums[runi] = sums[runi - 1] + remaining
        if sol[runi - 1] > 1:
            # only save current index to stack if it has remaining
            stack[stacksize] = runi - 1
            stacksize += 1
    break