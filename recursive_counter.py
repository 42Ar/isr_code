#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 19:44:23 2022

@author: frank
"""

def doit(leading, remaining, sol, sols):
    if remaining == 0:
        sols.append(sol)
    else:
        for digit in range(1, min(leading, remaining) + 1):
            doit(leading, remaining - digit, sol + [digit], sols)

def solve_recursive(L):
    sols = []
    for leading in range(1, L):
        doit(leading, L - leading, [leading], sols)
    return sols + [[L]]