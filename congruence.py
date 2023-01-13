#!/usr/bin/env python
# coding=utf-8
import numpy as np

def find_solution(a, p):
    a_double = np.double(a)
    p_double = np.double(p)
    for i in range(0, 100000):
        temp_1 = np.double(i) * p_double + np.double(1)
        if (np.mod(temp_1, a_double) == np.double(0)):
            x = temp_1 / a_double
            break
    return x

def legendre(a, p):
    a_double = np.double(a)
    p_double = np.double(p)
    iter = 100000
    if (np.mod(a_double, p_double) == np.double(0)):
        return np.double(0)
    else:
        for i in range(0, iter):
            temp_1 = np.double(i) * p_double + a_double
            x = np.sqrt(temp_1)
            if (np.mod(x, np.double(1.0)) == np.double(0)):
                return np.double(1)
                break
        if (i == iter - 1):
            return np.double(-1)

def compute_left(p = 11):
    temp_2 = np.zeros(p)
    for m in range(0, p):
        temp_1 = np.zeros(p, dtype=np.complex)
        for a in range(0, p):
            frac_up = np.double(m) * np.power(np.double(a), 3, dtype = np.double) + np.power(np.double(a), 2, dtype = np.double)
            frac_down = np.double(p)
            frac = np.double(frac_up / frac_down)
            temp_1[a] = np.exp(np.double(2) * np.pi * 1j * frac)
        temp_sum = temp_1.sum()
        temp_abs = np.absolute(temp_sum)
        temp_2[m] = np.power(temp_abs, 4)
    left = temp_2.sum()
    #print('p = {}, left = {:.4f}'.format(p, left))
    return left

def compute_right(p = 11):
    temp_1 = np.double(2) * np.power(p, 3)
    temp_2 = np.power(p, 2)
    temp_3 = np.double(3)
    temp_4 = np.zeros(p)
    for a in range(1, p):
        frac_up = np.double(a) - np.double(1) + np.double(find_solution(a, p))
        frac_down = np.double(p)
        temp_4[a] = legendre(frac_up, frac_down)
    temp_5 = temp_3 + temp_4[1:].sum()
    temp_6 = temp_5 * temp_2
    right = temp_1 - temp_6
    #print('right = {:.4f}'.format(right))
    return right


p_list = [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103, 107, 127, 131, 139, 151, 167, 179]
for p in p_list:
    left = compute_left(p)
    right = compute_right(p)
    print('p = {}, left = {:.4f}, right = {:.4f}, equal = {}'.format(p, left, right, np.round(left) == np.round(right)))
