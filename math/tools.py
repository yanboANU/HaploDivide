import os
import sys

def get_factorial(n):
    if n <= 1:
        return 1
    else:
        return n * get_factorial(n-1)


def get_combination_number(n,m):
    first = get_factorial(n)
    second = get_factorial(m)
    third = get_factorial(n-m)
    return first/(second*third)
