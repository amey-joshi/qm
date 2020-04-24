import math
import numpy as np

# Values of constants in CGS units.
e = 4.770e-10 # in e.s.u
a_0 = 0.5282e-8 # in cm

def get_alphas(Z):
    a1 = 0
    x1 = 5*Z - 3/2 
    x2 = math.sqrt(Z**2 - 3*Z + 45/8)

    a2 = (x1 - x2)/6
    a3 = (x1 + x2)/6

    return np.array([a1, a2, a3])

print(f'alphas for Z = 2: {get_alphas(2)}')
print(f'alphas for Z = 3: {get_alphas(3)}')
print(f'alphas for Z = 4: {get_alphas(4)}')

def f(Z):
    alphas = get_alphas(Z)

    x1 = 2 * np.multiply(np.power(alphas, 2), np.power(Z - alphas, 2))
    x2 = 9/64 * np.power(alphas, 2)
    x3 = np.multiply(np.power(alphas, 2), Z - alphas)

    return x1 - x2 - x3

K = e**4/a_0**2

Delta_sq_He = K * f(2)
Delta_sq_Li = K * f(3)
Delta_sq_Be = K * f(4)

print(f'Delta_sq_He = {Delta_sq_He}')
print(f'Delta_sq_Li = {Delta_sq_Li}')
print(f'Delta_sq_Be = {Delta_sq_Be}')

Delta_He = math.sqrt(Delta_sq_He[1])
Delta_Li = math.sqrt(Delta_sq_Li[1])
Delta_Be = math.sqrt(Delta_sq_Be[1])

# Ground state of Hydrogen atom.
W_H = e**2/(2 * a_0)

def E(Z):
    alpha = get_alphas(Z)[1]
    x1 = -2 * alpha**2
    x2 = 4 * alpha * (alpha - Z)
    x3 = 5 * alpha/4

    return (x1 + x2 + x3) * W_H

E_He = E(2)
E_Li = E(3)
E_Be = E(4)

# Convert ergs to eV.
to_eV = 0.6285E12
print(f'Energy estimate for He is {E_He * to_eV} eV.')
print(f'Energy estimate for Li is {E_Li * to_eV} eV.')
print(f'Energy estimate for Be is {E_Be * to_eV} eV.')

bounds_He = [E_He - Delta_He, E_He + Delta_He]
bounds_Li = [E_Li - Delta_Li, E_Li + Delta_Li]
bounds_Be = [E_Be - Delta_Be, E_Be + Delta_Be]

print(f'Bounds for He are [{E_He - Delta_He}, {E_He + Delta_He}].')
print(f'Bounds for Li are [{E_Li - Delta_Li}, {E_Li + Delta_Li}].')
print(f'Bounds for Be are [{E_Be - Delta_Be}, {E_Be + Delta_Be}].')

bounds_He = np.array([E_He - Delta_He, E_He + Delta_He]) * to_eV
bounds_Li = np.array([E_Li - Delta_Li, E_Li + Delta_Li]) * to_eV
bounds_Be = np.array([E_Be - Delta_Be, E_Be + Delta_Be]) * to_eV

print(f'Bounds for He are {bounds_He}.')
print(f'Bounds for Li are {bounds_Li}.')
print(f'Bounds for Be are {bounds_Be}.')
