import math
import numpy as np

# The two roots of the quadratic in $\alpha$ are
roots = np.roots([13, 98, 21])

# Calculate normalization constants. First the factors in them with
# \alpha^2 a in it.
factors = np.power(roots, 2) + 6*roots + 21
# The two normalization constants are:
N = np.sqrt(315/(16 * factors))

# Calculate E(\alpha).
E = 0.75 * np.divide(11*np.power(roots, 2) + 14*roots + 35, np.power(roots, 2) + 6 * roots + 21)

# Calculate D(\alpha).
D = np.multiply(np.power(N, 2), 42*np.power(roots, 2)/5 + 4*roots + 2)

# Calculate \Delta(\alpha).
Delta = np.sqrt(D - np.power(E, 2))

# Ground state Delta
Delta_g = np.min(Delta)

# Ground state energy
E_g = np.min(E)

# Energy bounds
bounds = [E_g - Delta_g, E_g + Delta_g]

# Width as percentage
E_0 =  math.pi**2/8 # True energe in units of \hslash^2/(ma^2).
pct = Delta_g/E_0

