import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import pandas as pd

# Paramètres physiques
gamma = 1.4
R     = 287.0

# Conditions initiales gauche
rho_L = 1.0
T_L   = 1000.0
u_L   = 0.0

# Conditions initiales droite
rho_R = 0.125
T_R   = 300.0
u_R   = 0.0

# Calcul de PL et PR

PL = rho_L*R*T_L
PR = rho_R*R*T_R

# Calcul de CL et CR

CL = np.sqrt(gamma*PL/rho_L)
CR = np.sqrt(gamma*PR/rho_R)

# definition d'une focntion f = differences de U* = 0

def f (P_star):
    u_star_gauche = u_L + (2*CL)/(gamma -1)*(1 - (P_star/PL)**((gamma - 1)/(2*gamma)))
    u_star_droite = u_R + (P_star - PR) * np.sqrt((2/(gamma+1))/(rho_R*(P_star + (gamma - 1)/(gamma +1)*PR)))
    return u_star_gauche - u_star_droite

#résolution

P_star = brentq(f, PR, PL)
u_star = u_L + (2*CL)/(gamma-1) * (1 - (P_star/PL)**((gamma-1)/(2*gamma)))

print("P_star =", P_star)
print("u_star =", u_star)

print("f(P_star) =", f(P_star))
print("f(PR) =", f(PR))
print("f(PL) =", f(PL))

#calcul des etats des zones

#zone 2 apres détente, avant contact
rho_2 = rho_L* (P_star/PL)**((1/gamma))

#zone 3 après contact, avant choc
rho_3 = rho_R * (P_star/PR) * (((gamma-1)/(gamma+1)) + P_star/PR) / (1 + ((gamma-1)/(gamma+1)) * P_star/PR)

print("rho_2 =", rho_2)
print("rho_3 =", rho_3)

S_head = u_L - CL
c_2 = np.sqrt(gamma*P_star/rho_2)
S_tail = u_star - c_2
S_contact = u_star
S_shock = u_R + CR *np.sqrt((gamma+1)/(2*gamma) * P_star/PR + (gamma-1)/(2*gamma))

print("S_head =", S_head)
print("S_contact =", S_contact)
print("S_shock =", S_shock)
print("S_tail =", S_tail)


def solution_analytique(x, t):
    xi = (x - 1.0) / t  # vitesse de similarité

    if xi < S_head:
        # Zone 1 — conditions initiales gauche
        return rho_L, u_L, PL

    elif xi < S_tail:
        # Zone détente — profil continu
        U =  2 / (gamma +1) * (xi + CL)
        C = CL - (gamma - 1 )/2 * U
        p1 = PL * (C/CL)**((2*gamma)/(gamma - 1 ))
        rho1 = rho_L * (p1/PL)**(1/gamma)
        return rho1, U, p1

    elif xi < S_contact:
        # Zone 2
        return rho_2, u_star, P_star

    elif xi < S_shock:
        # Zone 3
        return rho_3, u_star, P_star

    else:
        # Zone 4 — conditions initiales droite
        return rho_R, u_R, PR

t = 0.002  # ton t_max

x_test = [0.5, 0.8, 1.0, 1.2, 1.5, 1.8]
for x in x_test:
    rho, u, p = solution_analytique(x, t)
    print(f"x={x:.1f} | rho={rho:.3f} | u={u:.1f} | p={p:.1f}")


print("Position S_head   :", 1.0 + S_head * t)
print("Position S_tail   :", 1.0 + S_tail * t)
print("Position S_contact:", 1.0 + S_contact * t)
print("Position S_shock  :", 1.0 + S_shock * t)


# Charge le dernier CSV de ta simulation
df = pd.read_csv("../cmake-build-debug/resultat_220.csv")
# Construit la solution analytique sur tout le domaine
x_ana = np.linspace(0, 2, 1000)
rho_ana, u_ana, p_ana = [], [], []

for x in x_ana:
    r, u, p = solution_analytique(x, t)
    rho_ana.append(r)
    u_ana.append(u)
    p_ana.append(p)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

ax1.plot(x_ana, rho_ana, label='analytique')
ax1.plot(df['x'], df['rho'], label='numérique')

ax2.plot(x_ana, u_ana, label='analytique')
ax2.plot(df['x'], df['vitesse'], label='numérique')

ax3.plot(x_ana, p_ana, label='analytique')
ax3.plot(df['x'], df['pression'], label='numérique')

ax1.set_ylabel("Densité (kg/m³)")
ax2.set_ylabel("Vitesse (m/s)")
ax3.set_ylabel("Pression (Pa)")
ax3.set_xlabel("Position x (m)")
ax1.legend()
ax1.set_title(f"Validation Sod — t = {t} s")

plt.show()

rho_ana_num = []
for x in df['x'].values:
    r, u, p = solution_analytique(x, t)
    rho_ana_num.append(r)

rho_ana_num = np.array(rho_ana_num)
erreur_rho = np.sqrt(np.mean((df['rho'].values - rho_ana_num)**2))
print(f"Erreur L2 densité : {erreur_rho:.4f}")
