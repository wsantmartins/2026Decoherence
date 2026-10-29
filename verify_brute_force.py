"""
Checks the closed-form results of decoherence_model.py against brute-force
construction and diagonalization of the qubit-fragment state rho_SF.

  - I_SF_interacting  (Eq. 33, any J, contiguous fragment)
  - I_SF_noninteracting (Eq. 34, J = 0)
  - the bound I(S:F)/S(rho_S) <= f for J = 0 (Eq. 35)

The bath is diagonal in the configuration basis, so rho_SF is built exactly by
summing over all 2^N bath configurations. Run: python verify_brute_force.py
"""
import itertools
import numpy as np
from decoherence_model import (I_SF_interacting, I_SF_noninteracting,
                               Gamma_noninteracting, S2)


def vn_entropy(rho):
    w = np.linalg.eigvalsh(rho)
    w = w[w > 1e-14]
    return -np.sum(w*np.log(w))


def I_SF_brute(t, beta, h, J, alpha, N, k, a2=0.5):
    """Mutual information between the qubit and spins 0..k-1 of a periodic chain."""
    confs = np.array(list(itertools.product([1, -1], repeat=N)))
    E = -J*np.sum(confs*np.roll(confs, -1, axis=1), axis=1) - h*confs.sum(axis=1)
    p = np.exp(-beta*(E - E.min()))
    p /= p.sum()
    phase = np.exp(-2j*alpha*t*confs.sum(axis=1))

    # index of the fragment configuration (spins 0..k-1)
    idx = ((1 - confs[:, :k])//2) @ (2**np.arange(k)[::-1])
    D = 2**k
    pF = np.zeros(D)
    cF = np.zeros(D, complex)
    np.add.at(pF, idx, p)
    np.add.at(cF, idx, p*phase)

    b2 = 1 - a2
    ab = np.sqrt(a2*b2)
    rhoSF = np.zeros((2*D, 2*D), complex)
    rhoSF[:D, :D] = np.diag(a2*pF)
    rhoSF[D:, D:] = np.diag(b2*pF)
    rhoSF[:D, D:] = np.diag(ab*cF)
    rhoSF[D:, :D] = np.diag(ab*np.conj(cF))
    Gamma = np.sum(p*phase)
    rhoS = np.array([[a2, ab*Gamma], [ab*np.conj(Gamma), b2]])
    return vn_entropy(rhoS) + vn_entropy(np.diag(pF)) - vn_entropy(rhoSF)


if __name__ == "__main__":
    N, alpha = 8, 0.1
    worst = 0.0
    for beta, h, J in [(1, 1, 1), (0.5, 0.3, 1), (2, 0, 1), (1, 1, 0), (4, 0, 1)]:
        for t in [2.0, 7.85, 11.0]:
            for k in range(1, N):
                for a2 in [0.5, 0.2]:
                    exact = I_SF_interacting(t, beta, h, J, alpha, N, k, a2, 1 - a2)
                    brute = I_SF_brute(t, beta, h, J, alpha, N, k, a2)
                    worst = max(worst, abs(exact - brute))
                    if J == 0:
                        g = abs(Gamma_noninteracting(t, beta, h, alpha, N))
                        closed = I_SF_noninteracting(g, k/N, a2, 1 - a2)
                        worst = max(worst, abs(closed - brute))
    print(f"Eqs. (33)-(34) vs brute force, N={N}: max |difference| = {worst:.1e}")

    g = np.linspace(1e-4, 1 - 1e-7, 20001)
    ok = all(np.all(1 - S2(g**(1 - f))/S2(g) <= f + 1e-9) for f in np.linspace(0.01, 0.99, 99))
    print(f"Bound I/S(rho_S) <= f for J = 0 holds on the grid: {ok}")
