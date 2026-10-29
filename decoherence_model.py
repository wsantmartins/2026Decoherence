"""
Shared physics for the central-qubit / Ising-bath model.
Units: J = 1 (h, alpha, T in units of J; t in units of hbar/J), k_B = hbar = 1.
"""
import numpy as np


def Gamma_interacting(t, beta, h, J, alpha, N):
    """Decoherence function, periodic Ising chain (transfer matrix)."""
    C = np.cosh(beta*h - 2j*alpha*t)
    S = np.sqrt(np.exp(-4*beta*J) + np.sinh(beta*h - 2j*alpha*t)**2)
    C0 = np.cosh(beta*h)
    S0 = np.sqrt(np.exp(-4*beta*J) + np.sinh(beta*h)**2)
    return ((C+S)**N + (C-S)**N) / ((C0+S0)**N + (C0-S0)**N)


def Gamma_noninteracting(t, beta, h, alpha, N):
    """J = 0: Gamma = [cosh(beta h - 2i alpha t)/cosh(beta h)]^N."""
    return (np.cosh(beta*h - 2j*alpha*t) / np.cosh(beta*h))**N


def trace_distance(g):
    """BLP trace distance for |+>,|->:  D = |Gamma|."""
    return np.abs(g)


def S2(x, a2=0.5, b2=0.5):
    """Entropy of a dephased qubit with populations a2, b2 and |coherence| factor x."""
    x = np.clip(x, 0, 1)
    lam = 0.5*(1 + np.sqrt((a2-b2)**2 + 4*a2*b2*x**2))
    lam = np.clip(lam, 1e-15, 1-1e-15)
    return -lam*np.log(lam) - (1-lam)*np.log(1-lam)


def I_SF_noninteracting(gamma_abs, f, a2=0.5, b2=0.5):
    """J = 0 closed form: I(S:F) = S2(|Gamma|) - S2(|Gamma|^(1-f))."""
    return S2(gamma_abs, a2, b2) - S2(gamma_abs**(1-f), a2, b2)


def boundary_gammas(t, beta, h, J, alpha, L):
    """
    Boundary-conditioned decoherence functions of an open segment of L spins
    whose outer neighbours are fixed to s, s' = +-1:
        Gamma^{ss'}(t) = W_{ss'}(t) / W_{ss'}(0),   W = B (V B)^L,
    with B the bond matrix and V = diag(e^{x}, e^{-x}), x = beta h - 2i alpha t.
    Returns (Gamma matrix, W(0) matrix).
    """
    B = np.array([[np.exp(beta*J), np.exp(-beta*J)],
                  [np.exp(-beta*J), np.exp(beta*J)]], dtype=complex)

    def W(x):
        V = np.diag([np.exp(x), np.exp(-x)])
        M = B.copy()
        for _ in range(L):
            M = M @ V @ B
        return M

    W0 = W(beta*h + 0j)
    return W(beta*h - 2j*alpha*t) / W0, W0


def I_SF_interacting(t, beta, h, J, alpha, N, k, a2=0.5, b2=0.5):
    """
    Exact mutual information between the qubit and a contiguous fragment of k spins
    (periodic chain of N spins):
        I(S:F) = S2(|Gamma|) - sum_{ss'} P(s,s') S2(|Gamma^{ss'}_{rest}|),
    where s, s' are the two end spins of the fragment.
    """
    G, W0 = boundary_gammas(t, beta, h, J, alpha, N - k)
    B = np.array([[np.exp(beta*J), np.exp(-beta*J)],
                  [np.exp(-beta*J), np.exp(beta*J)]])
    V0 = np.diag([np.exp(beta*h), np.exp(-beta*h)])
    A = V0.copy()                      # fragment weight with end spins (s', s)
    for _ in range(k - 1):
        A = A @ B @ V0
    P = np.real(A * W0.T)
    P /= P.sum()
    I = S2(abs(Gamma_interacting(t, beta, h, J, alpha, N)), a2, b2)
    for sp in range(2):
        for s in range(2):
            I -= P[sp, s] * S2(abs(G[s, sp]), a2, b2)
    return I


def correlation_length(beta, J=1.0):
    """Zero-field Ising correlation length (lattice units)."""
    return 1.0 / np.log(1.0/np.tanh(beta*J))
