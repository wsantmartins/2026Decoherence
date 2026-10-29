"""Figure 5 (J = 1, h = 0): (a) normalized PIP at the first critical time t1 = pi/(4 alpha N);
(b) single-spin normalized redundancy vs xi/N for several N."""
import numpy as np
import matplotlib.pyplot as plt
import style
plt.rcParams.update({"font.size": 7})
from decoherence_model import Gamma_interacting, S2, I_SF_interacting, correlation_length

alpha, J, h = 0.1, 1.0, 0.0

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(3.4, 1.8))

N = 10
t1 = np.pi / (4*alpha*N)
betas  = [0.5, 1.0, 1.5, 2.0, 4.0]
colors = ['black', 'purple', 'blue', 'orange', 'red']
for beta, c in zip(betas, colors):
    SS = S2(np.abs(Gamma_interacting(t1, beta, h, J, alpha, N)))
    R = [0.0] + [I_SF_interacting(t1, beta, h, J, alpha, N, k) / SS for k in range(1, N+1)]
    ax1.plot(range(N+1), R, 'o-', color=c, lw=1.2, ms=2.2, label=rf'$\beta={beta:g}$')
ax1.plot([0, N], [0, 1], 'k--', lw=0.8)
ax1.set_xlim(0, N); ax1.set_ylim(0, 1.03); ax1.set_xticks(range(0, N+1, 2))
ax1.set_xlabel(r'$fN$'); ax1.set_ylabel(r'$\mathcal{I}(\mathcal{S}:\mathcal{F})/S(\rho_S)$')
ax1.legend(fontsize=5.5, frameon=False, loc='lower right')

beta_grid = np.linspace(0.3, 3.5, 120)
for N, c, ls in zip([6, 10, 20, 40], ['black', 'blue', 'orange', 'red'], ['-', '--', '-.', ':']):
    t1 = np.pi / (4*alpha*N)
    R1 = [I_SF_interacting(t1, b, h, J, alpha, N, 1) / S2(np.abs(Gamma_interacting(t1, b, h, J, alpha, N)))
          for b in beta_grid]
    ax2.semilogx(correlation_length(beta_grid) / N, R1, color=c, ls=ls, lw=1.1, label=rf'$N={N}$')
ax2.set_xlim(1e-2, 1e2); ax2.set_ylim(0, 1.03)
ax2.set_xlabel(r'$\xi/N$'); ax2.set_ylabel(r'$\mathcal{I}(\mathcal{S}:\mathcal{F}_1)/S(\rho_S)$')
ax2.legend(fontsize=5.5, frameon=False, loc='upper left')
ax1.text(0.04, 0.9, '(a)', transform=ax1.transAxes); ax2.text(0.88, 0.06, '(b)', transform=ax2.transAxes)
fig.tight_layout(pad=0.3, w_pad=0.8)
fig.savefig('figure5.pdf')
