"""Figure 4 (J = 0): (a) D(t) = |Gamma(t)|; (b) normalized PIP I(S:F)/S(rho_S) at t = tau/2."""
import numpy as np
import matplotlib.pyplot as plt
import style
plt.rcParams.update({"font.size": 7})
from decoherence_model import Gamma_noninteracting, S2, I_SF_noninteracting

N, alpha, h = 10, 0.1, 1.0
tau = np.pi / (2*alpha)
t = np.linspace(0, tau, 1500)
betas  = [0.1, 0.5, 1.0, 2.0, 4.0]
colors = ['black', 'purple', 'blue', 'orange', 'red']
k = np.arange(0, N+1)              # fragment sizes fN

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(3.4, 1.8))
for beta, c in zip(betas, colors):
    ax1.plot(t, np.abs(Gamma_noninteracting(t, beta, h, alpha, N)), color=c, lw=1.0, label=rf'$\beta={beta:g}$')
    g = np.abs(Gamma_noninteracting(tau/2, beta, h, alpha, N))
    R = I_SF_noninteracting(g, k/N) / S2(g)
    ax2.plot(k, R, 'o-', color=c, lw=1.2, ms=2.2)
ax2.plot([0, N], [0, 1], 'k--', lw=0.8)          # proportional share R = f
ax1.set_xlim(0, tau); ax1.set_ylim(0, 1.03)
ax1.set_xlabel(r'$t$'); ax1.set_ylabel(r'$D(\rho_{+}^{S},\rho_{-}^{S})$')
ax2.set_xlim(0, N); ax2.set_ylim(0, 1.03); ax2.set_xticks(range(0, N+1, 2))
ax2.set_xlabel(r'$fN$'); ax2.set_ylabel(r'$\mathcal{I}(\mathcal{S}:\mathcal{F})/S(\rho_S)$')
ax1.legend(fontsize=5.5, frameon=False, loc='center', bbox_to_anchor=(0.5, 0.42))
ax1.text(0.9, 0.06, '(a)', transform=ax1.transAxes); ax2.text(0.04, 0.9, '(b)', transform=ax2.transAxes)
fig.tight_layout(pad=0.3, w_pad=0.8)
fig.savefig('figure4.pdf')
