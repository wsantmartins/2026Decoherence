"""Figure 2: trace distance D(t) = |Gamma(t)| for four (J, h) settings and four temperatures."""
import numpy as np
import matplotlib.pyplot as plt
import style
plt.rcParams.update({"font.size": 7})
from decoherence_model import Gamma_interacting, Gamma_noninteracting

N, alpha = 10, 0.1
tau = np.pi / (2*alpha)
t = np.linspace(0, tau, 3000)
betas  = [0.0, 0.75, 1.75, 4.0]
colors = ['black', 'blue', 'purple', 'red']
panels = [('(a)', 1.0, 1.0), ('(b)', 1.0, 0.1), ('(c)', 1.0, 0.0), ('(d)', 0.0, 1.0)]  # (label, J, h)

fig, axes = plt.subplots(1, 4, figsize=(7.0, 1.9), sharey=True)
for ax, (lab, J, h) in zip(axes, panels):
    for beta, c in zip(betas, colors):
        G = Gamma_interacting(t, beta, h, J, alpha, N) if J > 0 else Gamma_noninteracting(t, beta, h, alpha, N)
        ax.plot(t, np.abs(G), color=c, lw=1.0, label=rf'$\beta={beta:g}$')
    ax.set_xlim(0, tau); ax.set_ylim(0, 1.03); ax.set_xlabel(r'$t$')
    ax.text(0.03, 0.06, lab, transform=ax.transAxes, bbox=dict(fc='white', ec='none', pad=1))
    ax.set_title(rf'$J={J:g},\ h={h:g}$', fontsize=7)
axes[0].set_ylabel(r'$D(\rho_{+}^{S},\rho_{-}^{S})$')
axes[0].legend(fontsize=5.5, frameon=False, loc='center', bbox_to_anchor=(0.5, 0.45))
fig.tight_layout(pad=0.3, w_pad=0.5)
fig.savefig('figure2.pdf')
