# Central qubit in an Ising spin bath: non-Markovianity and partial information plots

Code accompanying the manuscript

> W. S. Martins and D. O. Soares-Pinto, *Decoherence function as a probe of non-Markovianity and environmental redundancy in an Ising spin bath* (submitted to Physics Letters A).

It computes the decoherence function, the BLP trace distance, and the exact qubit–fragment mutual information (partial information plot, PIP) for a central qubit dephasing under a periodic Ising chain in a Gibbs state, and reproduces Figs. 2, 4 and 5 of the manuscript.

## Model and units

- Bath: $H_B = -J\sum_i \sigma^z_i\sigma^z_{i+1} - h\sum_i \sigma^z_i$, periodic boundary conditions, Gibbs state at inverse temperature $\beta$.
- Coupling: $H_{SB} = \alpha\,\sigma^z\otimes\sum_i\sigma^z_i$ (pure dephasing).
- Units: $J = 1$, so $h$, $\alpha$ and $T = 1/\beta$ are in units of $J$ and $t$ is in units of $\hbar/J$; $k_B = \hbar = 1$. For $J = 0$ the same numerical values are kept.
- The decoherence function has period $\tau = \pi/(2\alpha)$.

## Files

| File | Content |
|---|---|
| `decoherence_model.py` | Shared physics. Decoherence function for $J>0$ (transfer matrix) and $J=0$; dephased-qubit entropy $\mathcal{S}_2$; closed-form PIP for $J=0$ (Eq. 34); boundary-conditioned decoherence functions (Eq. 32) and exact PIP of a contiguous fragment for any $J$ (Eq. 33); zero-field correlation length. |
| `plot_figure2.py` | Fig. 2: trace distance $D = \lvert\Gamma(t)\rvert$ for $(J,h) = (1,1), (1,0.1), (1,0), (0,1)$ and $\beta = 0, 0.75, 1.75, 4$. |
| `plot_figure4.py` | Fig. 4: $J=0$. (a) $\lvert\Gamma(t)\rvert$; (b) normalized PIP $\mathcal{I}(\mathcal{S}:\mathcal{F})/S(\rho_S)$ at $t=\tau/2$. |
| `plot_figure5.py` | Fig. 5: $J=1$, $h=0$. (a) normalized PIP at $t_1=\pi/(4\alpha N)$; (b) single-spin value versus $\xi/N$ for $N = 6, 10, 20, 40$. |
| `style.py` | Shared matplotlib style. |
| `verify_brute_force.py` | Checks Eqs. (33)–(34) against direct construction and diagonalization of $\rho_{\mathcal{SF}}$ for $N=8$, and the bound $\mathcal{I}/S(\rho_S)\le f$ for $J=0$ (Eq. 35). |

All parameters are set at the top of each plotting script and can be edited directly.

## Usage

```bash
pip install -r requirements.txt
python verify_brute_force.py      # consistency checks (prints max deviation)
python plot_figure2.py            # writes figure2.pdf
python plot_figure4.py            # writes figure4.pdf
python plot_figure5.py            # writes figure5.pdf
```

Each script runs in a few seconds on a laptop. No data files are needed: all figures are computed from the closed-form expressions.

## Requirements

Python ≥ 3.9, NumPy, Matplotlib (tested with Python 3.13, NumPy 2.5, Matplotlib 3.11).

## Citation

If you use this code, please cite the manuscript above and this Zenodo record.

## Contact

Wilson S. Martins — W.SantanaMartins@umb.edu
