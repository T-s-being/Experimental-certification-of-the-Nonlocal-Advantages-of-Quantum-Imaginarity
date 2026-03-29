# Experimental certification of the Nonlocal Advantages of Quantum Imaginarity
Matlab project relevant to the work reported in the preprint ["Experimental certification of the Nonlocal Advantages of Quantum Imaginarity"]()

The Matlab function [`compute_NAQI.m`](compute_NAQI.m) computes the quantum parameter for the NAQI criterion for an arbitrary two‑qubit state $\rho_{AB}$.

The project [`MainProject1.m`](MainProject1.m) generates a theoretical figure, in which the red line corresponds to the violation parameter $\Delta_{\mathrm{QI}}$ for $\rho_{\mathrm{B}}(p)$ with varying mixture parameter $p \in [0,1]$, and the dots represent our experimental choices. After the operation is completed, the optimal angles for the experimental dots are stored as rows in the variable `Experiment_x`. More precisely, the i‑th row `Experiment_x(i,:)` contains the optimal angle parameters $\lbrace \theta_1,\phi_1,\theta_2,\phi_2,\theta_3,\phi_3,\theta_4,\phi_4\rbrace$ (in radians) for the i‑th dot.

The project [`MainProject2.m`](MainProject2.m) generates a theoretical figure for $\rho_{\mathrm{W}}(p)$, which is analogous to `MainProject1`.

It is worth noting that, due to the existence of multiple optimal solutions, the results obtained by running the program may differ from those given in the paper, yet they can all achieve optimal outcomes.

Noting that because of the optimal solution 
## Dependencies
- MATLAB (tested on R2022b)
