| Parameters                      | Variables                    | Value | Units |
|:--------------------------------|:-----------------------------|------:|:-----:|
| Electrode                       | Anode                        | False |       |
| Physics parameters              | Temperature                  |    80 |   C   |
|                                 | Initial potential            |   0.1 |   V   |
|                                 | Final potential              |   0.9 |   V   |
|                                 | Step potential               | 0.001 |   V   |
| Continuous Stirred-Tank Reactor | Concentration = f(potential) | False |       |
|                                 | Catalyst Active surface area |   1.0 |  cm2  |
|                                 | Volumetric flux              |   1.0 |  L/s  |
| Rate Constants                  | --------------------------   |  ---- | ----  |
| Pre-exponential                 | A                            |     1 |       |
|                                 | j*                           |  True |       |
|                                 | j* (value)                   |  1000 | A/cm2 |
| Transition state theory         | kappa * k_B * T^m / h        | False |       |
|                                 | kappa                        |     1 |       |
|                                 | m                            |     1 |       |
| Experimental                    | k_f, k_b                     | False |       |
| Chemical                        | Ga, DG_reaction              |  True |       |
|                                 | DG_reaction                  | False |       |
|                                 | G_formation                  |  True |       |
