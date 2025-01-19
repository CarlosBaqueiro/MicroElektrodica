| Parameters                      | Variables                    | Value | Units |
|:--------------------------------|:-----------------------------|------:|:-----:|
| Electrode                       | Anode                        |  True |       |
| Physics parameters              | Temperature                  |    25 |   C   |
|                                 | Initial potential            |   0.0 |   V   |
|                                 | Final potential              |   1.0 |   V   |
|                                 | Step potential               | 0.001 |   V   |
| Continuous Stirred-Tank Reactor | Concentration = f(potential) | False |       |
|                                 | Catalyst Active surface area |   1.0 |  cm2  |
|                                 | Volumetric flux              |   1.0 |  L/s  |
| Rate Constants                  | --------------------------   |  ---- | ----  |
| Pre-exponential                 | A                            |     1 |       |
|                                 | j*                           | False |       |
|                                 | j* (value)                   |   500 | A/cm2 |
| Transition state theory         | kappa * k_B * T^m / h        | False |       |
|                                 | kappa                        |     1 |       |
|                                 | m                            |     1 |       |
| Experimental                    | k_f, k_b                     |  True |       |
| Chemical                        | Ga, DG_reaction              | False |       |
|                                 | DG_reaction                  | False |       |
|                                 | G_formation                  | False |       |

 