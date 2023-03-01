# Importance of modelling hERG binding

In this work, we study the impact of including state-dependent drug binding in a mathematical ion channel model of the human Ether-À-Go-Go-Related Gene (hERG) channel.
We compared the action potential predictions when modelling drug binding of hERG using a state-dependent model versus a conductance scaling model.

Source code associated with "Importance of modelling hERG binding in predicting drug-induced action potential prolongations for drug safety assessment" by Hui Jia Farm,
Michael Clerx,
Fergus Cooper,
Liudmila Polonchuk,
Ken Wang,
David J. Gavaghan and
Chon Lok Lei.

## Code
All scripts that generate the data and plot the figures are in the [scripts](./scripts) folder.

## Reproduce figures, including intermediate files
Figures can be reproduced with the command:
```bash
make FigureX
```
where `X` is the figure number. The list of commands to reproduce the figures are provided below.

### List of commands to reproduce the figures
1. `Figure1` 
    - to reproduce Figure 1 (background figure)
2. `Figure2` 
    - to reproduce Figure 2 (method description figure)
    - requires data generated from the script [binding_kinetics_comparison.py](scripts/binding_kinetics_comparison.py) for example drug N (also run in command `Figure4`)
3. `Figure3`
    - to reproduce Figure 3 (model comparison for example drug T)
4. `Figure4`
    - to reproduce Figure 4 (model comparison for example drug N)
5. `Figure5`
    - to reproduce Figure 5 (AP comparison at transient phase)
6. `Figure6`
    - to reproduce Figure 6 (protocol dependence of AP prolongation)
7. `Figure7`
    - to reproduce Figure 7 and Figure S12 (parameter space exploration)
8. `FigureS1_S10`
    - to reproduce Figure S1 to S10 (APD90, qNet and Hill curves for several drugs)
9. `FigureS11`
    - to reproduce Figure S11 (RMSD distribution with varying parameter N)
10. `FigureS13_S14`
    - to reproduce Figure S13 and S14 (RMSD and APD90 at different Vhalf-trap values)
    - requires data generated from the script [SA_param_space.py](scripts/SA_param_space.py) (also run in command `Figure7`)

Some figures were modified for better visualisation, i.e. Figure 1 to Figure 4.

## Acknowledging this work
If you publish any work based on the contents of this repository please cite ([CITATION file](CITATION)):

Farm, H. J., Clerx, M., Cooper, F., Polonchuk, L., Wang, K., Gavaghan, D. J. and Lei, C. L. (2023).
Importance of modelling hERG binding in predicting drug-induced action potentials for drug safety assessment.
_Frontiers in Pharmacology_, 14:1110555.
[doi:10.3389/fphar.2023.1110555](https://doi.org/10.3389/fphar.2023.1110555).

