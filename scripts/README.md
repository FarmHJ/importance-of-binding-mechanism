# Organisation of the scripts

The `fig_script` folder contains all scripts used to plot the figures in the paper.
The `checking_script` folder contains all scripts used for data analysis but not shown in the results.

## Main results
[background_trapping.py](./background_trapping.py) - Generates data to introduce the trapping mechanism.
[binding_kinetics_comparison.py](./binding_kinetics_comparison.py) - Calibrates the ionic conductance of the CS model from the SD model for a dofetilide-like drug and a verapamil-like drug, then compare the APD90s at steady state.
[AP_simulation.py](./AP_simulation.py) - Simulates action potentials of the ORd-SD model and the ORd-CS model at transient phase.
[protocol_dependence.py](./protocol_dependence.py) - Compares the Hill curve of the SD model for drugs with different protocols and the APD90 of the ORd-CS model when the ionic conductance is scaled with the Hill curves.
[SA_param_space.py](./SA_param_space.py) - Explore the parameter space of drug-related parameters (Vhalf, Kmax and Ku) and compute the APD90 differences between the ORd-SD model and the ORd-CS model for a given virtual drug.
[SA_curve.py](./SA_curve.py) - Compute the APD90 differences between the ORd-SD model and the ORd-CS model for the parameter space around the boundary surface where the APD90s are similar.
[SA_drugs.py](./SA_drugs.py) - Compute the APD90 differences between the two AP models for all synthetic drugs.
[combine_APD.py](./combine_APD.py) - Combine all simulated data of the parameter space with essential information for easy loading when plotting figures.

## Supplementary materials
[supp_comparison_drugs.py](./supp_comparison_drugs.py) - Compare the APD90 of the ORd-SD model and the ORd-CS model for each
synthetic drug and generate the Hill curves of the drug for different protocols.
[supp_SA_parameter.py](./supp_SA_parameter.py) - Compute the RMSD between APD90s of the ORd-SD model and the ORd-CS model when the Hill coefficient of each synthetic drug changes.
