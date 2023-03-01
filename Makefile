.PHONY: all clean

DRUG_LIST = bepridil terfenadine cisapride ranolazine quinidine sotalol chlorpromazine ondansetron diltiazem mexiletine

# Figure 1 (background figure)
Figure1:
	cd scripts/; \
	python3 background_trapping.py; \
	cd fig_script/; \
	python3 background.py; \

# Figure 2 (method description figure)
# Require data generated from the script binding_kinetics_comparison.py for verapamil
Figure2:
	cd scripts/fig_script/; \
	python3 method_description.py; \

# Figure 3 (APD90 and qNet comparison for example drug T (dofetilide))
Figure3:
	cd scripts/; \
	python3 binding_kinetics_comparison.py dofetilide; \
	cd fig_script/; \
	python3 AP_compare_drug.py dofetilide; \

# Figure 4 (APD90 and qNet comparison for example drug N (verapamil))
Figure4:
	cd scripts/; \
	python3 binding_kinetics_comparison.py verapamil; \
	cd fig_script/; \
	python3 AP_compare_drug.py verapamil; \

# Figure 5 (AP and APD90 comparison at transient phases)
Figure5:
	cd scripts/; \
	python3 AP_transient_simulation.py dofetilide; \
	python3 AP_transient_simulation.py verapamil; \
	cd fig_script/; \
	python3 AP_transient.py; \

# Figure 6 (Hill curves of the drug using different protocols)
Figure6:
	cd scripts/; \
	python3 protocol_dependence.py; \
	cd fig_script/; \
	python3 protocol_dependence.py; \

# Figure 7 (parameter exploration) and Figure S12
Figure7:
	cd scripts/; \
	python3 SA_param_space.py; \
	python3 SA_drugs.py; \
	python3 SA_curve.py; \
	cd fig_script/; \
	python3 SA_3D.py; \

# Figure S1 to S10 (APD90, qNet and Hill curve for several drugs)
FigureS1_S10:
	cd scripts/; \
	python3 supp_comparison_drugs.py; \
	cd fig_script/; \
	for drug in $(DRUG_LIST); do \
		python3 supp_drugs_comparison.py $$drug; \
	done

# Figure S11 (RMSD distribution with varying Hill coefficient)
FigureS11:
	cd scripts/; \
	python3 supp_SA_N.py; \
	cd fig_script/; \
	python3 supp_N_difference.py; \

# Figure S13 and S14 (RMSD and APD90 at different Vhalf-trap values)
# Require data generated from the script SA_param_space.py
FigureS13_S14:
	cd scripts/fig_script/; \
	python3 APD_Vhalf.py; \
