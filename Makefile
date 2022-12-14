.PHONY: all clean

DRUG ?= $(shell bash -c 'read -p "Enter drug name:" drug; echo $$drug')
PROTOCOL ?= $(shell bash -c 'read -p "Enter protocol name:" protocol; echo $$protocol')
DRUG_LIST = dofetilide bepridil terfenadine cisapride verapamil ranolazine
DRUG_LIST2 = quinidine sotalol chlorpromazine ondansetron diltiazem mexiletine

# steady state vs transient phase
state_dependent_all_drugs: 
	cd scripts/; \
	for drug in $(DRUG_LIST); do \
		echo $$drug; \
		python3 supp_comparison_drugs.py $$drug; \
		cd fig_script/; \
		python3 supp_drugs_comparison.py $$drug; \
		cd ..; \
	done

# Figure 1 (background figure)
Figure1:
	cd scripts/; \
	python3 background_trapping.py; \
	cd fig_script/; \
	python3 background.py; \

# Figure 3 (APD90 comparison for dofetilide-like drug)
Figure3:
	cd scripts/; \
	python3 binding_kinetics_comparison.py dofetilide; \
	cd fig_script/; \
	python3 AP_compare_drug.py dofetilide; \

# Figure 4 (APD90 comparison for verapamil-like drug)
Figure4:
	cd scripts/; \
	python3 binding_kinetics_comparison.py verapamil; \
	cd fig_script/; \
	python3 AP_compare_drug.py verapamil; \

# Figure 5 (AP comparison at initial phases)
Figure5:
	cd scripts/; \
	python3 AP_simulation.py; \
	cd fig_script/; \
	python3 AP_transient.py; \

# Figure 6 (Hill curves of the drug using different protocols)
Figure6:
	cd scripts/; \
	python3 protocol_dependence.py; \
	cd fig_script/; \
	python3 protocol.py; \

