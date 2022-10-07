.PHONY: all clean

DRUG ?= $(shell bash -c 'read -p "Enter drug name:" drug; echo $$drug')
PROTOCOL ?= $(shell bash -c 'read -p "Enter protocol name:" protocol; echo $$protocol')
DRUG_LIST = dofetilide bepridil terfenadine cisapride verapamil ranolazine
DRUG_LIST2 = quinidine sotalol chlorpromazine ondansetron diltiazem mexiletine

# first run for optimisation
# with check plot
first_run:
	cd modelling/scripts/; \
	python3 binding_kinetics_comparison.py $(DRUG) $(PROTOCOL) True True True

# steady state vs transient phase
state_dependent_all_drugs: 
	cd modelling/scripts/; \
	for drug in DRUG_LIST; do \
		python3 binding_kinetics_comparison.py $$drug $(PROTOCOL) True; \
		python3 binding_kinetics_comparison.py $$drug $(PROTOCOL) False; \
	done

state_dependent: steady transient

steady:
	# ifeq (check, $(filter check, $(MAKECMDGOALS)))
	# 	cd modelling/scripts/; \
	# 	python3 binding_kinetics_comparison.py $(DRUG) $(PROTOCOL) True False True
	# else
	cd modelling/scripts/; \
	python3 binding_kinetics_comparison.py $(DRUG) $(PROTOCOL) True
	# endif

transient:
	cd modelling/scripts/; \
	python3 binding_kinetics_comparison.py $(DRUG) $(PROTOCOL) False

# protocols
protocol_dependent:
	cd modelling/scripts/; \
	for prot in Milnes Pneg80 P0 P40; do \
		python3 binding_kinetics_comparison.py $(DRUG) $$prot True; \
		python3 binding_kinetics_comparison.py $(DRUG) $$prot False; \
	done

sensitivity_analysis:
	cd modelling/scripts/; \
	for drug in $(DRUG_LIST2); do \
		echo $$drug; \
		python3 sensitivity_analysis.py $$drug; \
		python3 fig_script/sensitivity_plot.py $$drug; \
	done

sensitivity_analysis_plot:
	cd modelling/scripts/; \
	for drug in $(DRUG_LIST2); do \
		echo $$drug; \
		python3 fig_script/sensitivity_plot.py $$drug; \
	done

sensitivity_analysis_plot_short:
	cd modelling/scripts/; \
	for drug in $(DRUG_LIST); do \
		echo $$drug; \
		python3 fig_script/sensitivity_plot.py $$drug True; \
	done

create_dir:
	cd figures/binding_kinetics_comparison; \
	for drug in $(DRUG_LIST2); do \
		mkdir $$drug/; \
		cd $$drug/; \
		mkdir Milnes/; \
		cd ..; \
	done

scripts:
	cd modelling/scripts/; \
	python3 drug_norm_check.py; \
	python3 SA_parameter.py; \
	python3 SA_drug.py; \