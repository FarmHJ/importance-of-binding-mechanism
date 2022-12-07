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
