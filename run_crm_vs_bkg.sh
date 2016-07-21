mkdir -p sampleOutput
mkdir -p sampleOutput/CRM_vs_bkg
perl script/IMMBoost.pl --task crm_vs_bkg --step 12345678 --ktimes 2 --nfolds 5 sampleData/CRMsetsList.txt sampleOutput/CRM_vs_bkg/ sampleData/CRMsets/ sampleData/CRMGroup.txt
