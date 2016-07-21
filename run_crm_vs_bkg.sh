mkdir -p sampleOutput
mkdir -p sampleOutput/CRM_vs_bkg
perl script/IMMBoost.pl --task crm_vs_bkg --step 12345678 --times 2 sampleData/CRMsetsList.txt sampleOutput/CRM_vs_bkg/ sampleData/CRMsets/ sampleData/CRMGroup.txt
