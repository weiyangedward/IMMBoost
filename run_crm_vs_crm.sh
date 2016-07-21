mkdir -p sampleOutput
mkdir -p sampleOutput/CRM_vs_CRM
perl script/IMMBoost.pl --task crm_vs_crm --step 345678 --times 2 sampleData/CRMsetsList.txt sampleOutput/CRM_vs_CRM/ sampleData/CRMsets/ sampleData/CRMGroup.txt
