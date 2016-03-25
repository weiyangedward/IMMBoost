Introduction
============

IMMBoost is a novel strategy to improve the 
performance in two classification tasks: 1) to 
distinguish CRMs in a particular expression domain 
from non-functional sequences; 2) to discriminate 
whether a CRM drives expression in one expression 
domain versus other domains.

Installation
============

1. License
-----------
::

Note that LIBSVM and LIBLINEAR are both used in 
this code, license related issue please see 
COPYRIGHT files. Also, source code IMM from SCRMshaw: http://veda.cs.uiuc.edu/SCRMshaw/ is used.

2. Compile source code
--------------------------
::

	>> cd src
		
	1) IMM:
	>> cd imm
	>> make clean
	>> make

	2) liblinear:
	>> cd liblinear-2.1
	>> make clean
	>> make
	>> cd python
	>> make

	3) libsvm:
	>> cd libsvm-3.21
	>> make clean
	>> make
	>> cd python
	>> make

3. Required packages and libraries
----------------------------------
::
	
	Several PERL modules and R packages are needed. 
	Please make sure to have them installed before you 
	run the code:

	1) BioPerl (PERL)
	2) randomForest (R)
	3) ROCR (R)

QuickStart
==========

1. Example
----------
::

	Example data files for both of CRMvsCRM and CRMvsBkg are provided:

	1) CRMvsCRM:
	>> sh run_crm_vs_crm.sh

	2) CRMvsBkg:
	>> sh run_crm_vs_crm.sh
	
	Output can be found in folder sampleOutput:
	
	1) CRM_vs_CRM/
			+
			+-- CRMset1/
					+
					+-- IMM.average.auc
					+-- IMM_SVM.average.auc
					+-- IMM_RF.average.auc
					+-- ensembleModel.average.auc
					+-- kmerSVM.average.auc
					+-- allData/
					+-- time1/
							+
							+-- fold1/
							+-- fold2/
							...
					+-- time2/
					...
			+-- CRMset2/
			...
			+-- summaryAUC_msIMMBaseline.txt
			+-- summaryAUC_IMM_RF.txt
			+-- summaryAUC_IMM_SVM.txt
			+-- summaryAUC_ensembleModel.txt
			+-- summaryAUC_kmerSVM.txt

	2) CRM_vs_bkg/ (same file structure as above)

	Detailed performance please see summaryAUC*txt 
	files. Each of these files corresponds to the 
	average AUC scores over 10trials x 5folds cross 
	validation. Each file has two columns, where the 
	first column has CRMset names, and the second 
	column has the average AUC scores. Note that since 
	sampleData is just a random subset of real data, 
	and therefore the performance in sampleOutput 
	might not be ideal.


2. Data Format
--------------
::
	
	Input files including:

	1) "CRMsetsList.txt" : a list of path to CRMsets. Each CRMset folder should have 
	2) "../sampleOutput/CRM_vs_CRM/" : a directory for output files
	3) "../sampleData/CRMsets/" : a directory for data input. Each sub directory should be a CRMset, which contains 
	4) "CRM.group.V3.txt" : a list of grouping of CRMsets, where each row is a group. How to define a "group" would be subjective to users or biological grouptruth in our case.







Invocation
==========

Running the following command will the available functions::

	$ IMMBoost -h

Gives::


License
============