Introduction
============

IMMBoost is a novel strategy to improve the performance in two classification tasks: 1) to distinguish CRMs in a particular expression domain from non-functional sequences; 2) to discriminate whether a CRM drives expression in one expression domain versus other domains.

Installation
============

1. License
-----------
::

Note that LIBSVM and LIBLINEAR are both used in this code, license related issue please see COPYRIGHT files. Also, source code IMM from SCRMshaw: http://veda.cs.uiuc.edu/SCRMshaw/ is used.

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
	
	Several PERL modules and R packages are needed. Please make sure to have them installed before you run the code:

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
					+-- 
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

	2) CRM_vs_bkg/ (same file structure as above)

	


2. Data Format
===========

Input files including:

1) 
2) 







Invocation
==========

Running the following command will the available functions::

	$ IMMBoost -h

Gives::


License
============