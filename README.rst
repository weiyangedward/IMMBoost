Introduction
============

IMMBoost is a novel strategy to improve the 
performance in two classification tasks: 1) to 
distinguish CRMs in a particular expression domain 
from non-functional sequences; 2) to discriminate 
CRMs drive expression in one expression 
domain versus the others.

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

	(Please un-tar the source code in a directory where the 
	path to this directory does not contain empty space. 
	Otherwise there will be errors while running the code.)

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

	Examples of input as well as command to run 
	the code are provided for both of CRM vs CRM and CRM vs 
	Bkg. Output can be found in folder 'sampleOutput'. Note
	that although these two pipelines have very similar names
	among scripts, the implementation can be very different.
	Please don't simply copy over or attempt to combine 
	scripts that have the same names in directory CRM_vs_CRM 
	and CRM_vs_bkg.(Please empty this folder before your new 
	runs, although the code will do this for you if it finds
	the direcotry is pre-existed.):

	1) CRM vs CRM:
	>> sh run_crm_vs_crm.sh

	2) CRM vs Bkg:
	>> sh run_crm_vs_bkg.sh
	
	Descriptions of directory and files are shown within 
	parentheses (Example has 10 runs of 5-fold cross validations.):
	
	1) CRM_vs_CRM/
	+
	+-- CRMset1/ (output directory of CRM set1)
		+
		+-- IMM.average.auc (average AUC from IMM prediction)
		+-- IMM_SVM.average.auc (average AUC from IMMBoost-SVM prediction)
		+-- IMM_RF.average.auc (average AUC from IMMBoost-RF prediction)
		+-- ensembleModel.average.auc (average AUC from IMMBoost-Ensemble prediction)
		+-- kmerSVM.average.auc (average AUC from kmer-SVM prediction)
		+-- allData/ (CRM and negative seq fasta files and 
			the trained IMM models using this data)
		+-- time1/ (1st 5-fold cross validation)
			+
			+-- fold1/ (results of 1st fold)
			...
			+-- fold5/
		...
		+-- time10/
	+-- CRMset2/
	...
	+-- summaryAUC_msIMMBaseline.txt (AUCs from IMM)
	+-- summaryAUC_IMM_RF.txt (AUCs from IMMBoost-RF)
	+-- summaryAUC_IMM_SVM.txt (AUCs from IMMBoost-SVM)
	+-- summaryAUC_ensembleModel.txt (AUCs from IMMBoost-Ensemble)
	+-- summaryAUC_kmerSVM.txt (AUCs from kmer-SVM)

	2) CRM_vs_bkg/ (same file structure as above)

	Detailed performance please see summaryAUC*txt 
	files. Each of these files corresponds to the 
	average AUC scores over 10trials x 5folds cross 
	validation using one model. Each file has two 
	columns, where the first column has CRMset names, 
	and the second column has the average AUC scores. 
	Note that since sampleData is just a random subset 
	of real data, and therefore the performance in 
	sampleOutput might not be ideal.


2. Data Format
--------------
::
	
	Input files including:

	1. CRMsetsList.txt : a list of path to CRMsets. Each 
	CRMset folder should have a sub-folder called "fasta", 
	inside which there are files: 

		1) CRM.fasta : CRM seq file. Each seq has a unique name.

		2) randomGenomicSeq.fasta : Random genomic seq.

		3) msCRM.fasta : msCRM seq file (if you don't have 
			msCRM seq file, you can copy CRM.fasta over and 
			change the seqID to a corresponding species_seqID, 
			e.g., Dmel_seqID).

		4) negCRM.fasta : Negative CRM seq for CRM vs CRM task.

		5) negmsCRM.fasta : Negative msCRM seq for CRM vs CRM task.

	2. sampleData/CRMsets/ : a directory of input data. 
		Each sub directory contains a CRMset files.
		When running this code on your own data set,
		please replace CRM sets in sampleData with your
		own CRM sets.


	Code to generate background sequences with similar
	length and GC content to query CRMs are included in
	'./tools/generate_background_seq/'. This code requires 
	input file such as preprocessed accessible genome 
	regions with exons and repeats masked. Dmel's accessible 
	genome with exon masked can be obtained from UCSC. To 
	mask repeat we used Tandem Repeat Finder: http://tandem.
	bu.edu/trf/trf.download.html. To run the code on a toy 
	example:

	>cd tools/generate_background_seq/
	>sh run.sh

3. Data Set
-----------
::

	Supplementary data (Data S1 and S2) mentioned in the main text of 
	IMMBoost paper is also inlcuded in this repository:

	Data S1: 38 CRM sets. These are the input data used to train and evaluate our 
	classifiers. 36 out of 38 CRM sets are borrowed from (Kazemian et al., 2014), 
	and we added two additional CRM sets: mapping3.adult and mapping3.larva. Each 
	CRM set contains CRM sequences of D. melanogaster at a specific expression 
	domain (CRM.fasta), multi-species CRM sequences of D. melanogaster and other 10 (
	D.sim, D.sec, D.yak, D.ere, D.ana, D.pse, D.per, D.moj, D.vir, D.gri) Drosophila 
	species at a specific expression domain (msCRM.fasta), random genomic sequences (
	randomGenomicSeq.fasta), CRMs sequences of D. melanogaster in other expression 
	domains (negCRM.fasta), and multi-species CRM sequences in other expression 
	domains (negmsCRM.fasta). 

	Data S2: Accessible inter-genic genomic regions at stage 5 of D. melanogaster 
	embryonic development. This data set was used to generate random genomic 
	sequences as negative training sets while training classifiers to distinguish 
	CRMs from random genomic regions. DNA-accessibility information was taken from 
	DNase I Hypersensitivity (DHS) peaks at stage 5 of D. melanogaster development. 
	We also masked repeats in these regions using Tandem Repeats Finder.



4. To Run
---------
::

	perl IMMBoost.pl [options] CRMList Outdir Datadir

     --task <str>      What task to perform? default=crm_vs_crm. There are two modes:
                          1) "--task crm_vs_bkg": classify CRM from background 
                                  genomic seq; 
                          2) "--task crm_vs_crm": classify CRM from other CRM seq
  
  	--step <str>      What steps to run? default=12345678.
	                        1. prepare data for n-fold cross validation
	                        2. generate IMM score feature
	                        3. IMM prediction
	                        4. IMM-SVM prediction
	                        5. IMM-RF prediction 
	                        6. generate kmer-SVM feature
	                        7. kmer-SVM prediction
	                        8. IMM-Ensemble prediction
  	--nfolds <int>    To perform n-fold cross validation. default=5.
  	--ktimes <int>    To repeat n-fold cross validation for k times. default=2.


Additional Information
======================
All questions please contact author Wei Yang throgh email: 
weiyang4 AT illinois DOT edu

