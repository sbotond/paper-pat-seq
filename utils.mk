
.PHONY: com


com:
	git ci -a

BASE		=$(shell pwd)
LOG			=$(BASE)/Log
ALN_DIR		=$(BASE)/Aln
ANN_DIR		=$(BASE)/Ann

CORES		=4
MEMORY		=20000
SUBMIT_JOB	=bsub -M $(MEMORY) -n $(CORES) 

PATSY_DIR	=$(BASE)/patsy
ALIGN		=$(PATSY_DIR)/patsy-align
PARSE		=$(PATSY_DIR)/patsy-parse
TEST		=$(PATSY_DIR)/patsy-test
BUILD_REF	=$(PATSY_DIR)/utils/get_ref_transcripts


