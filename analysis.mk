include datasets.mk
include downstream.mk

.PHONY: parse test lsf_test

transcriptome:
	@python scripts/filter_transcripts.py -r $(ANN_DIR)/S288C_reference_sequence_R63-1-1_20100105.fsa \
	-t $(ANN_DIR)/saccharomyces_cerevisiae_R63-1-1_20100109.gff -p $(ANN_DIR)/3Tfill_peaks.gff \
	-o $(ANN_DIR)/transcripts_filtered.fa -pk $(ANN_DIR)/peaks.pk -ot $(REF) -o $(PARSED_TRANSCRIPTS) \
	-os $(SPS_TRANSCRIPTS)
	@cat $(SPIKEINS) >> $(REF)

parse:
	@$(foreach DS,$(DATASETS), $(SUBMIT_JOB) "make parse -f align_parse.mk D=$(DS)";) # wondering if there is a nicer way to do this

# Parsed pickles:
WT1PKS	= $(LOG)/WT1A_parsed.pk $(LOG)/WT1B_parsed.pk $(LOG)/WT1C_parsed.pk $(LOG)/WT1D_parsed.pk
WT2PKS	= $(LOG)/WT2A_parsed.pk $(LOG)/WT2B_parsed.pk $(LOG)/WT2C_parsed.pk $(LOG)/WT2D_parsed.pk

MUT1PKS	= $(LOG)/MUT1A_parsed.pk $(LOG)/MUT1B_parsed.pk $(LOG)/MUT1C_parsed.pk $(LOG)/MUT1D_parsed.pk
MUT2PKS	= $(LOG)/MUT2A_parsed.pk $(LOG)/MUT2B_parsed.pk $(LOG)/MUT2C_parsed.pk $(LOG)/MUT2D_parsed.pk

# Minimum sample size for Mann-Whitney U tests:
MW_MIN	= 20
# Significance level:
ALPHA	= 0.05
# Test output files:
TEST1_PK		= $(LOG)/TEST_WT1_vs_MUT1.pk
TEST2_PK		= $(LOG)/TEST_WT2_vs_MUT2.pk

TEST1_TAB_TRS	= $(LOG)/TEST_WT1_vs_MUT1_trs.tab
TEST2_TAB_TRS	= $(LOG)/TEST_WT2_vs_MUT2_trs.tab

TEST1_TAB_GLOB		= $(LOG)/TEST_WT1_vs_MUT1_glob.tab
TEST2_TAB_GLOB		= $(LOG)/TEST_WT2_vs_MUT2_glob.tab

TEST1_TAB_TRUNS		= $(LOG)/TEST_WT1_vs_MUT1_truns_
TEST2_TAB_TRUNS		= $(LOG)/TEST_WT2_vs_MUT2_truns_

TEST1_TAB_TMEANS		= $(LOG)/TEST_WT1_vs_MUT1_tmeans.tab
TEST2_TAB_TMEANS		= $(LOG)/TEST_WT2_vs_MUT2_tmeans.tab

TEST1_REP		= $(LOG)/TEST_WT1_vs_MUT1.pdf
TEST2_REP		= $(LOG)/TEST_WT2_vs_MUT2.pdf

TEST1_LOG		= $(LOG)/TEST_WT1_vs_MUT1.log
TEST2_LOG		= $(LOG)/TEST_WT2_vs_MUT2.log

# Tail distribution tabular files:
TR_MUT1A		=	$(LOG)/TEST_WT1_vs_MUT1_truns_MUT1A.tab
TR_MUT1B		=	$(LOG)/TEST_WT1_vs_MUT1_truns_MUT1B.tab
TR_MUT1C		=	$(LOG)/TEST_WT1_vs_MUT1_truns_MUT1C.tab
TR_MUT1D		=	$(LOG)/TEST_WT1_vs_MUT1_truns_MUT1D.tab
TR_MUT1			= 	$(TR_MUT1A) $(TR_MUT1B) $(TR_MUT1C) $(TR_MUT1D)

TR_WT1A			=	$(LOG)/TEST_WT1_vs_MUT1_truns_WT1A.tab
TR_WT1B			=	$(LOG)/TEST_WT1_vs_MUT1_truns_WT1B.tab
TR_WT1C			=	$(LOG)/TEST_WT1_vs_MUT1_truns_WT1C.tab
TR_WT1D			=	$(LOG)/TEST_WT1_vs_MUT1_truns_WT1D.tab
TR_WT1			=	$(TR_WT1A) $(TR_WT1B) $(TR_WT1C) $(TR_WT1D)

TR_MUT2A		=	$(LOG)/TEST_WT2_vs_MUT2_truns_MUT2A.tab
TR_MUT2B		=	$(LOG)/TEST_WT2_vs_MUT2_truns_MUT2B.tab
TR_MUT2C		=	$(LOG)/TEST_WT2_vs_MUT2_truns_MUT2C.tab
TR_MUT2D		=	$(LOG)/TEST_WT2_vs_MUT2_truns_MUT2D.tab
TR_MUT2			=	$(TR_MUT2A) $(TR_MUT2B) $(TR_MUT2C) $(TR_MUT2D)

TR_WT2A			=	$(LOG)/TEST_WT2_vs_MUT2_truns_WT2A.tab
TR_WT2B			=	$(LOG)/TEST_WT2_vs_MUT2_truns_WT2B.tab
TR_WT2C			=	$(LOG)/TEST_WT2_vs_MUT2_truns_WT2C.tab
TR_WT2D			=	$(LOG)/TEST_WT2_vs_MUT2_truns_WT2D.tab
TR_WT2			=	$(TR_WT2A) $(TR_WT2B) $(TR_WT2C) $(TR_WT2D)

# The SAM file with spikein reads
SPIKEIN_NVTR = nvtr_aln/nvtr_reads_sort.bam
SPIKEIN_NAMES = pGIBS-LYS pGIBS-PHE pGIBS-THR
SPIKEIN_CORE = $(ALN_DIR)/spikeins
SPIKEIN_SAM = $(SPIKEIN_CORE)_sort.sam
ERROR_THR = $(LOG)/error_thr.txt
ERROR_THR_PDF = $(LOG)/error_thr.pdf

test:
	@echo testing WT1 agains MUT1:
	@(cd patsy; ./patsy-test -t -a $(WT1PKS) -na WT1 -b $(MUT1PKS) -nb MUT1 -M $(MW_MIN) -s $(ALPHA) -op $(TEST1_PK) -ot $(TEST1_TAB_TRS) -og $(TEST1_TAB_GLOB) -otr $(TEST1_TAB_TRUNS) -orr $(TEST1_TAB_TMEANS) -r $(TEST1_REP) > $(TEST1_LOG) )
	@echo testing WT2 agains MUT2:
	@(cd patsy; ./patsy-test -t -a $(WT2PKS) -na WT2 -b $(MUT2PKS) -nb MUT2 -M $(MW_MIN) -s $(ALPHA) -op $(TEST2_PK) -ot $(TEST2_TAB_TRS) -og $(TEST2_TAB_GLOB) -otr $(TEST2_TAB_TRUNS) -orr $(TEST2_TAB_TMEANS) -r $(TEST2_REP) > $(TEST2_LOG) )

JOB_DEP_CMD = echo $(DATASETS) | sed 's/ /"\) \&\& done("/g'
JOB_DEP = -w 'done("$(shell $(JOB_DEP_CMD))")'

$(SPIKEIN_SAM):
	@echo Parsing out the spikein reads.
	@$(foreach DS,$(DATASETS), $(SUBMIT_JOB) -J $(DS) "samtools view -b -h $(ALN_DIR)/$(DS)/$(SPIKEIN_NVTR) $(SPIKEIN_NAMES) -o $(ALN_DIR)/$(DS)/spikeins.bam";)
	$(SUBMIT_JOB) -J header $(JOB_DEP) "samtools view -H $(ALN_DIR)/`echo $(DATASETS) | cut -d\  -f 1`/spikeins.bam > $(SPIKEIN_CORE)_header.sam; \
	samtools merge -h $(SPIKEIN_CORE)_header.sam -f $(SPIKEIN_CORE).bam $(ALN_DIR)/*/spikeins.bam; \
	samtools sort $(SPIKEIN_CORE).bam $(SPIKEIN_CORE)_sort_tmp; \
	samtools index $(SPIKEIN_CORE)_sort_tmp.bam; \
	samtools sort -n $(SPIKEIN_CORE)_sort_tmp.bam $(SPIKEIN_CORE)_sort; \
	samtools view $(SPIKEIN_CORE)_sort.bam > $(SPIKEIN_SAM)"

parse_spikeins: $(SPIKEIN_SAM)

error_model: $(SPIKEIN_SAM)
	@(cd patsy; $(SUBMIT_JOB) "python patsy-spike -f $(SPIKEINS) -n $(SPIKEIN_SAM) -r $(ERROR_THR_PDF) > $(ERROR_THR)")

lsf_test:
	@$(SUBMIT_JOB) "make test"
