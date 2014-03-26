
include datasets.mk

GTAIL_SIG	= 14
MIN_QUAL	= 30
ERR_TOL		= 1

A			= $(ALN_DIR)/$(D)
GTAIL_SAM	= $(A)/gtail_aln/gtail_reads.sam
NVTR_SAM	= $(A)/nvtr_aln/nvtr_reads.sam
PPK			= $(LOG)/$(D)_parsed.pk
PARSE_REPORT	= $(LOG)/$(D)_parse.pdf
PARSE_LOG		= $(LOG)/$(D)_parse.log
ALIGN_REPORT = $(LOG)/$(D)_align.pdf

align:
	@echo Aligning $(D)
	@[ -d $(A) ] || (cd $(PATSY_DIR); $(ALIGN) -p $(CORES) -f $(REF) -1 $($(D)_FQ1) -2 $($(D)_FQ2) -o $(A) -s $(LOG)/$(D)_aln.pk \
	-r $(ALIGN_REPORT) > $(LOG)/$(D)_aln.log)

parse: align
	@echo Parsing $(D)
	(cd $(PATSY_DIR); $(PARSE) -d $(D) -g $(GTAIL_SAM) -n $(NVTR_SAM) -f $(REF) -l $(GTAIL_SIG) -e $(ERR_TOL) -i $(PARSED_TRANSCRIPTS) -q $(MIN_QUAL) \
	-o $(PPK) -t -r $(PARSE_REPORT) > $(PARSE_LOG) )

