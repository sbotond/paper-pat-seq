include utils.mk

GENOME					=saccharomyces_cerevisiae
REF_BASE				=$(ALN_DIR)/reference
REF         				=$(REF_BASE)/$(GENOME).fas

PARSED_TRANSCRIPTS			=$(REF_BASE)/genes_single_iso.tab	# Parsed transcripts (coding with no intron)
SPS_TRANSCRIPTS				=$(REF_BASE)/genes_single_iso_peak.tab	# Transcripts with a single poly(A) site
SPIKEINS				=$(BASE)/ExpInfo/spikeins.fas

DD					=$(BASE)/RawData

# Wild type datasets:

WT1DS		= WT1A WT1B WT1C WT1D
WT2DS		= WT2A WT2B WT2C WT2D
WTDS		= $(WT1DS) $(WT2DS)

# WT Replicate #1

WT1A_FQ1	= $(DD)/lane1polyAGWT1a_1_sequence.txt.gz
WT1A_FQ2	= $(DD)/lane1polyAGWT1a_2_sequence.txt.gz

WT1B_FQ1	= $(DD)/lane1polyAGWT1b_1_sequence.txt.gz
WT1B_FQ2	= $(DD)/lane1polyAGWT1b_2_sequence.txt.gz

WT1C_FQ1	= $(DD)/lane1polyAGWT1c_1_sequence.txt.gz
WT1C_FQ2	= $(DD)/lane1polyAGWT1c_2_sequence.txt.gz

WT1D_FQ1	= $(DD)/lane1polyAGWT1d_1_sequence.txt.gz
WT1D_FQ2	= $(DD)/lane1polyAGWT1d_2_sequence.txt.gz

# WT Replicate #2

WT2A_FQ1	= $(DD)/lane1polyAGWT2a_1_sequence.txt.gz
WT2A_FQ2	= $(DD)/lane1polyAGWT2a_2_sequence.txt.gz

WT2B_FQ1	= $(DD)/lane1polyAGWT2b_1_sequence.txt.gz
WT2B_FQ2	= $(DD)/lane1polyAGWT2b_2_sequence.txt.gz

WT2C_FQ1	= $(DD)/lane1polyAGWT2c_1_sequence.txt.gz
WT2C_FQ2	= $(DD)/lane1polyAGWT2c_2_sequence.txt.gz

WT2D_FQ1	= $(DD)/lane1polyAGWT2d_1_sequence.txt.gz
WT2D_FQ2	= $(DD)/lane1polyAGWT2d_2_sequence.txt.gz


# Mutant datasets:

MUT1DS		= MUT1A MUT1B MUT1C MUT1D
MUT2DS		= MUT2A MUT2B MUT2C MUT2D
MUTDS		= $(MUT1DS) $(MUT2DS)

# MUT replicate #1

MUT1A_FQ1	= $(DD)/lane1polyAGMUT1a_1_sequence.txt.gz 
MUT1A_FQ2	= $(DD)/lane1polyAGMUT1a_2_sequence.txt.gz

MUT1B_FQ1	= $(DD)/lane1polyAGMUT1b_1_sequence.txt.gz
MUT1B_FQ2	= $(DD)/lane1polyAGMUT1b_2_sequence.txt.gz

MUT1C_FQ1	= $(DD)/lane1polyAGMUT1c_1_sequence.txt.gz
MUT1C_FQ2	= $(DD)/lane1polyAGMUT1c_2_sequence.txt.gz

MUT1D_FQ1	= $(DD)/lane1polyAGMUT1d_1_sequence.txt.gz
MUT1D_FQ2	= $(DD)/lane1polyAGMUT1d_2_sequence.txt.gz

# MUT replicate #2

MUT2A_FQ1	= $(DD)/lane1polyAGMUT2a_1_sequence.txt.gz
MUT2A_FQ2	= $(DD)/lane1polyAGMUT2a_2_sequence.txt.gz

MUT2B_FQ1	= $(DD)/lane1polyAGMUT2b_1_sequence.txt.gz
MUT2B_FQ2	= $(DD)/lane1polyAGMUT2b_2_sequence.txt.gz

MUT2C_FQ1	= $(DD)/lane1polyAGMUT2c_1_sequence.txt.gz
MUT2C_FQ2	= $(DD)/lane1polyAGMUT2c_2_sequence.txt.gz

MUT2D_FQ1	= $(DD)/lane1polyAGMUT2d_1_sequence.txt.gz
MUT2D_FQ2	= $(DD)/lane1polyAGMUT2d_2_sequence.txt.gz

# All datasets:

DATASETS	= $(WTDS) $(MUTDS)

