include utils.mk

GENOME					=saccharomyces_cerevisiae
REF_BASE				=$(ALN_DIR)/reference
REF         				=$(REF_BASE)/$(GENOME).fas

PARSED_TRANSCRIPTS			=$(REF_BASE)/genes_single_iso.tab	# Parsed transcripts (coding with no intron)
SPS_TRANSCRIPTS				=$(REF_BASE)/genes_single_iso_peak.tab	# Transcripts with a single poly(A) site
SPIKEINS				=$(BASE)/ExpInfo/spikeins.fas

DD					=$(BASE)/RawData


fetch:
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472943/ERR472943_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472943/ERR472943_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472944/ERR472944_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472944/ERR472944_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472935/ERR472935_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472935/ERR472935_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472945/ERR472945_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472945/ERR472945_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472940/ERR472940_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472940/ERR472940_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472933/ERR472933_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472933/ERR472933_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472934/ERR472934_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472934/ERR472934_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472937/ERR472937_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472937/ERR472937_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472936/ERR472936_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472936/ERR472936_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472938/ERR472938_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472938/ERR472938_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472941/ERR472941_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472941/ERR472941_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472939/ERR472939_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472939/ERR472939_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472946/ERR472946_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472946/ERR472946_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472948/ERR472948_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472948/ERR472948_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472942/ERR472942_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472942/ERR472942_2.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472947/ERR472947_1.fastq.gz;
	wget -P RawData  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR472/ERR472947/ERR472947_2.fastq.gz;


MUT1A_FQ       =   $(DD)/ERR472943_1.fastq.gz
MUT1A_FQ       =   $(DD)/ERR472943_2.fastq.gz
MUT1B_FQ       =   $(DD)/ERR472944_1.fastq.gz
MUT1B_FQ       =   $(DD)/ERR472944_2.fastq.gz
MUT1C_FQ       =   $(DD)/ERR472935_1.fastq.gz
MUT1C_FQ       =   $(DD)/ERR472935_2.fastq.gz
MUT1D_FQ       =   $(DD)/ERR472945_1.fastq.gz
MUT1D_FQ       =   $(DD)/ERR472945_2.fastq.gz
MUT2A_FQ       =   $(DD)/ERR472940_1.fastq.gz
MUT2A_FQ       =   $(DD)/ERR472940_2.fastq.gz
MUT2B_FQ       =   $(DD)/ERR472933_1.fastq.gz
MUT2B_FQ       =   $(DD)/ERR472933_2.fastq.gz
MUT2C_FQ       =   $(DD)/ERR472934_1.fastq.gz
MUT2C_FQ       =   $(DD)/ERR472934_2.fastq.gz
MUT2D_FQ       =   $(DD)/ERR472937_1.fastq.gz
MUT2D_FQ       =   $(DD)/ERR472937_2.fastq.gz
WT1A_FQ       =   $(DD)/ERR472936_1.fastq.gz
WT1A_FQ       =   $(DD)/ERR472936_2.fastq.gz
WT1B_FQ       =   $(DD)/ERR472938_1.fastq.gz
WT1B_FQ       =   $(DD)/ERR472938_2.fastq.gz
WT1C_FQ       =   $(DD)/ERR472941_1.fastq.gz
WT1C_FQ       =   $(DD)/ERR472941_2.fastq.gz
WT1D_FQ       =   $(DD)/ERR472939_1.fastq.gz
WT1D_FQ       =   $(DD)/ERR472939_2.fastq.gz
WT2A_FQ       =   $(DD)/ERR472946_1.fastq.gz
WT2A_FQ       =   $(DD)/ERR472946_2.fastq.gz
WT2B_FQ       =   $(DD)/ERR472948_1.fastq.gz
WT2B_FQ       =   $(DD)/ERR472948_2.fastq.gz
WT2C_FQ       =   $(DD)/ERR472942_1.fastq.gz
WT2C_FQ       =   $(DD)/ERR472942_2.fastq.gz
WT2D_FQ       =   $(DD)/ERR472947_1.fastq.gz
WT2D_FQ       =   $(DD)/ERR472947_2.fastq.gz

# Wild type datasets:

WT1DS		= WT1A WT1B WT1C WT1D
WT2DS		= WT2A WT2B WT2C WT2D
WTDS		= $(WT1DS) $(WT2DS)

# Mutant datasets:

MUT1DS		= MUT1A MUT1B MUT1C MUT1D
MUT2DS		= MUT2A MUT2B MUT2C MUT2D
MUTDS		= $(MUT1DS) $(MUT2DS)

# All datasets:

DATASETS	= $(WTDS) $(MUTDS)

