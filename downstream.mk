
# Filter transcripts by G-tail read coverage.
MINCOV_WT				= 1000
MINCOV_MUT				= 1
GTAIL_COV_FILTER_REP			= $(LOG)/gtail_cov_filter.pdf
GCF_TEST1				= $(LOG)/TEST_WT1_vs_MUT1_gcf.tab
GCF_TEST2				= $(LOG)/TEST_WT2_vs_MUT2_gcf.tab
gtail_cov_filter:
	@echo Filtering test results by G-tail read coverage:
	@scripts/01-gtail_cov_filter.R --mc_wt $(MINCOV_WT) --mc_mut $(MINCOV_MUT) --test1 $(TEST1_TAB_TRS) --test2 $(TEST2_TAB_TRS) --test1_out $(GCF_TEST1) --test2_out $(GCF_TEST2) --rep $(GTAIL_COV_FILTER_REP)

TCOR_WT1	= $(LOG)/tech_cor_WT1.pdf
TCOR_WT2	= $(LOG)/tech_cor_WT2.pdf
TCOR_MUT1	= $(LOG)/tech_cor_MUT1.pdf
TCOR_MUT2	= $(LOG)/tech_cor_MUT2.pdf

gtail_tech_corr:
	@echo Plotting correlation between technical replicates:	
	@scripts/011-gtail_tech_corr.R $(TCOR_WT1) $(TR_WT1) 
	@scripts/011-gtail_tech_corr.R $(TCOR_WT2) $(TR_WT2) 
	@scripts/011-gtail_tech_corr.R $(TCOR_MUT1) $(TR_MUT1) 
	@scripts/011-gtail_tech_corr.R $(TCOR_MUT2) $(TR_MUT2) 

# Classify tail run distributions:
CLS_MINCOV	= 1000

CLS_REP_WT1	= $(LOG)/CLS_WT1.pdf
CLS_REP_WT2	= $(LOG)/CLS_WT2.pdf

CLS_REP_MUT1	= $(LOG)/CLS_MUT1.pdf
CLS_REP_MUT2	= $(LOG)/CLS_MUT2.pdf

classify_tail_dists:
	@scripts/02-classify_tail_dist.R WT1 $(CLS_MINCOV) $(CLS_REP_WT1) $(TR_WT1)
	@scripts/02-classify_tail_dist.R WT2 $(CLS_MINCOV) $(CLS_REP_WT2) $(TR_WT2)
	@scripts/02-classify_tail_dist.R MUT1 $(CLS_MINCOV) $(CLS_REP_MUT1) $(TR_MUT1)
	@scripts/02-classify_tail_dist.R MUT2 $(CLS_MINCOV) $(CLS_REP_MUT2) $(TR_MUT2)

PAL_TOTAL	= $(BASE)/ExtData/PAL_Cerevisiae_total.tab
PAL_CYT		= $(BASE)/ExtData/PAL_Cerevisiae.tab
PASTA_DATA	= $(BASE)/ExtData/PASTA_polyAranking.tab

COR_PAL_TOTAL_WT1	= $(LOG)/PAL_total_vs_WT1.pdf
COR_PAL_TOTAL_WT2	= $(LOG)/PAL_total_vs_WT2.pdf
COR_PAL_CYT_WT1		= $(LOG)/PAL_cyt_vs_WT1.pdf
COR_PAL_CYT_WT2		= $(LOG)/PAL_cyt_vs_WT2.pdf
COR_PAL_CYT_PASTA	= $(LOG)/PAL_cyt_vs_PASTA.pdf
COR_PAL_TOTAL_PASTA	= $(LOG)/PAL_total_vs_PASTA.pdf
COR_WT1_PASTA		= $(LOG)/WT1_vs_PASTA.pdf
COR_WT2_PASTA		= $(LOG)/WT2_vs_PASTA.pdf

corr_with_studies:
	@scripts/03-corr_with_study.R --d1 $(GCF_TEST1) --name1 PAT_WT1 --df1 runs_ma --tf1 transcripts --d2 "ExtData/PASTA_polyAranking.tab" --name2 PASTA --df2 weighted_length --tf2 transcripts --rep test1.pdf
	@scripts/03-corr_with_study.R --d1 $(GCF_TEST2) --name1 PAT_WT2 --df1 runs_ma --tf1 transcripts --d2 "ExtData/PASTA_polyAranking.tab" --name2 PASTA --df2 weighted_length --tf2 transcripts --rep test2.pdf
	@scripts/03-corr_with_study.R --d1 $(GCF_TEST1) --name1 PAT_WT1 --df1 runs_ma --tf1 transcripts --d2 $(COR_PAL_CYT_PASTA) --name2 PASTA --df2 Mean.TL --tf2 Transcript.ID --rep $(COR_PAL_TOTAL_WT1)
	@scripts/03-corr_with_study.R --d1 $(GCF_TEST2) --name1 PAT_WT2 --df1 runs_ma --tf1 transcripts --d2 $(PAL_TOTAL) --name2 PAL_total --df2 Mean.TL --tf2 Transcript.ID --rep $(COR_PAL_TOTAL_WT2)
	@scripts/03-corr_with_study.R --d1 $(GCF_TEST1) --name1 PAT_WT1 --df1 runs_ma --tf1 transcripts --d2 $(PAL_CYT) --name2 PAL_cyt --df2 Mean.TL --tf2 Transcript.ID --rep $(COR_PAL_CYT_WT1)
	@scripts/03-corr_with_study.R --d1 $(GCF_TEST2) --name1 PAT_WT2 --df1 runs_ma --tf1 transcripts --d2 $(PAL_CYT) --name2 PAL_cyt --df2 Mean.TL --tf2 Transcript.ID --rep $(COR_PAL_CYT_WT2)
	@scripts/03-corr_with_study.R --d1 $(PASTA_DATA) --name1 PASTA_weighted_length --df1 weighted_length --tf1 transcripts --d2 $(PAL_CYT) --name2 PAL_cyt --df2 Mean.TL --tf2 Transcript.ID --rep $(COR_PAL_CYT_PASTA)
	@scripts/03-corr_with_study.R --d1 $(PASTA_DATA) --name1 PASTA_weighted_length --df1 weighted_length --tf1 transcripts --d2 $(PAL_TOTAL) --name2 PAL_total --df2 Mean.TL --tf2 Transcript.ID --rep $(COR_PAL_TOTAL_PASTA)
	@scripts/03-corr_with_study.R --d1 $(PASTA_DATA) --name1 PASTA_weighted_length --df1 weighted_length --tf1 transcripts --d2 $(GCF_TEST1) --name2 PAT_WT1 --df2 runs_ma --tf2 transcripts --rep $(COR_WT1_PASTA)
	@scripts/03-corr_with_study.R --d1 $(PASTA_DATA) --name1 PASTA_weighted_length --df1 weighted_length --tf1 transcripts --d2 $(GCF_TEST2) --name2 PAT_WT2 --df2 runs_ma --tf2 transcripts --rep $(COR_WT2_PASTA)

TECH_DESIGN_FILE=technical_replicates.tab
BIO_DESIGN_FILE=biological_replicates.tab
SUMMARY_REPORT=Replicate_correlations.pdf

corr_summary:
	@(cd Log; ../scripts/04-corr_summary.R --tech_replicates $(TECH_DESIGN_FILE) --bio_replicates $(BIO_DESIGN_FILE) --rep $(SUMMARY_REPORT))
