.PHONY: com ta tp tt

com:
	git ci -a

# Test patsy-align:
ta:
	@rm -fr ./test/t_aln;./patsy-align -1 test/t1.fq -2 test/t2.fq -f test/t.fas -o test/t_aln

# Test patsy-parse:
tp:
	@./patsy-parse -t -d name -g test/t_aln/gtail_aln/gtail_reads.sam -n test/t_aln/nvtr_aln/nvtr_reads.sam -f test/t.fas
# Test patsy-test:
tt:
	@./patsy-test -na WT -nb MUT -a test/WT1A_parsed.pk -b test/MUT1A_parsed.pk

