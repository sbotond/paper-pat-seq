A pipeline to analyse PAT-seq data
==================================

Reference
---------

- TODO

Using the pipeline
------------------

The pipeline can be used by invoking the following `make` targets:

- Fetch raw data from ENA: `make fetch`
- Generate transcriptome from SGD annotation: `make transcriptome`
- Align and parse reads: `make parse`
- Test for differential polyadenylation: `make test` or `make lsf_test`
- Parse spike-in reads (`make parse_spikeins`) and build "error model" (`make error_model`)
- Filter test results by G-tail coverage: `make gtail_cov_filter`
- Plot correlation between technical replicates: `make gtail_tech_corr`
- Plot and cluster tail length distributions: `make classify_tail_dists`
- Correlate thresholded tail lengths with PASTA and PAL-seq: `make corr_with_studies`

Index of selected raw results
-----------------------------

### Quantifying tail length slippage using spike-in standards

- Tail run lengths until the first 1-5 non-A bases in reads mapped to spike-in poly(A) tracts [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/error_thr.pdf)

#### Alignment

|   Sample  |   Alignment log    |      Alignment report       |
|:---------:|:------------------:|:---------------------------:|
| MUT1A     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1A_aln.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1A_align.pdf) |
| MUT1B     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1B_aln.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1B_align.pdf) |
| MUT1C     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1C_aln.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1C_align.pdf) |
| MUT1D     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1D_aln.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1D_align.pdf) |
| MUT2A     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2A_aln.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2A_align.pdf) |
| MUT2B     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2B_aln.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2B_align.pdf) |
| MUT2C     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2C_aln.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2C_align.pdf) |
| MUT2D     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2D_aln.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2D_align.pdf) |

#### Parsing alignments

|   Sample  |   Parse log    |     Parse report       |
|:---------:|:------------------:|:---------------------------:|
| MUT1A     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1A_parse.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1A_parse.pdf) |
| MUT1B     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1B_parse.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1B_parse.pdf) |
| MUT1C     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1C_parse.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1C_parse.pdf) |
| MUT1D     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1D_parse.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT1D_parse.pdf) |
| MUT2A     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2A_parse.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2A_parse.pdf) |
| MUT2B     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2B_parse.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2B_parse.pdf) |
| MUT2C     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2C_parse.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2C_parse.pdf) |
| MUT2D     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2D_parse.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/MUT2D_parse.pdf) |

#### Testing differences between wild type and mutant tail runs

|   Comparison  |   Test log    |     Test report       |   Results |
|:---------:|:------------------:|:---------------------------:|:------------------------------:|
| WT1 vs. MUT1     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/TEST_WT1_vs_MUT1.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/TEST_WT1_vs_MUT1.pdf) | [CSV](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/TEST_WT1_vs_MUT1_trs.tab) |
| WT2 vs. MUT2     |  [TXT](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/TEST_WT2_vs_MUT2.log) | [PDF](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/TEST_WT2_vs_MUT2.pdf) | [CSV](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/TEST_WT2_vs_MUT2_trs.tab) |

#### Tail run distributions from all transcripts with G-tail coverage > 1000

- [WT1](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/CLS_WT1.pdf)
- [WT2](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/CLS_WT2.pdf)
- [MUT1](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/CLS_MUT1.pdf)
- [MUT2](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/CLS_MUT2.pdf)

#### Cross-study correlation

- [PAL_total vs. WT1](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/PAL_total_vs_WT1.pdf)
- [PAL_total vs. WT2](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/PAL_total_vs_WT2.pdf)
- [PAL_total vs. PASTA](http://www.ebi.ac.uk/goldman-srv/pat-seq/Log/PAL_total_vs_PASTA.pdf)

Dependencies
------------

- [Platform LSF](http://en.wikipedia.org/wiki/Platform_LSF)
- [Python](http://www.python.org/) 2.x
- [numpy](https://pypi.python.org/pypi/numpy) >= 1.6.2
- [matplotlib](https://pypi.python.org/pypi/matplotlib) >= 1.1.0
- [scipy](https://pypi.python.org/pypi/scipy) >= 0.10.1
- [biopython](https://pypi.python.org/pypi/biopython) >= 1.60
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) >= 2.1.0
- [samtools](http://samtools.sourceforge.net) >= 0.1.19+
- [wget](https://www.gnu.org/software/wget/)

Using the analysis tools
-----------------------

The analysis tool can be found under `patsy/`:

### *patsy-align* - classify read pairs and align them using Bowtie2

```
usage: patsy-align [-h] -1 fq1 -2 fq2 -f ref [-o outdir] [-s stats_pickle]
                   [-l gtail_sig] [-G gtag_min] [-N max_N] [-I min_fsize]
                   [-X max_fsize] [-p nr_threads] [-r report]

Align PAT-seq reads using Bowtie2 (1.1).

optional arguments:
  -h, --help       show this help message and exit
  -1 fq1           First FASTQ file.
  -2 fq2           Second FASTQ file.
  -f ref           Reference fasta.
  -o outdir        Output directory.
  -s stats_pickle  Stats pickle file.
  -l gtail_sig     Portion of read start/end used for G-tail classification
                   (14).
  -G gtag_min      Minimum G-tag length(3).
  -N max_N         Maximum number of Ns in the first -l bases (6).
  -I min_fsize     Minimum fragment size (0).
  -X max_fsize     Maximum fragment size (500).
  -p nr_threads    Number of threads to use (1).
  -r report        Report PDF.
```

### *patsy-parse* - parse classified and aligned PAT-seq read pairs

```
usage: patsy-parse [-h] -g gtail_sam -n nvtr_sam -d dataset_id -f ref
                   [-l gtail_sig] [-G gtag_min] [-N max_N] [-e err_tol]
                   [-o out_pickle] [-i tr_list] [-q min_q] [-r report] [-t]

Parse classified and aligned PAT-seq read pairs (1.1).

optional arguments:
  -h, --help     show this help message and exit
  -g gtail_sam   SAM file containing G-tail alignments.
  -n nvtr_sam    SAM file containing NVTR alignments.
  -d dataset_id  Dataset identifier.
  -f ref         Reference fasta.
  -l gtail_sig   Portion of read start/end used for G-tail classification.
  -G gtag_min    Minimum G-tag length(3).
  -N max_N       Maximum number of Ns in the first -l bases (6).
  -e err_tol     Number of errors tolerated in the tail.
  -o out_pickle  Output pickle file.
  -i tr_list     List of transcripts considered.
  -q min_q       Mapping quality treshold (30).
  -r report      Report PDF.
  -t             Plot per-transcript coverage reports.
```

### *patsy-test* - test for differential polyadenylation in PAT-seq data

```
usage: patsy-test [-h] -a [a_pickles [a_pickles ...]] -na a_name -b
                  [b_pickles [b_pickles ...]] -nb b_name [-i lrt_list]
                  [-P lik_penalty] [-M min_size_U] [-s sig_level]
                  [-op out_pickle] [-ot out_trs] [-og out_glob]
                  [-otr out_runs_prefix] [-orr out_rep_prefix] [-r report]
                  [-t]

Test for differential polyadenylation in PAT-seq data (1.1).

optional arguments:
  -h, --help            show this help message and exit
  -a [a_pickles [a_pickles ...]]
                        Parsed read pickles - group A.
  -na a_name            Name of data group A.
  -b [b_pickles [b_pickles ...]]
                        Parsed read pickles - group B.
  -nb b_name            Name of data group B.
  -i lrt_list           Transcripts to be tested with anchors LRT.
  -P lik_penalty        Log-likelihood penalty for data points outside valid
                        range.
  -M min_size_U         Minimum sample size when performing Mann-Whitney U
                        test (30).
  -s sig_level          Significance level.
  -op out_pickle        Output pickle file.
  -ot out_trs           Output tabular file: transcript properties.
  -og out_glob          Output tabular file: global results.
  -otr out_runs_prefix  Output tabular file: tail runs prefix.
  -orr out_rep_prefix   Output tabular file: tail means per replicate.
  -r report             Report PDF.
  -t                    Plot reports for all transcripts.
```

### *patsy-spike* - estimate the number of sequencing errors in runs of bases.

```
usage: patsy-spike [-h] -n spike_sam -f ref [-w window] [-m max_errors_plot]
                   [-o out_pickle] [-q min_q] [-r report] [-l read_len]
                   [-pk pickle]

Estimate the number of sequencing errors in runs of bases (1.0).

optional arguments:
  -h, --help          show this help message and exit
  -n spike_sam        SAM file containing NVTR alignments.
  -f ref              Reference fasta with the spike-in sequences.
  -w window           Size of the flanking sequence around the run of As.
  -m max_errors_plot  Maximum number of errors for which to plot the length
                      distribution.
  -o out_pickle       Output pickle file.
  -q min_q            Mapping quality treshold (30).
  -r report           Report PDF.
  -l read_len         Read length.
  -pk pickle          Result pickle file.
```

