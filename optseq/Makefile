# Generation of optimal sequences for the simplecomposition
# Time-stamp: <2013-01-12 17:24 christophe@pallier.org>


.PHONY: all10 all8


# generate lists of 10 conditions of 8 stims, min SOA=8

all10: # 10 conditions
	rm -f simpcomp*
	optseq2 --ntp 512 --tr 1.5 --o simcomp \
		--psdwin 0 20 .5 --nkeep 10 --nsearch 1000  --sumdelays \
		--ev cond1 6 8 --ev cond2 6 8 --ev cond3 6 8 --ev cond4 6 8 \
		--ev cond5 6 8 --ev cond6 6 8 --ev cond7 6 8 --ev cond8 6 8 \
		--ev cond9 6 8 --ev cond10 6 8  \
		--evc -2 -1 0  1 2 -2 -1 0 1 2

# generate lists of 8 conditions of 8 stims, min SOA=8

all8: # 8 conditions
	rm -f simpcomp*
	optseq2 --ntp 410 --tr 1.5 --o simcomp \
		--psdwin 0 20 .5 --nkeep 10 --nsearch 1000  --sumdelays \
		--ev cond1 6 8 --ev cond2 6 8 --ev cond3 6 8 --ev cond4 6 8 \
		--ev cond5 6 8 --ev cond6 6 8 --ev cond7 6 8 --ev cond8 6 8 \
		--evc -3 -1 1 3 -3 -1 1 3
#	Rscript stats_seq.R # generate 'durations_of_null_events.pdf'
