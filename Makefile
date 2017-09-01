run:
	./preprocess.r
	./analysis.r oos
	./analysis.r cv

clean:
	rm -rvf *.csv  *~ *.txt data/* plots/*
