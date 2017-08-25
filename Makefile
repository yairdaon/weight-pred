run:
	./preprocess.r
	./analysis.r

clean:
	rm -rvf *.csv  *~ *.txt data/*
