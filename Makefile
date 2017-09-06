run:
	./preprocess.r
	./weights.r oos
	./weights.r cv

clean:
	rm -rvf *.csv  *~ *.txt data/* plots/*
