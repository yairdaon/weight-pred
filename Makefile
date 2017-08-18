run:
	./process.r
	./lag_data.r

clean:
	rm -rvf *.csv  *~ *.txt data/*
