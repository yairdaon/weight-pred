run:
	make clean
	./predictors.r oos
	./predictors.r cv

big:
	make cleaner
	./preprocess.r
	./calcs.r oos 
	./predictors.r oos
	./calcs.r cv
	./predictors.r cv

clean:
	rm -rvf *.csv  *~ *.txt plots/*
cleaner:
	rm -rvf *.csv  *~ *.txt data/* plots/*
