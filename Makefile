test:
	make plots/exp_time_series_test.svg
	make plots/precision_weighted_predictions_test.pdf
	make plots/var_chl_test.pdf


oos:
	make plots/exp_time_series_oos.svg
	make plots/precision_weighted_predictions_oos.pdf
	make plots/var_chl_oos.pdf

loocv:
	make plots/exp_time_series_loocv.svg
	make plots/precision_weighted_predictions_loocv.pdf
	make plots/var_chl_loocv.pdf

all_pred:
	make 
	make analysis/all_predictions.r

everything:
	make all_pred
	make test
	make oos
	make loocv

clean:
	rm -rvf *.csv *.txt *~ analysis/*~ process/*~  helpers/*~  plots/*  

## Run code to make the plots
plots/precision_weighted_predictions_test.pdf: analysis/predictors.r data/test_full_run.Rdata
	./analysis/predictors.r 
plots/precision_weighted_predictions_oos.pdf: analysis/predictors.r data/oos_full_run.Rdata
	./analysis/predictors.r oos
plots/precision_weighted_predictions_loocv.pdf: analysis/predictors.r data/loocv_full_run.Rdata
	./analysis/predictors.r loocv

plots/exp_time_series_test.svg: analysis/predictors.r data/test_full_run.Rdata
	./analysis/predictors.r 
plots/exp_time_series_oos.svg: analysis/predictors.r data/oos_full_run.Rdata
	./analysis/predictors.r oos 
plots/exp_time_series_loocv.svg: analysis/predictors.r data/loocv_full_run.Rdata
	./analysis/predictors.r loocv


plots/var_chl_test.pdf: analysis/vars.r data/test_full_run.Rdata
	./analysis/vars.r
plots/var_chl_oos.pdf: analysis/vars.r data/oos_full_run.Rdata
	./analysis/vars.r oos
plots/var_chl_loocv.pdf: analysis/vars.r data/loocv_full_run.Rdata
	./analysis/vars.r loocv



analysis/all_predictions.r: data/processed_block.Rdata data/oos_serial_days.Rdata FORCE
	./analysis/all_predictions.r


## Make the dates files
data/test_serial_days.Rdata:	process/dates.r
	./process/dates.r
data/oos_serial_days.Rdata:	process/dates.r
	./process/dates.r
data/loocv_serial_days.Rdata:	process/dates.r
	./process/dates.r

## Pre-process
data/processed_block.Rdata: process/preprocess.r
	./process/preprocess.r

## Run a bunch of s-maps.
data/test_full_run.Rdata: process/calcs.r data/processed_block.Rdata data/test_serial_days.Rdata
	./process/calcs.r
data/oos_full_run.Rdata: process/calcs.r data/processed_block.Rdata data/oos_serial_days.Rdata
	./process/calcs.r oos
data/loocv_full_run.Rdata: process/calcs.r data/processed_block.Rdata data/loocv_serial_days.Rdata
	./process/calcs.r loocv



FORCE: ;
