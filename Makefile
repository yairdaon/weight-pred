now:
	make leftover

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

leftover:
	make data/leftover_predictions_bloom_loocv.Rdata
	./analysis/sum.r bloom
	make data/leftover_predictions_full_loocv.Rdata
	./analysis/sum.r full

all_pred:
	make analysis/all_predictions.r
	make analysis/predictors.r	

everything:
	make all_pred
	make test
	make oos
	make loocv
	make leftover

vars: ./analysis/vars.r data/avg_var_df_loocv.Rdata data/processed_block.Rdata
	./analysis/vars.r loocv

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

data/leftover_predictions_full_loocv.Rdata: process/leftover.r data/processed_block.Rdata
	./process/leftover.r 

data/leftover_predictions_bloom_loocv.Rdata: process/leftover.r data/processed_block.Rdata
	./process/leftover.r bloom

## Pre-process
data/processed_block.Rdata: process/preprocess.r data/solar_radiation.Rdata
	./process/preprocess.r

## Run a bunch of s-maps.
data/test_full_run.Rdata: process/calcs.r data/processed_block.Rdata data/test_serial_days.Rdata
	./process/calcs.r
data/oos_full_run.Rdata: process/calcs.r data/processed_block.Rdata data/oos_serial_days.Rdata
	./process/calcs.r oos
data/loocv_full_run.Rdata: process/calcs.r data/processed_block.Rdata data/loocv_serial_days.Rdata
	./process/calcs.r loocv

data/solar_radiation.Rdata: process/solar.r
	./process/solar.r


FORCE: ;
