# COVID19 during an emergency evacuation

The results in this study [1] were obtained in two steps:

## Run an updated version of the social force model

Run different scenarios using an updated version of PanicPackage code from [2]. Within this: 
 1. Create a subdirectory called 'dump'
 2. update the sd_lib.c file according to the diff file sd_lib-diff.txt and compile the application
 3. update the sd.par file to assign different parameter values
 4. parse_files.sh to have multiple runs store datafiles into the 'dump' folder

## Run post_analysis_sd.m which has functions that must be run in the order they are listed. 

% coarse_sd3(); % 1 coarse the dataset so that we sample at 0.01 s
% get_evac_times(); % 2 extract evacuation times for 90% of the crowd
% replay_a_trial(); % for verification
% plot_traj(); % for verification
% calculate_distmat(); % 3 this takes time! go for a coffee break!
% calculate_exposure_from_IA(16); % 4 takes the tau value as input
composite_plot(); % 5 draw the composite plot for paper

## References

1. Butail S and Porfiri M (2021) The Effect of An Emergency Evacuation on the Spread of COVID19.  *Front. Phys.* 8:631264. 
2. Helbing D, Farkas I, Vicsek T. Simulating dynamical features of escape panic. *Nature* (2000) 407:487–90. http://angel.elte.hu/panic/
