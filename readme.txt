The analysis is run as follows:

1) Run different scenarios using a slightly updated version of PanicPackage code from [1].  Within this: 
	a) update the sd_lib.c file and compile the application
	b) update the sd.par file to assign different parameter values, and then use the 	parse_files.sh to have multiple runs store datafiles into the 'dump' folder

2) Run post_analysis_sd.m which has functions that must be run in the order they are listed. 

% coarse_sd3(); % 1 coarse the dataset so that we sample at 0.01 s
% get_evac_times(); % 2 extract evacuation times for 90% of the crowd
% replay_a_trial(); % for verification
% plot_traj(); % for verification
% calculate_distmat(); % 3 this takes time! go for a coffee break!
% calculate_exposure_from_IA(16); % 4 takes the tau value as input
composite_plot(); % 5 draw the composite plot for paper


[1] http://angel.elte.hu/panic/