%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%   ADSHAM v.1.0  Main Program
%%        created by Sooncheol Hwang on 2022. 04. 05.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%
%%   READ ME First !!!!!!
%%
%%  This is a main script to initiate ADSHAM by setting up prescribed test cases.
%%  To get the simulation started, you'll have to enter the case name 
%%  of tests into the filename. Two test cases are available for now as shown below. 
%%
%%  1) 'Convergence_n***' : convergence test to investigate numerical order of accuracy
%%     In convergence case, 6 different cell numbers (or cell sizes) 
%%     including reference cell number (n=12,800) are to be simulated.
%%     For example, for the test with 800 cells, just enter 
%%     the filename as 'Convergence_n800'.
%%
%%  2) 'PureDiffusion' : A pure scalar diffusion of Gaussian distribution on 2D space
%%     Please input the filename 'PureDiffusion' for this simulation, then the test 
%%     outputs with figures will be produced.


clear all; close all; clc;

filename = 'Convergence_n200';
run('Engine.m')
