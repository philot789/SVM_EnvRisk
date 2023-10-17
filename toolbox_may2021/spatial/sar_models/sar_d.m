% PURPOSE: An example of using sar() 
%          spatial autoregressive model
%          (on a small data set)                  
%---------------------------------------------------
% USAGE: sar_d (see also sar_d2 for a large data set)
%---------------------------------------------------

clear all;

load anselin.ford;

W = anselin;
% W-matrix from Anselin's neigbhorhood crime data set

load anselin.dat; % standardized 1st-order spatial weight matrix
% 5 columns:
% column1 = crime
% column2 = household income
% column3 = house values
% column4 = latitude coordinate
% column5 = longitude coordinate

n = length(anselin);
y = anselin(:,1);
x = [ones(n,1) log(anselin(:,2:3))];
vnames = strvcat('crime','constant','income','hvalue');


info.lflag = 0; % use full lndet no approximation
result0 = sar(y,x,W,info);
prt(result0,vnames);


