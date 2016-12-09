%Examples for how to run the functions described in the main text
clear all %clean up work space

%%%READ ME!
%Unless otherwise noted, this example file, and all included code, are
%covered under the GPL, version 3, available here:
%(http://www.gnu.org/licenses/gpl-3.0.en.html)

%One function needed to run this code (mvnpdf_log.m) cannot be distributed
%because it is a modified version of a proprietary MATLAB function.
%This function can be recreated by taking the standard MATLAB "mvnpdf.m"
%function, and removing "exp" from the command in the final line (l. 201).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Categorial traits:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Trees vs. Lianas, nonparametric and likelihood methods
load('E.mat'); %load wavelet variance for all species - columns are species, rows are spatial scales
load('dispbin.mat'); %load matrix identifying species as trees or lianas
load('Nsave_lianatree.mat'); %list of population sizes for each species, for weighting
dispnamesbin = ['tree '; 'liana']; %column 1 of dispbin is trees, column 2 is lianas
load('scale.mat'); %load list of spatial scales corresponding to E
cats(1).a=[0]; cats(2).a=[1 2]; %define likelihood tests - compare no groups ("0") to trees vs. lianas (columns 1&2 of dispmat)

niteruse=100; %number of bootstrapping iterations
weighttype='none'; %other options are 'log' or 'lin' for log or proportional weights, respectively
subs2=ones(1,size(E,2)); %optional masking vector for excluding species from analysis
%should be "1" for include species, "0" for don't include species
[dLL] = do_lik_boot_cat(E, dispbin, dispnamesbin, scale, Nsave_lianatree, niteruse, weighttype, subs2);
%make figures and extrace difference in log likelihood between H0 and HA

%%%%%%Trees vs. Lianas, fitting analytical model
load('Awfit.mat') %load covariance between scales and V for analytical model fitting

hold all;
doplot = 1; %plot outputs of fitted models
[totmat] = getcatfit(E, Nsave_lianatree, dispbin, scale, doplot, Awfit, weighttype, subs2); %fit analytical models I and II to data
%totmat output: rows show parameter estiamtes for cD, cK, and PI1.
%First two columns show mean paramter fits for model I for trees and lianas, respectively
%Columns 3-4 show standard deviation for parameters for model I
%Columns 5-6 and 7-8 show mean parameters and standard deviations for model II
%Columns 9 shows AIC for trees for model I and model II.
%Column 10 shows AIC for lianas.

figure(2); %make new figure
hist(dLL) %plot histogram of differences in log likelihood. LL>0 indicates significant difference among groups

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Continuous traits:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Tree DBH, nonparametric and likelihood methods
load('dbhquant'); %load standardized log maximum DBH for tree species
Nsave_tree = Nsave_lianatree(dispbin(:,1)==1); %extract N for all trees
Etree = E(:,dispbin(:,1)==1); %extract V for all trees

plotxname='standardized maximum dbh.'; %name for x axis on plot
unlog=1; %indicates that trait variable is log-scaled
h0=2; %kernel size - 0 uses "rule of thumb"

figure(3);

subs2=ones(1,size(Etree,2)); %use all dbh data
%subs2=((dbhquant>=quantile(dbhquant, 0.05)) & (dbhquant<=quantile(dbhquant, 0.95)))'; %use only mid-90% of DBH
[modfit]=do_lik_boot_cont(Etree, Nsave_tree, h0, dbhquant, niteruse, plotxname, unlog, weighttype, subs2); %make figures and extrace log likelihood
likeout = modfit.logLlstBOOT(1)-modfit.logLlstPERM; %difference in log likelihood between HA (boot) and H0 (perm)
figure(4);
hist(likeout);

%%%%%%Tree DBH, fitting analytical models
[modfit]=do_lik_boot_cont_fitting(Etree, Nsave_tree, h0, dbhquant, niteruse, Awfit, weighttype, subs2);
%modfit.totfit
%mean parameters, standard deviation, and AIC for models I and II for each
%trait value, formatted following outline in totmat above


%END