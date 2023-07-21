%% Main CHAPTER ************ 3 ********************
% Available at: https://github.com/TobiHasse/MRLATS.
% Purpose:  This script will call other functions to create the suite of 
%           meandering river simulations and step through all the
%           simulation and analysis of Chapter 3 of Tobias Hasse's 
%           Dissertation.
% Attribution: The meandering river simulation code has been edited from
%           the original form which was published by Jon Shwenk:
%           http://onlinelibrary.wiley.com/doi/10.1002/2014JF003252/full
%           along with other functions as supplemetary info for the 
%           article: Schwenk, J., Lanzoni, S., & Foufoula?Georgiou, E. 
%           (2015). The life of a meander bend: Connecting shape and 
%           dynamics via analysis of a numerical model. 
%           Journal of Geophysical Research: Earth Surface.
%
% File repository: Some of the steps in the simulation take a lot of
%           computational time, particularly due to high RAM demands.  To
%           shortcut this, the user is encouraged to download select output
%           files from the file repository:
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     2015 - October 2021

%% generate parameter input files
% rather than putting parameters directly in *.mat files, I have created
% commented code which will generate input parameter files.  
% This allows you to see commentary and variable definitions, or
% make adjustements for your simulations

close all
clear
dir_run_output = 'C:\Users\thasse\Documents\MATLAB\test\Ch3 diss';
% dir_run_output = 'C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch3\other';
cd(dir_run_output)

save_params_meander()               % parameters for meander model
% save_params_meander_Schwenk()     % as used by Schwenk (with commentary)

save_initial_planform()             % initial random planform 
                                    % included for (reproducibility)

%% initialize variables for meander model

do_you_have_stats_toolbox = 0;  % enter 1 if you have the Statistics 
                                % Toolbox, else enter 0
sim_time_ky   = .5; % 10        % Adjust this to change simulation length
sim_time      = sim_time_ky * 1000;
dt_save       = 12;             % save every ## iterations
dt_save_years = 2;
dt            = dt_save_years / dt_save;
                                    
%% Preallocate variables for output
my_data.A_input   = [];
my_data.Cfo       = [];
my_data.H_mean    = [];
my_data.U_mean    = [];
my_data.Ub_median = [];
my_data.tor       = [];
my_data.lam_curv  = [];
my_data.lam_crtsn = [];

% Cfo and A parameter space for Tobias Hasse's dissertation Chapter 3.
% They must be vectors of the same length, and data pairs represent 
% settings for a meandering river simulation
CFO =[ .0036 .024  .03   .014 .01  .007  .005 .0036 .0025 .00125 .04  ...
       .02   .014  .01   .007 .005 .0018 .03  .02   .014  .01    .007 ...
       .005  .0036 .0023 .044 .032 .019  .014 .01   .007  .005   .003 ...
       .046  .03   .02   .014 .01  .005  .046 .03   .02   .01 ];
Ain =[ 10     3    16    16   16   16    16   16    16    16     10   ...
       10    10    10    10   10   10     6    6     6     6      6   ...
        6     6     6     3    3    3     3    3     3     3      3   ...
        1     1     1     1    1    1     0    0     0     0 ];
%% Run migration model
CFO = [ .04 .024 .01 ]; % temporary for code checking small parameter space
Ain = [ 10  3    1   ]; %

linear_algorithm = true;
for ii = 1:2 % both algorithms
    start_new_algorithm = memory
    alg = 'Do';               % first time, use linear parameters algorithm
    desc = sprintf(strcat('Results from linear parametric testing.',...
        ' low value runs completed using \nrecursive planforms and',...
        ' higher Eo. Double check: algorithm = %s'), alg);
    scaleup = 5.7;      % reset in case it was changed to 1 on prev loop
    if isequal(ii,2)
        linear_algorithm = 0; % second time, switch to sinuous algorithm
        alg = 'Dch';
    desc = sprintf(strcat('Results from sinuous parametric testing.',...
        ' low value runs completed using \nrecursive planforms and',...
        ' higher Eo. Double check: algorithm = %s'), alg);
    end
for i=1:length(CFO) % parameter space
    start_new_model = memory
    outfile = sprintf('Hasse_param_A%s_Cfo%s_testing_%s_ka_%s',...
        num2str(Ain(i)),num2str(CFO(i)),num2str(sim_time_ky),alg)
    if Ain(i)>20
            Ub_hist_max = 5
    elseif Ain(i)>9
            Ub_hist_max = 2
    elseif Ain(i)>2.9
            Ub_hist_max = 1
    else
            Ub_hist_max = 0.001  %0.01
    end
    scaleup = 5.7;
    % WARNING without changing spacing thresholds, scalup can remove the
    % entire stream centerline
    % for the tightest bends, especially using the linear algorithm, I
    % suggest using a smaller scaleup factor. The code below was not used
%     if linear_algorithm && CFO(i)*Ain(i) > 0.138 % this is for the linear 
%         scaleup = 1;           % which makes such tight bends that the 
%     end                        % node spacing needs to be closer together
% initialplanform available at: https://doi.org/10.5281/ZENODO.5651841.
load initialplanform 
inplanformname = 'initialplanform';
X = scaleup*Xo;      
Y = Yo;              
    
[river, nodecount, B, absub,bin_centers, freq_mig, out] = ...
    migration_model_TRH_Ch3(CFO(i), Ain(i), do_you_have_stats_toolbox,...
    outfile,sim_time,dt,dt_save,Ub_hist_max,X,Y,inplanformname,...
    linear_algorithm,scaleup); %
% When restarting appending data to previous output is useful:
% my_data.A_input   = [ my_data.A_input   , out.A_input   ];
% my_data.Cfo       = [ my_data.Cfo       , out.Cfo       ];
% my_data.H_mean    = [ my_data.H_mean    , out.H_mean    ];
% my_data.U_mean    = [ my_data.U_mean    , out.U_mean    ];
% my_data.Ub_median = [ my_data.Ub_median , out.Ub_median ];
% my_data.tor       = [ my_data.tor       , out.tor       ];
% my_data.lam_curv  = [ my_data.lam_curv  , out.lam_curv  ];
% my_data.lam_crtsn = [ my_data.lam_crtsn , out.lam_crtsn ];
my_data.A_input(i)   = out.A_input   ;
my_data.Cfo(i)       = out.Cfo       ;
my_data.H_mean(i)    = out.H_mean    ;
my_data.U_mean(i)    = out.U_mean    ;
my_data.Ub_median(i) = out.Ub_median ;
my_data.tor(i)       = out.tor       ;
my_data.lam_curv(i)  = out.lam_curv  ;
my_data.lam_crtsn(i) = out.lam_crtsn ;
cell2mat(struct2cell(my_data))
i
%
figure(99)
    set(gcf,'position',[10 40 560 420])
    scatter(my_data.Cfo,my_data.A_input+2,my_data.Ub_median*1000,...
        log(my_data.Ub_median),'filled')
    ch = colorbar;
    ylabel(ch,'log( median excess velocity )')
    ylabel('True A')
    xlabel('Cfo')
    title(sprintf('Median Ub. %d out of %d runs complete',i,length(CFO)))
    
%     lab = num2str([my_data.Ub_median]');
%     lab = cellstr(lab(:,2:4));
%     strfind(lab,'.')
    la = [];
    % this relooping is inefficient but allows for restarts
    for j = 1:numel(my_data.Cfo) 
       lstr = num2str(my_data.Ub_median(j));
       la = [la, cellstr(lstr(2:5))];
    end
% Hack to fix the only label greater than 1
%     lstr = num2str(my_data.Ub_median(3));
%     la(3) = cellstr(lstr(1:4));

    text(my_data.Cfo,my_data.A_input+2,la)
    set(gca,'xscale','log','yscale','log')
    xlim([0.001 0.05])
    ylim([1 20])
    drawnow
%
end % i parameter space

%% save my data
% dir_run_output = 'C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch3\other';
% % change to new directory for files to use next time
% cd 'C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch3\data';

% Note the "_" at the start of the "Parametric testing" file avoids
% overwriting important output.  Used for testing only.  If you leave this
% here by mistake, find the file and edit the name by deleting "_"
outfile = sprintf('_Parametric testing %s.mat',alg)
desc 
save(outfile,'my_data','desc','-v7.3')

% change back to directory for all run output files
cd(dir_run_output)

end % ii algorithms

% outfile = 'double check 2016 params'
% [river, B, mxub, freq_mig] = check2016run( do_you_have_stats_toolbox,...
%     outfile,sim_time,dt,dt_save); %

% i=1; CFO = 0.0036; Ain = 10;
% outfile = 'Schwenk low slope'

% i=1; CFO = 0.024; Ain = 3;
% outfile = 'Hasse 2016 low slope'
% [river, nodecount, B, mxub, freq_mig] = migration_model_TRH(CFO(i),...
%     Ain(i), do_you_have_stats_toolbox,outfile,sim_time,dt,dt_save); %
%% Call scripts to make figures for Chapter 2
% these scripts will clear variables and read data in from file
% To make publication figures I saved each figure to .svg format and made
% additional formatting using the open source software: Inkscape

% Plot the reach averaged down valley wavelength of all the simulations
Ch3fig1_lambda_bubbles

% Show the 3 parameter fit models 
coefficients = 3; %can choose 0, 1, 2, or 3
[fobj,rs]=Ch3_lambda_fit('linear' ,coefficients) 
[fobj,rs]=Ch3_lambda_fit('sinuous',coefficients)

coefficients = 0;   % to show the quality of fit for the model with no 
                    % adjustable fit coefficients
[fobj,rs]=Ch3_lambda_fit('linear' ,coefficients)
[fobj,rs]=Ch3_lambda_fit('sinuous',coefficients)

% Plot the predicted vs observed wavelengths
Ch3fig2_lambda_pred_obs

% Show the 3 parameter fit models 
coefficients = 3; %can choose 0, 1, 2, or 3
[fobj,rs]=Ch3_Ub_fit('linear' ,coefficients) 
[fobj,rs]=Ch3_Ub_fit('sinuous',coefficients)

coefficients = 0;   % to show the quality of fit for the model with no 
                    % adjustable fit coefficients
[fobj,rs]=Ch3_Ub_fit('linear' ,coefficients)
[fobj,rs]=Ch3_Ub_fit('sinuous',coefficients)

% Plot the predicted vs observed excess velocity
Ch3fig3_U_pred_obs

%% Plot the contours of the storage timescale Ts
% You may want to run this multiple times and check the command window for
% contour label values
Ch3fig4_Ts_contours_U

%% Plot the absolute value of the dimensionless Ub distribution for one
% simulation.  Currently depends on simulation completed with Cfo = 0.24
% and Ain = 3, but could be adjusted to read in other files
Ch3fig5_Ub_distribution

%% Open each output file and grab a parameter to compare between 
% simulations, much faster than opening it manually.  You will have to edit
% the script to select the variable of interest
Ch3Check_output_files

