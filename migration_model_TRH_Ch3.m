function [river, nodecount, B, absub,bin_centers, freq_mig, out] = ...
    migration_model_TRH_Ch3(Cfo, alphaao, statstoolbox,outfilename,...
    sim_time_yrs,dt_yrs,save_dt,Ub_hist_max,X,Y,inplanformname,...
    linear_algorithm,scaleup)
disp('welcome to Schwenks model modefied by Hasse')
% This function was written by Jon Schwenk 2014-2015 
% and is available here 
% http://onlinelibrary.wiley.com/doi/10.1002/2014JF003252/full
% along with other functions as supplemetary info for the article:
% Schwenk, J., Lanzoni, S., & Foufoula?Georgiou, E. (2015). The life of a 
% meander bend: Connecting shape and dynamics via analysis of a numerical 
% model. Journal of Geophysical Research: Earth Surface.
% This is the main function for computing the migration of a meandering
% channel
% Updates by Tobias Hasse, tobiack@udel.edu, April 2015 include:
% 1) calling updated functions that are further optimized (for speed &/or 
%    memory) from the Schwenk functions, 
% 2) commenting out unused variables using "%TRH" to differentiate from 
%    Schwenk comments
% 3) adding some variables to help the new functions be more efficient (
%    e.g.) cell_row_offset
% 4) and iserting timers used to time different sections of the loop.  
%Updated functions have the suffix _TRH and are available from Tobias Hasse

% Additional updated by Tobias Hasse tobiack@udel.edu 2020-2021 to
% save additional parameters as part of Tobias Hasse's dissertation
% research for Chapter 2

% Tobi's last name changed during this dissertation.  It is possible that
% references to the initials TRA refer to files or functions labeled with
% the initials TRH, same with the name change from Ackerman to Hasse

% Note that the save_nodecount C section has been changed so that duplicate
% indecies are not saved.  This requires the updated atom_tracking_TRH
% function to be called from Model_Extract_and_Plot.
clock_time_start = clock;
%% Load and format coordinates, parameters
load params_meander.mat; % load reachwide variables, parameters
% S_valo = 0.000188; %0.0001506;  % parametric testing these slopes
% replicates D_ch etc with sending Do to flowfield rather than D_ch to 
% flowfield
% You can overwrite these parameters which are in the params*.mat file
% Cfo = .024;  % TRH update Beatton River CA data based on criteria k ~= 
%                epsilon o from Johanesson & Parker 1989 
% alphaao = 3;  % TRH ***** WARNING ***** alphaao has 1 added within 
%                 flowfield for unknown reasons (ref: Johannesson & Parker 
%                 1985). TRH: Schwenk used 10 not 3
if strcmp(outfilename(end-8:end),'recursive')
%     inplanformname = 'initialplanform'
% if a recursive initial planform is needed migration rates are very slow
    Eo = Eo*10;   
    disp(sprintf('increasing Eo %f',Eo))
    Eo
end
% Eo = Eo*2;   % TRH the new Cfo, A values, slow down the migration rate
% **** The updated Eo is in the params_meander.mat file ****
% Generate initial centerline coordinates
% X = [0:0.5:150]'*B;
% Y = randn(numel(X),1);
x_max_orig = X(end);
%% Boolean speed considerations: parallel and saving--what & how often?
% par_flow = 1;           % the flowfield can be parallelized
par_cut = 0;            % the intersections search can be parallelized 
                        % required by cutoff.m (no stats toolbox)
%TRH save_stats = 0;    % saves various statistics in structure called stat
save_riv_vars = 0;      % saves river variables Xcl, Ycl, mig_Xcl, mig_Ycl 
                        % (centerlines and migrations at each time step)
save_nodecount = 0;     % saves variables needed to perform node tracking 
                        % - used for atom identification
nodecount = 0;          %if not saving nodecount, gives something to return
% save_dt = 20;         % to save every nth time step
appx_conv_int = 1       % approximate convolution integral by truncating 
                        % its computation after some threshold distance 
                        % (previously hard-coded to be 100B)

save_mem = 0;        % uses single precision instead of double if activated

%TRH save_points = 1;   % number of points throughout simulation to save 
                        % data (adds time but prevents loss of data if the 
                        % run fails for some reason)
%TRH save_number = 0;   % used for accounting; do not change this variable 
                        % from 0.

%% Ch1 Simulation variables  % timing commented out bc read into function
% sim_time_yrs  = 100000;  %(TRH passed into function) length of 
                           % simulation, years
% dt_yrs = 0.1;            % time step, years, larger dt causes problems 
                           % in enforce_spacing.m
sim_time = ceil(sim_time_yrs/dt_yrs); % number of time steps, no dimensions 
dt = dt_yrs*365.25*24*3600; % time step in seconds
% disp_progress =  10000;  % display t every disp_progress time step
% disp_progress = floor( sim_time / 10 ); % show progress every 1/10th
disp_progress = floor( sim_time / 2 );    % only show 50% complete model

%% Prescribe boundary, initial conditions, constants, and parameters
% Constants
rho = 998.1;    % density of water, kg/m^3
g = 9.81;       % gravity, m/s^2

% Solve initial flow variables and parameters
% tortuosity calculated assuming straight valley
[~, tortuosity] = quicktor_TRH(X,Y);
Cf_ch = Cfo;                        % initial friction factor
S_val = S_valo;                     % initial valley slope
alphaa = alphaao;                   % alphaa doesn't change in time
W = B*2;                         % channel width, m, doesn't change in time
So = S_valo/tortuosity(1);          % initial stream slope
% initial reach-averaged depth 
% (TRH is this only truely Do because tortuosity ? 1?)
Do = (Qo/2/B)^(2/3)*(Cfo/g/So)^(1/3); 
Uo = Qo/(Do*W);                     % initial reach-averaged velocity
% half-width:depth ratio for straight channel characterized by valley slope
%TRH betaao = B/Do;                 % TRH (apparently unused)     
% initial (Froude #)^2 for straight channel characterized by valley slope
F2o = (Uo/sqrt(g*Do))^2;           

%% Spacing criteria - adjusting these could cause instabilities, other bugs

% ****** Given Chapter 2 of my (Tobias Hasse's) dissertation there is
% better theory to manage all of the spacing thresholds ********

% TRH search_excl_range speeds up remove_cutoffs_TRH.m
% search_excl_range is a number of nodes greater than search_radius
% donwstream, but less than the smallest cutoff from Schwenk 2015, fig 7
% TRH Here we have a bit of a rabbit hole:  
% Schwenk 2015, fig 7 is nondimensionalized using the length scale Lo
% Lo = Do/2Cfo, for Schwenk 2015 Cfo = 0.0036 and 
% Do = 2.282 using the equation above (1.0582 in his input params.mat file
% This means that Lo = 317, but Schwenk used 147 (personal comm Jan 2020)
% Using Lo = 317, and Schwenk figure 7, minimum Lcut is 150 B which is why
% I made the max search threshold much lower than this (at 50) to be 
% conservative 
% Using Lo = 147, minimum Lcut is 70B, however, based on the
% adjustments made by Tobias Hasse (table 1 of dissertation) to create
% smaller bends minimum Lcut might be 46. 
% Fortunately the search exlcusion range was 22.8 which is still less than
% 46, but not as conservative.
% Thresholds in params*.m file assume scaleup is 5.7, 
% since scaleup can vary then the spacing thresholds must be scaled also
dS_spacing_thresh = dS_spacing_thresh * scaleup / 5.7;
too_close_thresh  = too_close_thresh  * scaleup / 5.7;

% number of channel half-widths between centerline to detect cutoffs; 
% default value (2) is width of channel
search_radius = B * 2; 
search_excl_range = ceil(search_radius/too_close_thresh +2);
cell_row_offset = num2cell([search_excl_range:(numel(X)*4+...
    search_excl_range)]'); 
    % make the offsets vector big enough for some increse in nodes
if search_excl_range*dS_spacing_thresh/B > 30;  %50
    sprintf(strcat('remove_cutoffs.m is skipping %d channel half',... 
        'widths downstream, cutoffs can be as small as 46 half widths'),...
        ceil(search_excl_range*dS_spacing_thresh/B)) 
end
% TRH this approximates conv_int_thresh = 100 in flowfield_TRH
conv_node_thresh = ceil( 100 / dS_spacing_thresh * B * 7/4);  
%% Former spacing criteria - adjusting these could cause instabilities and
% % % % other bugs.  Included here, commented out, for completeness.  
% % % % Should give same results as the code above which replaces this 
% % % % (depends on inputs in params*.mat file)
% % % % TRH NOTE: the spacing thresholds are dimensionalized 15 lines below
% % % dS_spacing_thresh = scaleup*2/3;% when (CL node spacing)/B is larger 
% % %                                 % than this value, interpolate a node 
% % %                                 % there
% % % too_close_thresh = scaleup*0.1; % if two nodes are shorter than this 
% % %                                 % distance apart (in X and Y 
% % %                                 % directions independently) then 
% % %                                 % remove one of them
% % % % Cutoff threshold distance
% % % cutoff_thresh = max(2,dS_spacing_thresh/1.9); % number of channel 
% % %                                 % half-widths between centerline to 
% % %                                 % detect cutoffs; default value is 
% % %                                 % width of channel
% % % search_radius = cutoff_thresh*B;
% % % % TRH search_exl_range speeds up remove_cutoffs_TRH.m
% % % % TRH cell_row_offset is a number of nodes greater than 
% % % % search_radius donwstream, but less than the smallest cutoff from 
% % % % Schwenk 2015, fig 7
% % % search_excl_range = ceil(cutoff_thresh/too_close_thresh +2);
% % % cell_row_offset = num2cell([search_excl_range:(numel(X)*4+...
% % %     search_excl_range)]');      % make the offsets vector big enough
% % %                                 % for some increse in nodes
% % % if search_excl_range*dS_spacing_thresh > 50;
% % %    sprintf(strcat('remove_cutoffs.m is skipping %d channel half-',...
% % %         'widths downstream, cutoffs can be as small as 135'),...
% % %         ceil(search_excl_range*dS_spacing_thresh)) 
% % % end
% % % %TRH this approximates conv_int_thresh = 100 in flowfield_TRH
% % % conv_node_thresh = ceil( 100 / dS_spacing_thresh * 7/4);  
% % % dS_spacing_thresh = dS_spacing_thresh*B;  
% % % too_close_thresh = too_close_thresh*B;   
%% Curvature boundary conditions, smoothing
bc1o = 0;           % upstream initial and boundary condition for curvature
bcset = 1;          % 1 for periodic, 2 for C(1)=0, 3 for small random 
                    % Gaussian, 4 for linear interpolation
num_C_smooth = 1;   % number of times to smooth curvature each iteration

%% Initialize model variables
D_ch  = nan(1,sim_time);    D_ch(1) = Do;
S_ch  = nan(1,sim_time);    S_ch(1) = S_valo;
U_ch  = nan(1,sim_time);    U_ch(1) = Uo;
F2_ch = nan(1,sim_time);    F2_ch(1) = F2o; 
%TRH optoinal for measuring reach wide parameters at each step
chan_len = nan(1,sim_time); bends = nan(1,sim_time); 
% % % %%TRH DEBUG variables
% % % mig_mag_max = 3;
% % % freq_mig = zeros( 1,31 );
bc1   = nan(1,sim_time);   bc1(1) = bc1o;
%% Implement single precision - quicker runs, smaller output files
% TRH this might not matter on 64 bit machines.  I think 64 bit machines 
% use the same amount of RAM for single and double precision
if save_mem == 1
    X = single(X);
    Y = single(Y);
    U_ch = single(U_ch);
    D_ch = single(D_ch);
    Cf_ch = single(Cf_ch);
    S_val = single(S_val);
    Qo = single(Qo);
    B = single(B);
    alphaa = single(alphaa);
    bc1(1) = single(bc1(1));
end
    
%% Initialize saved variables' structure arrays
if save_nodecount
    save_dt = 1;
   disp(sprintf('Changed to save every step due to nodecount requirement'))
    nodecount=[]; % make sure nodecount is blank before reassigning
    nodecount( ceil( sim_time / save_dt ) ).A_cutoff_rem = [];
    nodecount( ceil( sim_time / save_dt ) ).B_spacing_ins = [];
    nodecount( ceil( sim_time / save_dt ) ).C_duplicate_rem = [];
%TRH     cut_idx = cell(1,sim_time);
end

if save_riv_vars
    river( ceil( sim_time / save_dt ) ).Xmig = [];
    river( ceil( sim_time / save_dt ) ).Ymig = [];
end
river( ceil( sim_time / save_dt ) ).Xcl = [];
river( ceil( sim_time / save_dt ) ).Ycl = [];
% TRH optional for measuring wavelength of each bend only saved time steps
% this method was developed after chapter 3 was completed and is included
% as an additional (and perhaps superior) method
waves(numel(river)).length=[]; 

%% Begin modelling
% Initialize centerline  
%TRH if not saving riv_vars (L95) then these are not preallocated
save_t = 1;
river(save_t).Xcl(:,1) = X;
river(save_t).Ycl(:,1) = Y;
%% migration rate and velocity variables
mig_mag_max = 3;
freq_mig = zeros( 1,31 );
absub = zeros(1,100); 
%% DEBUG variables
% clock_curv = 0; clock_flow = 0; clock_angle = 0; clock_cutoff = 0; 
% clock_smooth = 0; clock_space = 0; clock_update = 0; tic
%% additional setup options
meas_ub_migration = 0;  % dont measure migration at beginning of simulation
% load 120kplanform     % optional load developed planform
algorithm_start = clock 

for t = 1:sim_time   % ****************************************************
%     line_266_t_is = t 
% keyboard % continue with dbcont, quit with dbquit
if isequal(t, ceil(sim_time*.4))    % if 40% through the simulation
    if ~meas_ub_migration           % and no prior cutoffs
        meas_ub_migration =1;       % force measuring ub and migration rate
        strt = t;                   % record when measurements start
        disp(sprintf('forcing excess velocity measurements at t = %d', t))
        % if next is still false, then no cutoffs yet & no prior recursions
%         keyboard
        if ~strcmp(outfilename(end-8:end),'recursive')
            % for this method to work well, the sequence of simulations
            % must go systematically from simulations with small bends to 
            % simulations with large bends
            disp('loading prior planform, reset t =1')
            load outplanform.mat % this file is overwritten by each sim
            % save the recursive starting planform
            save([outfilename,' in planform.mat'],'X','Y',...
                'outplanformname','-v7.3') 
            % rename variables for files to document the recursive sim.
            inplanformname = outplanformname 
            outfilename = [outfilename,' recursive'];
            figure(102)             % show the recursive input planform
               plot(X,Y)
              title(sprintf('loading initial planform: %s',inplanformname))
               xlabel(sprintf('Current Cfo is %f, A is %d',Cfo, alphaao))
               axis equal
               drawnow
            % recursive simulations have very low ub, so reduce the max
            Ub_hist_max = Ub_hist_max /2; 
            % call migration_model_TRH_Ch3 recursively 
            [river, nodecount, B, absub,bin_centers, freq_mig, out] = ...
                migration_model_TRH_Ch3(Cfo, alphaao, statstoolbox,...
                outfilename,sim_time_yrs,dt_yrs,save_dt,Ub_hist_max,X,Y,...
                inplanformname,linear_algorithm,scaleup)
            return              % after the recursive function call, return
        end
    end
end
tic
%% Solve hydrodynamics
% (dimensionalized) curvature and related variables
[C, S_cum, dS] = curvars_TRH(X,Y,bc1(t),bcset, num_C_smooth); 
% if bcset == 1, the upstream curvature at the next time step is equal to 
% the downstream-most calculable (non-interpolated) curvature value, 
% a la a periodic boundary condition
bc1(t+1) = C(end-1); 
% clock_curv = clock_curv + toc; tic
chan_len(t) = S_cum(end); 
% Schwenk 2014 flowfield solution
% [ub,Lo_in_Bs,tf1,tf2] = flowfield_Schwenk(B,C,D_ch(t),S_cum,alphaa,...
%     F2_ch(t),Cf_ch,appx_conv_int);

% Hasse 2015 flowfield solution (based on Schwenk's function)
if linear_algorithm
    % TRH, correct method using hypothetical straightened reach parameters
    % this closeley approximates Parker & Andrews 1985 Beatton River
    % ******** THIS [ub] FUNCTION CALL RETURNS CORRECT RESULTS *********
    [ub] = flowfield_TRH(B,C,Do     ,S_cum,dS,alphaa,F2o     ,Cf_ch,...
        appx_conv_int,conv_node_thresh); 
else
    % this is from Schwenk
    % TRH, Schwenk's method using sinuous channel parameters for simulation 
    % This method was used in Chapter 1 of Tobias Hasse's dissertation
    % This algorithm returns extra erroneously large bends due to the use
    % of Dch and F2_ch from the sinuous reach
    [ub] = flowfield_TRH(B,C,D_ch(t),S_cum,dS,alphaa,F2_ch(t),Cf_ch,...
        appx_conv_int,conv_node_thresh); 
end

% TRH in the past I measured meander wavelength using fft on the curvature
% series, concatenating many streamlines together to make a long series the
% method below measures bends without fft, but rather uses ub
% ub is a smothed version of curvature, some small bends may be included
% bends is tracking the number of meanders at each t, this will be
% converted to meander wavelength later in this function
temp=zeros(size(ub)); 
temp(ub>0)=1; 
bends(t) = sum(abs(diff(temp)))/2; %meander bends,(pair: 1 left, 1 right)

% clock_flow = clock_flow + toc; tic %DEBUG Lo_mn = min(Lo_mn,Lo_in_Bs); 
% % uncomment these lines if not running both flowfields
% Ro = min(abs(1./C)); % minimum radius of curvature along reach
% nu0 = B/Ro;
% ub = flowfield_Zolezzi_u(X, Y, C, S_cum, dS, B, D_ch(t), Cf_ch(t),...
%     S_ch(t), nu0);
%% Perform migration given flow field

Eps_cl  = Eo*ub; % migration rate along the centerline

% Compute erosion along outer banks (n=1)
% angles computed under assumption of constant valley direction found via 
% linear regression
angle = cl_angles(X,Y); 
% TRH multiplying by B below is in error and is what was done by Schwenk
% and by Hasse in dissertation chapter 2 & 3. 
% These lines remain for reproducibility of results

% ********** THE RESULT OF THIS ERROR IS THAT ACTUAL Eo IS B*Eo **********

dX_cl = -Eps_cl.*sin(angle)*U_ch(t)*dt*B; % dimensionalized dx
dY_cl = Eps_cl.*cos(angle)*U_ch(t)*dt*B; % dimensionalized dy

% ******** USE THE TWO LINES BELOW, NOT ABOVE, FOR CORRECT RESULTS *******

% dX_cl = -Eps_cl.*sin(angle)*U_ch(t)*dt; % dimensionalized dx
% dY_cl = Eps_cl.*cos(angle)*U_ch(t)*dt; % dimensionalized dy

if meas_ub_migration %TRH measure migration rate of nodes
    % this appears to be much greater than the migration rate orthogonal 
    % to streamflow (TRH June 2016)
    mig_mag = sqrt( dX_cl.^2 + dY_cl.^2 );     
    fh = figure('visible','off'); % figure(100)
        mh = histogram([mig_mag',mig_mag_max],length(freq_mig));
        freq_mig = freq_mig + mh.Values;
        close(fh)
    fh = figure('visible','off'); % figure(100)
        mh = histogram(abs([ub',Ub_hist_max]),100);
        % notice the last element in absub should be deleted later since it
        % is accumulating the frequency of times the histogram is created.
        % I did not do this.  
        % This *****ERROR***** is about 1/number of nodes in channel and is
        % minimal since there are hundreds to thousands of nodes pr channel
        absub=absub + mh.Values;
        close(fh)

end
Xf = X; Yf=Y; %TRH for saving final X,Y coordinates that match ub
mig_Xcl = X + dX_cl; % migrated X coordinates
mig_Ycl = Y + dY_cl; % migrated Y coordinates
% clock_angle = clock_angle + toc; tic
%% Locate and perform cutoffs
Npre = numel(mig_Xcl); 
if statstoolbox
%     [X, Y, cut_idx] = remove_cutoffs(mig_Xcl, mig_Ycl, cutoff_thresh, B);
    [X, Y, cut_idx,cell_row_offset] = remove_cutoffs_TRH(mig_Xcl, ...
        mig_Ycl, search_radius, search_excl_range, cell_row_offset);
else
%     disp(strcat('WARNING: The function cutoff smooths cutoffs inside',...
%         ' but using different parameters than smooth_after_cutoffs below'))
    % 2023 cutoff is cutting off when it should not. removing large
    % sections of river centerline
%     [X,Y,n_cutoffs,Xl,Yl,Xr,Yr,cut_idx,ncut] = cutoff(mig_Xcl,mig_Ycl,...
%         B,par_cut);
    % make TRH_rm_cutoffs call here for no stats pack, pbbly faster than 
    % cutoff b.c. similar speed to remove_cutoffs
    [X, Y, cut_idx] = TRH_rm_cutoffs(mig_Xcl, mig_Ycl, search_radius,...
        search_excl_range);

end

% smooth around cutoffs
% clock_cutoff = clock_cutoff + toc; tic
if isempty(cut_idx) == 0
    [X, Y, ~] = smooth_after_cutoffs(X, Y, cut_idx, Npre);
    % TRH Chapter 3. If previously not measuring ub and migration rates, 
    % now that there has been a cutoff, start measuring those variables
    if ~meas_ub_migration 
        meas_ub_migration = 1; 
        strt = t;
    end
    if save_nodecount
        nodecount(t).A_cutoff_rem = cut_idx;
    end
end    
% clock_smooth = clock_smooth + toc; tic
%% Enforce maximum streamwise spacing threshold 
% TRH combine enforce spacing [X,Y,interp_idcs,cut_idcs] = 
%                     enforce_spc(X,Y,dS_max_spc,dS_min_spc,save_nodecount)
[X, Y, interp_idcs] = enforce_spacing(X, Y, dS_spacing_thresh,...
    save_nodecount);
if save_nodecount
    nodecount(t).B_spacing_ins = [interp_idcs]; 
end
%% Enforce minimum streamise spacing threshold
dS_ibn = sqrt(diff(X).^2+diff(Y).^2); 
% find indices where centerline nodes are too close
remove_duplicates = find(dS_ibn < too_close_thresh); 
if save_nodecount
    if isempty(remove_duplicates) == 0
        %TRH repmat is not needed for the updated atom_tracking, which
        %doesn't need two copies of the removed node at lines 132 & 135
        nodecount(t).C_duplicate_rem =  remove_duplicates; 
        %TRH sort(repmat(remove_duplicates,1,2)); 
        % the sort(repmat()) is because each adjustment of nodes has a 
        % starting and ending node, so if it's only one node being removed 
        % there should be two identical entries in the structure
    end
end
X(remove_duplicates)=[];
Y(remove_duplicates)=[];
% clock_space = clock_space + toc;tic
%% Store variables
if ~mod(t,save_dt) %TRH saves river planform, but not every time step
    save_t = save_t + 1;  % increment river save counter
    river(save_t).Xcl = X;
    river(save_t).Ycl = Y;
    % TRH length of each bend to left or right ÷ half width and tortuosity 
    % results in a struct called 'waves' which contains the straight down
    % the valley wavelength of each bend in the channel 
    waves(t).length=(diff(S_cum(abs(diff(temp))==1))/(B*tortuosity));
    if save_riv_vars == 1
        river(save_t).Xmig = dX_cl/dt; % migration distance
        river(save_t).Ymig = dY_cl/dt; % migration distance
    end
end
%% Compute new tortuosity and channel slope
% tortuosity calculated assuming straight valley
% keyboard
[~, tortuosity] = quicktor_TRH(X,Y);    
S_ch(t+1) = S_val/tortuosity;
%% Update flow variables
D_ch(t+1) = (Qo/W*sqrt(Cf_ch/g/S_ch(t+1)))^(2/3);
U_ch(t+1) = Qo/(D_ch(t+1)*W); 
% Froude no^2 for straight channel characterized by valley slope
F2_ch(t+1) = (U_ch(t+1)/sqrt(g*D_ch(t+1)))^2;  
% clock_update = clock_update+toc; tic
%% Display progress. Most of the code below here is written by TRH
if rem(t,disp_progress) == 0
    % TRH The river seems to extend off the downstream end of the valley 
    % over long simulation times, this truncates to original valley length
    % occasionally only at disp_progress interruptions
    new_end = find( X > x_max_orig , 1 ,'first');
    [new_end numel(X)]
    if new_end                       % not empty
        X(new_end:end)=[];           % remove the tail of the simulation
        Y(new_end:end)=[];
    end

    % Display progress:
    pct = t/sim_time*100;
    elapsed_algorithm_time = datenum(clock) - datenum(algorithm_start);
    disp([num2str(pct), ' % finished.',...
        sprintf(' Current Time:  %s \nest completion %s',...
        datestr(datenum(clock)) , datestr( datenum( algorithm_start ) + ...
        elapsed_algorithm_time * 100 / pct ) ) ] )
    memry = memory
    figure(1)
        plot(X,Y)
        title(sprintf('Final planform t is %d cfo is %f A is %d',t,Cfo,...
            alphaao))
        xlabel('Distance (m)'); ylabel('Distance (m)');
        axis equal
        drawnow
    figure(2)  
        bar(absub(1:end-1))
        title('Excess velocity (scaled)')
        xlabel(sprintf('Excess velocity (scaled %d ÷ %0.1f) (m/s)',...
            length(absub)-1,Ub_hist_max))
        ylabel('Frequency')
        drawnow
    figure(3)
        plot(chan_len)
        title('channel length')
        xlabel('Model time steps'); ylabel('Length (m)');
        drawnow
%% If all the timer variables and tic toc are uncommented, (above) then you 
% can display them with the code below including labels
%     a = [clock_curv clock_flow clock_angle clock_cutoff clock_smooth...
%         clock_space clock_update];
%      timers = strcat('clock_curv clock_flow clock_angle clock_cutoff',...
%          ' clock_smooth clock_space clock_update total numel(riv)')
%     [a sum(a) numel(river)]

%     save('par_x_y.mat','X','Y','t','-v7.3'); %optional save output
end % if dislay progress
end % for t
%%
%% clean up some stuff for parametric testing runs
% if run was recursive, inplanformname is already assigned, otherwise it
% should be assigned here
if ~strcmp(outfilename(end-8:end),'recursive')
    inplanformname = 'initialplanform'
end
% strt is the time t @ first cutoff (or 40 % of simulation whichever is
% less [see line 235]). We want a long data set, but want to avoid boundary
% condition effects of model spinup.  If the first cutoff happened early we
% can afford to push back the strt time for calculating simulation
% parameters
if strt*5<t
    strt = strt*2;  % strt = t 
end

%     figure('visible','off') % figure(100)
figure(100)
    % must remake the histogram to get the bin centers
    mh = histogram(abs([ub',Ub_hist_max]),100); 
    % save bin centers to plot histogram with bar chart
    bin_centers = mh.BinEdges + mh.BinWidth/2; bin_centers(end)=[]; 
    csv = cumsum(absub)/sum(absub); % cumulative distribution
    % the idx method was attempted to find the median, but can find only
    % the median bin (if searching for 0.5) and since the histogram is made
    % of continuous data, linear interpolation between bin centers seems
    % appropriate. The idx method is maintained because interp1.m must have
    % a strictly monotonic increasing vector
    idx= find(csv>.6,1,'first');    
    median_ub = interp1(csv(1:idx),bin_centers(1:idx),.5);


% reach averaged meander wavelength measured along the sinuous channel for
% each channel planform at each time step in the simulation
lambda_c = chan_len./bends/70;          % wavelength along sinuous channel
lam_c = mean(lambda_c(strt:end));       % average sinuous wavelength
tortuosity= S_val./S_ch;                % tortuosity at every timestep
tor = mean(tortuosity(strt:end));       % average over latter part of sim
lam = lam_c/tor;                        % average wavelength down valley
% standard deviation of the reach averaged wavelength (the time series for
% sinuous wavelength is used, but converted to down valley wavelength using
% the variable 'tor'.  This is included here for informational purposes,
% but not really relevant, because what we would REALLY want is the
% standard deviation of the size of individual wavelengths (bend pairs) and
% this is only given via the variable 'waves' (see below)
lam_std = std(lambda_c(strt:end))/tor;  % not used or returned from fuction
U = mean(U_ch(strt:end));               % average downstream velocity
H = mean(D_ch(strt:end));               % average channel depth
% [Cfo alphaa lam_c 0 U H tor lam];     % DEBUG print the variables
% put the variables in a struct to return from the function
out.lam_curv=lam_c;
out.lam_crtsn = lam;
out.tor=tor;
out.Ub_median = median_ub;
out.U_mean = U;
out.H_mean = H;
out.Cfo = Cfo;
out.A_input = alphaao;

% the outplanform is overwritten and used only as an input planform for
% recursive runs.  This was done in an attempt to force simulations to
% meander when they would otherwise straighten out the minimal random noise
% in the initial planform.  The recursive method requires incrementally
% stepping through the A (alphaao) and Cfo variables from simulations with
% small bends to simulations with large bends.  This did not work very well
% and was not used in Chapter 2 of Tobias Hasse's dissertation
outplanformname=outfilename;
save('outplanform.mat','X','Y','outplanformname','-v7.3')
%% finalize figures
figure(1)
    plot(X,Y)
    title(sprintf('Final planform %s',outfilename))
    xlabel('Down valley distance (m)')
    ylabel('Cross valley distance (m)')
    axis equal
    print('-painters', '-dpng', '-r600', sprintf('%s planform.png',...
        outfilename))

figure(2)
    bar(bin_centers(1:end-1),absub(1:end-1))
    yl=ylim;
    text(median_ub,yl(2)/1.2,'+')   % plot the median value as a + symbol
    title('Frequency distribution of excess velocity |ub|') 
    ylabel( sprintf( 'Frequency, cfo = %0.4f, A = %d, Eo = %1.2e',...
        Cfo,alphaao,Eo ) )
    xlabel(sprintf( 'Excess velocity (m/s), median (+) is %f m/s',...
        median_ub))
%     print('-painters', '-dpng', '-r600',... 
%         sprintf('%ska Ub',num2str(floor(sim_time_yrs/1000))))
    print('-painters', '-dpng', '-r600', sprintf('%s Ub.png',outfilename))

figure(4)
    bar(freq_mig)
    title('Nodal migration rate')
    ylabel('Frequency')
    xlabel(sprintf('Meters per %1.2f years',dt_yrs ) )
    drawnow
%     print('-painters', '-dpng', '-r600',...
%        sprintf('%ska nodal migration',num2str(floor(sim_time_yrs/1000))))
    print('-painters', '-dpng', '-r600', ...
        sprintf('%s nodal migration.png',outfilename))
    %%
figure(5) % reach averaged wavelength at all time steps
    wavelength_straight = chan_len .* S_ch( 1:length( bends )) ./ ...
        ( bends * S_val * 2 * B );
    plot([1:length(bends)],bends,[1:length(bends)],wavelength_straight)
    title(sprintf(strcat('Number of bends, reach averaged',...
        ' wavelength: %1.1f ± %0.2f of time average'),...
            mean(wavelength_straight,'omitnan'),...
            std(wavelength_straight,'omitnan')))
%         nanmean(wavelength_straight(strt:end)),...
%         nanstd(wavelength_straight(strt:end))))
    xlabel('Model time steps');
    ylabel(sprintf('Down valley (2B), Bends'));
    lh = legend('Number of bends in reach',...
       'Reach averaged wavelength (in channel widths)','location','north');
    set(lh,'color','none');
    drawnow
%     print('-painters', '-dpng', '-r600', ...
%         sprintf('%s reach avearged wavelength.png',outfilename))   

figure(6)   % histogram: down valley wavelength of each bend from the time 
            % steps selected by save_dt this could be undersampling of very
            % short simulations
%     all_waves(1:10000)=[];  %suggest trimming off early measurements
    all_waves=vertcat(waves((ceil(strt/save_dt)):end).length);
    mh = histogram(all_waves);
    % compute bin centers to find median
    bin_centers = mh.BinEdges + mh.BinWidth/2; bin_centers(end)=[]; 
    wave_cdf = cumsum(mh.Values)/sum(mh.Values); % cumulative distribution
    % the idx method was attempted to find the median, but can find only
    % the median bin (if searching for 0.5) and since the histogram is made
    % of continuous data, linear interpolation between bin centers seems
    % appropriate. The idx method is maintained because interp1.m must have
    % a strictly monotonic increasing vector
    idx= find(wave_cdf>.6,1,'first');    
    median_allwave_lam = interp1(wave_cdf(1:idx),bin_centers(1:idx),.5);

    yl=ylim;
    text(median_allwave_lam,yl(2)/1.2,'+')   % plot median value (+ symbol)
    title(sprintf(strcat('Down valley wavelength of each bend avg:',...
        ' %0.2f ± %0.3f of bend average \n Compare to: %0.2f',...
        ' the average reach averaged wavelength'),...
            mean(all_waves,'omitnan'),std(all_waves,'omitnan')))
%         nanmean(all_waves),nanstd(all_waves),lam))
    ylabel( sprintf( 'Frequency, cfo = %0.4f, A = %d, Eo = %1.2e',...
        Cfo,alphaao,Eo ) )
    xlabel(sprintf(strcat('Down valley wavelength (2B), median (+) is',...
        ' %0.2f m/s'), median_allwave_lam))
    drawnow
%     print('-painters', '-dpng', '-r600', ...
%         sprintf('%s meander bend wavelegth histogram.png',outfilename))
%%
clock_time_end = clock;
end_time_is = num2str(clock_time_end)
duration = clock_time_end - clock_time_start

%% save run parameters
% taking penultimate coordinates which match final ub calculation
X=Xf; Y=Yf; 
save(sprintf('%s_other.mat',outfilename),'scaleup','save_dt','dt_yrs',...
    'sim_time_yrs','Cfo','alphaao','Eo','-v7.3')
save(sprintf('%s_other_more.mat',outfilename),'scaleup','save_dt',...
    'dt_yrs','sim_time_yrs','Cfo','alphaao','Eo','Cf_ch','S_val',...
    'alphaa','So','Do','Uo','F2o','D_ch','S_ch','U_ch','F2_ch',...
    'bends','all_waves','chan_len','absub','X','Y','ub','out',...
    'inplanformname','-v7.3')
end_model = memory
beep on % tell me you're finished
beep;pause(1);beep;
% to DEBUG and play with the variables at the end of the function
% keyboard % to quit: dbquit, to continue: dbcontinue
end % function 

