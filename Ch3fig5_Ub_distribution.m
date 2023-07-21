% Create example excess velocity plot for Ub figure in chapter 2 discussion
% select a file and extract the Ub data and make a bar chart
% Tobias Hasse October 30, 2020
% *************Figure 5 chapter 2*****************
% also see Ch3_lambda_fit.m, Ch3_Ub_fit.m, an_parm_tests.m, Main_Ch3.m,
% migration_model_TRH_Ch3.m

clear
% Ub_hist_max is hard coded because it was not saved in the output file.
% look this value up in Migration_model_TRH or look at the .png file for
% the Ub plot saved during the run to interpret the maximum bin of the
% histogram

Ub_hist_max = 1; % **** caution HARD CODED *******
% cd 'C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch3\other'
cd 'C:\Users\thasse\Documents\MATLAB\test\Ch3 diss'
% load 'Hasse_param_A3_Cfo0.024_testing_10_ka Do_other_more.mat'
% % load 'Hasse_param_A3_Cfo0.024_testing_10_ka Do_other.mat'
load 'Hasse_param_A3_Cfo0.024_testing_0.5_ka_Do_other_more.mat'

cur_dir = pwd;

% clean up the absub vector
if absub(end)>0
    if absub(end-1)>0
        disp(sprintf(strcat('The maximum of the histogram might have',...
            ' been too low and some of the absolute\n values of median',...
            ' velocity might be larger than the maximum histogram bin',...
            '\n Histogram max is %f end bins are'),Ub_hist_max))
        absub(end-10:end) % show the values
    else
        absub(end)=0; % delete the last value which is the histogram count  
                % (really it is the number of times Ub_hist_max was counted
                % and as long as the absub distribution reaches 0 near the
                % end of the vector there are no velocities recorded here
    end
end

% make a histogram and calculate the mean and median
figure(100)
    % must make the histogram to get the bin centers
    mh = histogram(abs([ub',Ub_hist_max]),100); 
    % save bin centers to plot histogram with bar chart
    bin_centers = mh.BinEdges + mh.BinWidth/2; bin_centers(end)=[]; 
    close 100 % close the figure
    pause(0.1)
    ub_pdf = absub/sum(absub); 
    ub_cdf = cumsum(ub_pdf);
    % the idx method was attempted to find the median, but can find only
    % the median bin (if searching for 0.5) and since the histogram is made
    % of continuous data, linear interpolation between bin centers seems
    % appropriate. The idx method is maintained because interp1.m must have
    % a strictly monotonic increasing vector
    idx= find(ub_cdf>.6,1,'first');    
    median_ub = interp1(ub_cdf(1:idx),bin_centers(1:idx),.5);
    mean_ub = sum(absub.*bin_centers)/sum(absub);
    
%% plot the histogram tops as circles, include the mean and median
figure(16)
    plot(bin_centers,ub_pdf,'color','k','linestyle','none',...
        'marker','.','markersize',8)
    yl=ylim;
    hold on
    plot( [median_ub median_ub], yl(2)./[20 1.1],'marker','d','color','k')
    plot( [mean_ub   mean_ub  ], yl(2)./[20 1.1],'marker','o','color','k')
    hold off
    title( sprintf( 'cfo = %0.4f, A = %d, Eo = %1.2e', Cfo,alphaao,Eo ) )
    ylabel('Probability per node')
    xlabel('Absolute value of Dimensionless Excess Velocity (U_b)')
    legend('Probability per node',...
        sprintf('Median %0.4f',median_ub),sprintf('Mean %0.3f',mean_ub),...
        'location','northeast')
    set(gcf,'color','w','position',[100 100 752 350]) 
        pause(0.1)

%% save figures
%     fig2svg('Ub 3 24 with mean.svg')
    print('-painters', '-dpng', '-r600', 'Ub 3 24 with mean.png')

