% create storage timescale Ts field plot with countours and simulation dots
% from the equation for Ts in Chapter 2, solving for A
% Tobias Hasse Dec 4, 2020 
% *************Figure 4 chapter 2*****************
% also see Ch3_lambda_fit.m, Ch3_Ub_fit.m, an_parm_tests.m

clear
% cd 'C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch3\other'
cd 'C:\Users\thasse\Documents\MATLAB\test\Ch3 diss'
load params_meander.mat
% load 'Parametric testing Do.mat'
load 'Parametric testing Dch.mat'
cur_dir = pwd;

g = 9.81;
Cfo = [ my_data.Cfo];
% depth of hypothetical 'straight' reach
Ho = (Qo/2/B)^(2/3)*(Cfo/g/S_valo).^(1/3); 
% must add +2 because Schwenk adds +1 in flowfield and should have -1.
% This is based on an error in Parker & Johannesson 1985 equations 14 & 15
A_effective = [ my_data.A_input ] + 2;     
%% set the models with recursive planform
% ******* this is hard coded and should be updated for other runs ********
% Random initial planform ************
rip = Ch3_rand_initial_planform(Cfo,A_effective);

secs_yr = 60*60*24*365.25;
a1_b1 = 1.9/0.7; % 1 for Do, 1.9/.7 for Dch a1/b1 coefficient
% show the characteristic timescale related to storage for Schwenk & Hasse
cfo=.0036; A = 12; eo = Eo/2; % Schwenk
    % depth of hypothetical 'straight' reach
    Ho = (Qo/2/B)^(2/3)*(cfo/g/S_valo).^(1/3); 
    Ts = a1_b1 * pi*Ho / (sqrt(2) * eo * cfo.^(4/3) * (A-1).^(5/2));
    Ts_Schwenk = Ts/secs_yr
cfo=.024; A = 5; eo = Eo;   % Hasse
    % depth of hypothetical 'straight' reach
    Ho = (Qo/2/B)^(2/3)*(cfo/g/S_valo).^(1/3);
    Ts = a1_b1 * pi*Ho / (sqrt(2) * eo * cfo.^(4/3) * (A-1).^(5/2));
    Ts_Hasse = Ts/secs_yr

ts_step = log10( a1_b1); %a1/b1 ratio for Dch sets the contour spacing
Ts_yr = 10.^sort([2:-ts_step:-0,2+ts_step:ts_step:3.5 ]);
Ts = Ts_yr * secs_yr;
% range of friction factors to use
cfo = 0.001:0.0001:0.06; % cfo = 0.00125:0.0001:0.046;
% depth of hypothetical 'straight' reach (as a vector)
Ho = (Qo/2/B)^(2/3)*(cfo/g/S_valo).^(1/3); 
% A as a function of Ts and cfo, 
Ats = zeros(numel(Ts),numel(cfo));
for i = 1:numel(Ts)
   % equation for Ts but solved for A using dimensional U in chapter 2
    Ats(i,:)= ( a1_b1 * pi*Ho  ./ (Ts(i) * sqrt(2)* Eo * ...
        cfo.^(4/3) ) ).^(2/5) +1; 
end
figure(15)
h=plot(cfo,Ats);
    xlim([.001 .06])
    ylim([1 30])
    box on
    axis square
    xlabel('Friction factor C_f_o','fontsize',12)
    ylabel('Cross stream factor A_e_f_f_e_c_t_i_v_e','fontsize',12)
legend(num2str(round([Ts_yr 1 0 2 Ts_yr])'),'location','southwest')
set(gca,'xscale','log','yscale','log','xgrid','on','ygrid','on')
 
lh = length(h);
for i = 1:lh
    set(h(i),'color',i/1.2./[lh lh lh],'linewidth',mod(i+1,2)+1)
end
hold on
    scatter(Cfo(rip), A_effective(rip), 50,'ko') % succsessful sims
    scatter(Cfo(~rip),A_effective(~rip),50,'kx') % failed sims
    scatter([.0036 .024],[12 5],'ko','filled')   % Schwenk & Hasse sims
hold off
set(gcf,'color','w',...
        'position',[200 200 752 564])
set(gca,'xticklabels',[0.001 0.01],...
    'yticklabels',[1 10])

lab = {'T_1','T_2','T_3','T_4','T_5','T_6','T_7','T_8'};
legend(lab,'location','southwest')

% display model names and contour labels to the command window
spc = {' '  ,' '  ,' '  ,' '  ,' '  ,' '  ,' '  ,' '  }; % spaces vector   
contour_models = 'contour, sinuous 1.85, sinuous 3.7, linear 3.7'
contour_labels = [char(lab),char(spc),...
    strcat(num2str(round(Ts_yr*2    )'),' yr '),char(spc),...
    strcat(num2str(round(Ts_yr      )'),' yr '),char(spc),...
    strcat(num2str(round(Ts_yr/a1_b1*10)'/10),' yr ')]
title('See command window for contour labels, press any key to continue')
disp('Please press any key')
waitforbuttonpress
%% modification for legend vs contour labels
for i = 1:lh
    set(h(i),'color','k','linewidth',1)
end
title('')
legend('Contour labels are T_s (years)',' for models:',...
    'Sinuous E_o = 1.85x10^-^8','Sinuous E_o = 3.7x10^-^8',...
    'Linear E_o = 3.7x10^-^8','location','southwest')

%% save figure
% must save within figure menu due to scatter plot
% fig2svg('Ts contours U contour labels.svg') 
print('-painters', '-dpng', '-r600','Ts contours U contour labels.png')

