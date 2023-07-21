% analyze length and timescales in 'my_data' for the parametric tests run
% using the sinuous method of Schwenk to calculate excess velocity Ub (Dch)
% and using the more typical hypothetical straightened reach parameters
% (Do)
% Tobias Hasse July 24, 2020
% *************Figure 2 chapter 2*****************
% also see Ch3_lambda_fit.m, Ch3_Ub_fit.m, an_parm_tests.m

clear
% cd 'C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch3\other'
cd 'C:\Users\thasse\Documents\MATLAB\test\Ch3 diss'

load params_meander.mat
load 'Parametric testing Do.mat'

cur_dir = pwd;
%%
g = 9.81;
Cfo = [ my_data.Cfo];
% depth of hypothetical 'straight' reach
Ho = (Qo/2/B)^(2/3)*(Cfo/g/S_valo).^(1/3); 
% must add +2 because Schwenk adds +1 in flowfield and should have -1.
% This is based on an error in Parker & Johannesson 1985 equations 14 & 15
A_effective = [ my_data.A_input ] + 2;     

lam_obs = my_data.lam_crtsn;

%% set the models with recursive planform
% ******* this is hard coded and should be updated for other runs ********
rip = Ch3_rand_initial_planform(Cfo,A_effective);
%% ************* Enter fit equation here ********************
lam_pred = 1.2 * 2 * pi * Ho ./ (Cfo .* (2 * (A_effective -1 ) )...
    .^ 0.5) / (2*B);
lmo = lam_obs(rip);
lmp = lam_pred(rip);
aef = A_effective(rip);
ua = unique(aef);
% figure for lambda
figure(13)
clf
for plt = 1:2
    if isequal(plt,2)
        % if the next file is in a different directory
%         cd 'C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch3\other'
        cd 'C:\Users\thasse\Documents\MATLAB\test\Ch3 diss'
        load 'Parametric testing Dch.mat'
        cd(cur_dir)
        Cfo = [ my_data.Cfo];
        Ho = (Qo/2/B)^(2/3)*(Cfo/g/S_valo).^(1/3); 
        A_effective = [ my_data.A_input ] + 2;     
        lam_obs = my_data.lam_crtsn;
        % set the models with recursive planform
        % **** this is hard coded and should be updated for other runs ****
        % sims with random initial planform 
        rip = Ch3_rand_initial_planform(Cfo,A_effective);
% ************* Enter fit equation here ********************
        lam_pred = 1.9 * 2 * pi * Ho ./ (Cfo .* (2 * (A_effective -1 ) )...
            .^ 0.5) / (2*B);
        lmo = lam_obs(rip);
        lmp = lam_pred(rip);
        aef = A_effective(rip);
        ua = unique(aef);
    end
subplot(1,2,plt)
% separate the data into vectors for each value of A_effective, and plot
for i = 1 : length(ua)
    arr = [];
    v1 = false(size(aef));
    v1(aef == ua(i)) = true;
    arr(1,:)= lmo(v1);
    arr(2,:)= lmp(v1);
    myplot(i).data = arr;
    hold on
    h(i) = plot(arr(1,:),arr(2,:));
end
hold off
% set marker and linestyle
mk_sz=[6 8 6 6 4]; 
for i = 1:length(h)
   set(h(i),'linestyle','none','color','k')
   if i == 1
       set(h(i),'marker','^','markerfacecolor','k','markersize',mk_sz(i))
   end
   if i == 2
       set(h(i),'marker','>','markerfacecolor',[.8 .8 .8],...
           'markeredgecolor','k','markersize',mk_sz(i))
   end
   if i == 3
       set(h(i),'marker','s','markerfacecolor','k','markersize',mk_sz(i))
   end
   if i == 4
       set(h(i),'marker','d','markerfacecolor',[.8 .8 .8],...
           'markeredgecolor','k','markersize',mk_sz(i))
   end
   if i == 5
       set(h(i),'marker','o','markerfacecolor','k','markersize',mk_sz(i))
   end
   if i == 6
       set(h(i),'marker','<','markerfacecolor','k')
   end
   if i == 7
       set(h(i),'marker','>','markerfacecolor','k')
   end
end

%% plot the diagonal 1:1 line
hold on
plot([0,45],[0,45],'color','k','linestyle','-')
hold off
axis square
box on
ylim([0 45])
xlim([0 45])
%
set(gcf,'color','w',...
    'position',[100 100 700 564])
set(gca,'xtick',[0 10 20 30 40],...
    'ytick',[0 10 20 30 40])
xlabel('Simulated characteristic meander wavelength l_v_c/2B',...
    'fontsize',12)
ylabel(sprintf('Predicted characteristic meander wavelength l_v_c/2B'),...
    'fontsize',12)
lab = strcat('A= ', num2str(ua'));
lab(length(lab)+1,:)='1:1 ';
legend( lab,'location','southeast')
    if isequal(plt, 1)
        text(2,42,'A: Linear Method','fontsize',12)
    else
        text(2,42,'B: Sinuous Method','fontsize',12)
    end
end %plott

%% save figure
% fig2svg('lam obs pred for Do Dch_scale.svg')
print('-painters', '-dpng', '-r600','lam obs pred for Do Dch_scale.png')

