% analyze length and timescales in 'my_data' for the parametric tests run
% using the sinuous method of Schwenk to calculate excess velocity Ub (Dch)
% and using the more typical hypothetical straightened reach parameters
% (Do)
% Tobias Hasse July 24, 2020
% Tobias Hasse Dec 4, 2020 FIXING to use dimensional value
% *************Figure 3 chapter 2*****************
% also see Ch3_lambda_fit.m, Ch3_Ub_fit.m, an_parm_tests.m
clear
% cd 'C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch3\other'
cd 'C:\Users\thasse\Documents\MATLAB\test\Ch3 diss'
load params_meander.mat
load 'Parametric testing Do.mat'

cur_dir = pwd;

g = 9.81;
Cfo = [ my_data.Cfo];
% depth of hypothetical 'straight' reach
Ho = (Qo/2/B)^(2/3)*(Cfo/g/S_valo).^(1/3); 
% must add +2 because Schwenk adds +1 in flowfield and should have -1.
% This is based on an error in Parker & Johannesson 1985 equations 14 & 15
A_effective = [ my_data.A_input ] + 2;     

% second version dimensionalizing U
ub_obs = [my_data.Ub_median].*[my_data.U_mean]*B;   

%% set the models with recursive planform
% ******* this is hard coded and should be updated for other runs ********
rip = Ch3_rand_initial_planform(Cfo,A_effective);
%% quick plot of the data
% figure(4)
%     scatter(Cfo,A_effective,my_data.Ub_median*1000,...
%         log(my_data.Ub_median),'filled')
%     colormap(jet)
%     lab = num2str(~rip');%[A_effective]');%Ub_median]');
%     lab = cellstr(lab);%(:,2:4));
%     text(Cfo,A_effective,lab)
%     set(gca,'xscale','log','yscale','log')
%% Pretty figure, 
% ************* Enter fit equation here ********************
% ub_pred = 0.2* (A_effective-1).^1.5 .* Cfo.^(2/3);
ub_pred = 1.2* (A_effective-1).^1.5 .* Cfo.^(1/3);

% split up data by A_effective value
ubo = ub_obs(rip);
ubp = ub_pred(rip);
aef = A_effective(rip);
ua = unique(aef);
% figure for ub
figure(14)
clf
for plt = 1:2
    if isequal(plt,2)
%         cd 'C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch3\other'
        cd 'C:\Users\thasse\Documents\MATLAB\test\Ch3 diss'
        load 'Parametric testing Dch.mat'
        cd(cur_dir)
        Cfo = [ my_data.Cfo];
        Ho = (Qo/2/B)^(2/3)*(Cfo/g/S_valo).^(1/3); 
        A_effective = [ my_data.A_input ] + 2;     
        % second version dimensionalizing U
        ub_obs = [my_data.Ub_median].*[my_data.U_mean]*B;
        % set the models with recursive planform
        % **** this is hard coded and should be updated for other runs ****
        % sims with random initial planform 
        rip = Ch3_rand_initial_planform(Cfo,A_effective);
% ************* Enter fit equation here ********************
        ub_pred = 0.7* (A_effective-1).^1.5 .* Cfo.^(1/3);
        ubo = ub_obs(rip);
        ubp = ub_pred(rip);
        aef = A_effective(rip);
        ua = unique(aef);
    end
subplot(1,2,plt)
for i = 1 : length(ua)
    arr = [];
    v1 = false(size(aef));
    v1(aef == ua(i)) = true;
    arr(1,:)= ubo(v1);
    arr(2,:)= ubp(v1);
    myplot(i).data = arr;
    hold on
    h(i) = plot(arr(1,:),arr(2,:));
end
hold off
%
mk_sz=[5 4 6 4 3]; % original sizes [6 5 7 5 4]
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
       set(h(i),'marker','o','markerfacecolor','k',...
           'markeredgecolor','k','markersize',mk_sz(i))
   end
   if i == 6
       set(h(i),'marker','<','markerfacecolor','k')
   end
   if i == 7
       set(h(i),'marker','>','markerfacecolor','k')
   end
end
hold on
plot([.003,40],[.003,40],'color','k','linestyle','-')
hold off
axis square
box on
ylim([.1 40])
xlim([.1 40])

set(gca,'xscale','log','yscale','log','fontsize',8)
set(gcf,'color','w',...
        'position',[100 100 700 564]) 
set(gca,'xtick',[.1 1 10],...
    'xticklabels',[.1 1 10],...
    'ytick',[.1 1 10],...
    'yticklabels',[.1 1 10])
xlabel('Simulated characteristic excess velocity','fontsize',12)
ylabel('Predicted characteristic excess velocity','fontsize',12)
lab = strcat('A= ', num2str(ua'));
lab(length(lab)+1,:)='1:1 ';
legend( lab,'location','southeast')
    if isequal(plt, 1)
        text(.15,25,'A: Linear Method','fontsize',12)
    else
        text(.15,25,'B: Sinuous Method','fontsize',12)
    end
end %plt
%% save figure
% fig2svg('Ub obs pred Do Dch dimensioned.svg')
print('-painters', '-dpng', '-r600','Ub obs pred for Do Dch_scale.png')

