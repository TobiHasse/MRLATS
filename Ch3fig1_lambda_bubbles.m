% analyze length and timescales in 'my_data' for the parametric tests run
% using the sinuous method of Schwenk to calculate excess velocity Ub (Dch)
% and using the more typical hypothetical straightened reach parameters
% (Do)
% Tobias Hasse July 24, 2020
% *************Figure 1 chapter 2*****************
% also see Ch3_lambda_fit.m, Ch3_Ub_fit.m, an_parm_tests.m

% cd 'C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch3\other'
cd 'C:\Users\thasse\Documents\MATLAB\test\Ch3 diss'
for ii = 1:2 % both algorithms
    clearvars -except ii
    load params_meander.mat
    linear_algorithm = 1; % second time, switch to sinuous algorithm
    alg = 'Do';           % first time, use linear parameters algorithm
    load 'Parametric testing Do.mat'
    method = 'Linear';
    outfilename = sprintf('Parameter space %s legend',alg);
    if isequal(ii,2)
        linear_algorithm = 0; % second time, switch to sinuous algorithm
        alg = 'Dch';
        load 'Parametric testing Dch.mat'
        method = 'Sinuous';
        outfilename = sprintf('Parameter space %s legend',alg)
    end


g = 9.81;
Cfo = [ my_data.Cfo];
% depth of hypothetical 'straight' reach
Ho = (Qo/2/B)^(2/3)*(Cfo/g/S_valo).^(1/3); 
% must add +2 because Schwenk adds +1 in flowfield and should have -1.
% This is based on an error in Parker & Johannesson 1985 equations 14 & 15
A_effective = [ my_data.A_input ] + 2;     


%% set the models with recursive planform
% ******* this is hard coded and should be updated for other runs ********
rip = Ch3_rand_initial_planform(Cfo,A_effective);
% figure(4)
%     scatter(Cfo,A_effective,my_data.Ub_median*1000,...
%         log(my_data.Ub_median),'filled')
%     colormap(jet)
%     lab = num2str(recursive_planform');%[A_effective]');%Ub_median]');
%     lab = cellstr(lab);%(:,2:4));
%     text(Cfo,A_effective,lab)
%     set(gca,'xscale','log','yscale','log')
%% scatter plot of lambda
bn = colormap(bone(64));  % choose grayscale colormap

lam_obs = my_data.lam_crtsn;
figure(10+ii)
   scatter(Cfo(rip),A_effective(rip),lam_obs(rip)*40,lam_obs(rip),'filled')
    cx = caxis;
    caxis([-cx(2) cx(2)*1.1])
    hold on
    scatter(Cfo(~rip),A_effective(~rip),50,ones(size(Cfo(~rip)))/.1,'x')
    hold off
    colormap(flipud(bn(32:end,:)))
    set(gca,'xscale','log','yscale','log','fontsize',12)
    set(gcf,'color','w',...
        'position',[100 100 752 564])
    xlim([.001 .06])
    ylim([1 30])
    box on
    axis square
    title(sprintf('%s method, down valley meander wavelength',method))
    xlabel('Friction factor C_f_o','fontsize',12)
    ylabel('Cross stream factor A_e_f_f_e_c_t_i_v_e','fontsize',12)
    legend(sprintf('%s_v_c / 2B',char(955)),'Failed simulation',...
        'location','southwest','fontsize',12)
    labs = lam_obs;
    labs(labs>9.5) = round(labs(labs>9.5),0);
    labs(labs<9.5) = round(labs(labs<9.5),1);
    
    lab = num2str(labs(rip)');%[A_effective]');%Ub_median]');
    lab = cellstr(lab);%(:,2:4));
    text(Cfo(rip)/1.14,A_effective(rip),lab,'fontsize',12)
%%
% fig2svg('Parameter space Do legend.svg')
     print('-painters', '-dpng', '-r600',sprintf('%s.png',outfilename))
end 
