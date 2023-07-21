function [fobj, Rsquared] = Ch3_lambda_fit(algorithm,no_coeffs)
% Purpose:  Run fitting algorithms on the output data created for chapter 2
%           of Tobias Hasse's dissertation
% Author:   Tobias Hasse, tobiack@udel.edu
% Date:     June 8, 2021
% Inputs:   algorithm   ('linear' or 'sinuous')
%           no_coeffs   the number of adjustable coefficients to use in the
%                       fit model ( 0 to 3 )
% Outputs:  fobj        A fit object with the equation, coefficeints, etc
%           Rsquared    The R squared of the fit
% also see Ch3_lambda_fit.m, Ch3_Ub_fit.m, an_parm_tests.m

% load params_TRH.mat
switch algorithm
    case 'linear'
        load 'Parametric testing Do.mat'
        [Cfo, A_effective, Ho, lam_obs, B, piHonC]=calc_variables(my_data);
        lam_pred = 1.2*2* piHonC ./ (2 * (A_effective -1) ) .^ .5/(2*B);
    case 'sinuous'
        load 'Parametric testing Dch.mat'
        [Cfo, A_effective, Ho, lam_obs, B, piHonC]=calc_variables(my_data);
        lam_pred = 1.9*2* piHonC ./ (2 * (A_effective -1) ) .^ .5/(2*B);
    otherwise
        disp(strcat('Please choose the "linear" or "sinuous"',...
            ' algorithm. ~Ch3_lambda_fit.m'))
end

%recursive initial planform (HARD CODED)
rip = Ch3_rand_initial_planform(Cfo,A_effective); 
%% AUTOMATIC lambda best fit copied from curve_fits_array.m which was used
% to fit the storage time distributions
switch no_coeffs
    case 3 % best fit all 3 coefficients
        % x is A_effective, y is Cfo
        ft = fittype(' a * 2* piHonC ./ ((2 * (A_effective -b) ) .^ c)',...
            'independent',{'piHonC','A_effective'},'dependent','z'); 
        %      a b c
        sp = [ 1.2 1 .5];  % fit bounds fit all 3 Do & Dch
        up = [20  5  2];
        lo = [ 0 .5 .25];
    case 2 % fit 2 coefficients
        % x is A_effective, y is Cfo
        ft = fittype(' a * 2* piHonC ./ ((2 * (A_effective -1) ) .^ c)',...
            'independent',{'piHonC','A_effective'},'dependent','z'); 
        %      a    c
        sp = [ 1.2  .5];  % fit bounds A-1 Do & Dch
        up = [20    2];
        lo = [ 0  .25];
    case 1 % fit 1 coefficient
        % x is A_effective, y is Cfo
        ft = fittype(' a * 2* piHonC ./ ((2 * (A_effective -1) ).^ .5)',...
            'independent',{'piHonC','A_effective'},'dependent','z');
        %      a 
        sp = [ 1.2  ];  % fit bounds only scalar Do & Dch
        up = [20    ];
        lo = [ 0  ];
    case 0 % use bespoke model
        Rsquared = my_r_squared(lam_obs(rip),lam_pred(rip),0,0);
        fobj = 'none';
        return
    otherwise
        disp(strcat('Please choose 0 to 3 coefficients for the fit',...
            ' model ~Ch3_lambda_fit.m'))
end
[fobj,gof,out] = fit([piHonC(rip)',...
    (A_effective(rip))'],lam_obs(rip)'*2*B,ft,...
    'Startpoint',sp,'lower',lo,'upper',up,'tolfun',1e-10);
    %,'maxfunevals',2000);
% fobj
Rsquared=gof.rsquare;
    
end
function [Cfo, A_effective, Ho, lam_obs, B, piHonC]=calc_variables(my_data)
% This is a subfunction of Ch3_lambda_fit.m
% keyboard
% load params_TRH.mat
load params_meander.mat
g = 9.81;
Cfo = [ my_data.Cfo];
% depth of hypothetical 'straight' reach
Ho = (Qo/2/B)^(2/3)*(Cfo/g/S_valo).^(1/3); 
% must add +2 because Schwenk adds +1 in flowfield and should have -1.
% This is based on an error in Parker & Johannesson 1985 equations 14 & 15
A_effective = [ my_data.A_input ] + 2;     
piHonC = pi* Ho./Cfo; % combine some terms
lam_obs = my_data.lam_crtsn;
end
 
