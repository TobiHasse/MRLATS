% Script for checking parameters in simulation output files from Chapter 3
% Written by Tobias Hasse, June 2021

% This is just the skelaton, if you want to inspect a certain variable from
% each of the files this script will open each file and grab that variable
% allowing for a quick comparison

clear
% cd(strcat('C:\Users\User\Documents\MATLAB\Heap\Parametric testing',...
%     ' Do\Do Ub median slope 000 67'))
% %     ' Dch\Dch Ub median slope 000 67'))
cd 'C:\Users\thasse\Documents\MATLAB\test\Ch3 diss'
files = dir;
for i = 1:numel(files)
    infilename = files(i).name;
    if length(infilename)>10 % only check file names long enough (avoid err
        if strcmp(infilename(end-9:end),'_other.mat')
            load(infilename)
            % a series of if statements to check different things
            if Eo == 1.85*10^-8
%             if scaleup == 5.7
%             if strcmp(inplanformname,'initialplanform')

                % do something if true
%                 sprintf('regular Eo')
            else 
%                 sprintf('Eo has been changed to %1.2e',Eo)
                infilename
%                 inplanformname
%                 keyboard
            end
            
        end
    end
end  

