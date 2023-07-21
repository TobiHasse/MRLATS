function [rip] = Ch3_rand_initial_planform(Cfo,A_effective)
% This function is HARD CODED to return a vector with the recursive
% planforms of chapter 2 of Tobias Hasse's dissertation.
% Author: Tobias Hasse tobiac@udel.edu
% June 8 2021

% it would be better to save into the struct "my_data" whether the run was
% recursive
recursive_planform = false(size(Cfo));
recursive_planform(A_effective < 4) = true;
recursive_planform(34)=false;
recursive_planform(33)=true;
rip = ~recursive_planform; % Random initial planform ************

rip = ~~[1,0,1];  %% Hack for mini testing data set 2023

end
