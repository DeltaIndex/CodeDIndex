
% Delta index calculation example 
clear all;
close all;

% Expected percentage result
% Threshold dose 20%: 80%
% Threshold dose 50%: 57%
% Threshold dose 90%: 72%

% Open Example files
load('Example_Files/Dose.mat');
load('Example_Files/DDM.mat');
load('Example_Files/Structure.mat');
load('Example_Files/Voxelsize.mat');

% Function to calculate the Delta Index
[Delta_total_20, Percentage_20] = Delta_index(DDM, Dose,...
    Structure, 20, 60, 0.03, 3, Voxelsize(1), ...
    Voxelsize(2), Voxelsize(3));
[Delta_total_50, Percentage_50] = Delta_index(DDM, Dose,...
    Structure, 50, 60, 0.03, 3, Voxelsize(1), ...
    Voxelsize(2), Voxelsize(3));
[Delta_total_90, Percentage_90] = Delta_index(DDM, Dose,...
    Structure, 90, 60, 0.03, 3, Voxelsize(1), ...
    Voxelsize(2), Voxelsize(3));

% Example plot
figure(); 
imagesc(Delta_total_20(:,:,ceil(size(Delta_total_50,3)/2)));
title('Delta index for prescription dose higher than 20% of the dose');
colorbar;
