function [Delta_final, Percentage] = Delta_index(DDM, Dose, Structure, Percentage, Per_dose, Dose_tol, range, VoxelsizeX, VoxelsizeY, VoxelsizeZ)

%% Calculate the delta index - 05-03-2021
% 
% [Delta_total, Percentage] = Delta index(...) calculates the 3D matrix 
% with delta index values and the delta index passing rate for the given 
% inputs. 

% Input: 
%   - DDM = 3D matrix with DDM values
%   - Dose = 3D matrix with Dose values (size must match to DDM matrix)
%   - Structure = 3D matrix with the structure with values 1 and the part
%   outside the structure without values (NaN). Size must match to the DDM
%   and dose matrices. You can use this structure if you want to calculate 
%   the delta over a specified cropped region. If you enter an empty struct
%   [] the full 3D matrix range is used.
%   - Percentage = threshold percentages (commonly 20, 50 and 90% are used)
%   - Per_dose = Perscription dose (in Gy)
%   - Dose_tol = Dose tolerance criteria (mostly 3 or 2%)
%   - Range = To save calculation time only calculate the delta index in a
%   range of DDM (commonly 3x DDM is used)
%   - VoxelsizeX, VoxelsizeY, VoxelsizeZ = Contain the voxel resolution
%   sizes in mm  (example voxel sizes: 0.5872  0.5872  1.2000)
%
% Output: 
%   - Delta_total = delta index values (in 3D matrix)
%   - Percentage = delta index passing rate (in %)


% Dose is defined as a minimum threshold dose dependent on the prescription
% dose and the percentage 
Dose(Dose < (Per_dose*(Percentage*0.01))) = 0;

% Make sure the Dose, DDM and structure matrix sizes are the same size, 
% otherwise throw an error.

if ~isempty(Structure)
assert(isequal(size(DDM),size(Dose), size(Structure)),'Dose, DDM and/or structures have different matrix sizes')
end

assert(isequal(size(DDM),size(Dose)),'Dose and DDM have different matrix sizes')

Size_matrix = size(DDM);
Delta_total = zeros([Size_matrix(1) Size_matrix(2) Size_matrix(3)]); % Zero matrix

% Calculate the delta index for every voxel in the matrix
for x_ref = 1:Size_matrix(1)
    for y_ref = 1:Size_matrix(2) 
        for z_ref = 1:Size_matrix(3)
            Delta_total(x_ref,y_ref,z_ref) = NaN;
            DDM_Ref = DDM(x_ref,y_ref,z_ref);
            Dose_Ref = Dose(x_ref,y_ref,z_ref);
            % Define the distance of the neigbourhood for every voxel which
            % is used to define the range for the evaluating voxels.
            DDM_dist_voxels = (range*DDM_Ref);
            DDM_dist_voxels_ceilXY = floor((range*DDM_Ref)/VoxelsizeX); % In mm
            DDM_dist_voxels_ceilZ = floor((range*DDM_Ref)/VoxelsizeZ); % In mm
            if Dose_Ref > 0  % Only look at this voxel if the Reference dose is nonzero
                Delta = 0;
                for x_eval = (x_ref-DDM_dist_voxels_ceilXY):(x_ref+DDM_dist_voxels_ceilXY) % get all evaluating voxels in the specified range
                    for y_eval = (y_ref-DDM_dist_voxels_ceilXY):(y_ref+DDM_dist_voxels_ceilXY)
                        for z_eval = (z_ref-DDM_dist_voxels_ceilZ):(z_ref+DDM_dist_voxels_ceilZ)
                    % The voxels should not go out of range or to negative
                    % voxels, therefore the evaluating voxels on the edge
                    % of the matrix are adapted.
                        if x_eval < 1
                        x_eval = 1;
                        end
                        if y_eval < 1
                        y_eval = 1;
                        end
                        if z_eval < 1
                        z_eval = 1;
                        end
                        if x_eval > Size_matrix(1)
                        x_eval = Size_matrix(1);
                        end
                        if y_eval > Size_matrix(2)
                        y_eval = Size_matrix(2);
                        end
                        if z_eval > Size_matrix(3)
                        z_eval = Size_matrix(3);
                        end
                        % Determine distance between reference and evaluating voxel
                        Distance_Ref_Ev = sqrt(((x_eval-x_ref)*VoxelsizeX)^2 + ((y_eval-y_ref)*VoxelsizeY)^2+ ((z_eval-z_ref)*VoxelsizeZ)^2); 
                        % The distance should be smaller than the DDM*range and the dose in the evaluating voxel should be higher than 0. 
                        if Distance_Ref_Ev <= DDM_dist_voxels && Dose(x_eval,y_eval,z_eval)>0  
                        % Calculate the delta index for all evaluating
                        % voxels; save as temporary delta
                        temp_delta = abs(((Dose_Ref - Dose(x_eval,y_eval,z_eval))/(Dose_Ref*Dose_tol)))*(exp(-0.5*((Distance_Ref_Ev)^2/DDM_Ref^2)));
                            if isnan(temp_delta)
                            temp_delta = 0;
                            end 
                            % Select the largest value over all evaluating voxels
                            if temp_delta > Delta 
                            Delta = temp_delta;
                            end
                        end
                    end
                end
            end
            if ~isempty(Delta)
                Delta_total(x_ref,y_ref,z_ref) = max(Delta);
            end
            clear Delta
         end
    end
end
end

if ~isempty(Structure)
Delta_final = Structure.*Delta_total(:,:,1:Size_matrix(3));
elseif isempty(Structure)
Delta_final = Delta_total(:,:,1:Size_matrix(3));
end

% Calculate the delta index passing rate
Percentage = 100 - length(find(Delta_final>1))/length(find(Delta_final>=0)) * 100;

end

