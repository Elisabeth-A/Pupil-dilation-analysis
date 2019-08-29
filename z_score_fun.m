% z-score calculation

% input arguements:
% B...matrix containing raw intensity values of one ROI colums =
% trials, rows=frames
% f_base... vector containing start and end frame for calculation
%   
% output:
% 

function [M_ROI_z,ROIs_z]=z_score_fun(B,f_base) 
    
B_base=B(f_base(1):f_base(2),:);
    n_frames=length(B);
    mu=mean(B_base);
    S=std(B_base);
    for c=1:size(B,2);
        for i=1:n_frames;
        z=(B(i,c)-mu(:,c))/S(:,c);
         ROIs_z(i,c)=z;   
        end
        M_ROI_z = nanmean(ROIs_z, 2);
    end
    
end