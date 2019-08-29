%% calculating df/f
%
% input arguments:
% B...matrix containing raw intensity values of one ROI columns =
% trials, rows=frames
% f_base...vector with start and end frame for baseline 
 

% output
% matrix with df/f values of one ROI


function y= df_f_fun(B,f_base)
    n_frames=length(B);
    B_base=B(f_base(1):f_base(2),:);
        f0=mean(B_base,1);
        for c=1:size(B,2);
            for i=1:n_frames;
               df_f=(B(i,c)-f0(:,c))/f0(:,c);
                y(i,c)=df_f;   
            end
        end
end