%%% Pupil dilation analysis main %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

addpath(genpath('/Users/')); 

mouse_ID={}

% load file physiological data for info on trial identities
file_phys=char(strcat(mouse_ID,'_raw.mat'));

load(file_phys,'n_trials_hab', 'n_trials_ret',  'minus_trials', 'plus_trials', 'baseline_trials',...
                'minus_trials_ret', 'plus_trials_ret', 'baseline_trials_ret');

% set paramaters
bin_factor=2
n_frames=900
n_frames_ret=700
trial_length= 35

% calculate frame_rate
frame_rate=n_frames/trial_length;
frame_rate_ret=n_frames_ret/trial_length_ret;

% minimum number of frames a trial has to have to be included
N_values= ceil(n_frames/bin_factor)-2;
N_values_ret= ceil(n_frames_ret/bin_factor)-2;

% determine baseline 
d_base=[floor((5*frame_rate)/bin_factor), floor((20*frame_rate)/bin_factor)];
d_base_recall=[floor((5*frame_rate_ret)/bin_factor), floor((20*frame_rate_ret)/bin_factor)];

% determine end of integration window for response integral 
window_AUC_end_ret=ceil((30*frame_rate_ret)/bin_factor);
window_AUC_end=ceil((30*frame_rate)/bin_factor);

%% load pupil dilation .txt files
% habituation
 baseline_hab=[];
 CS_minus_hab=[];
 CS_plus_hab=[];
 S=[];
 
for i=1:n_trials_hab
    filename=char(strcat('/SOMATA/pupil_dilation/',mouse_ID,'/hab/',num2str(i),'.txt')); %intensitiy values rows = ROIs columns=frames
     if exist(filename) == 0;
        empty= NaN(n_frames,1);
        raw_hab=empty;
    else
    load(filename) 
    fRead = fopen(filename,'r');
    data = textscan(fRead, '%n %n %n %n %n','delimiter', '\t'); % read the whole file in the data variable and removes header (check help textscan)
    fclose(fRead);
    % check if too noisy -> not more than a certain number of frames in a
      if length(data{5})<(n_frames-4)
         empty= NaN(n_frames,1);
         data_in(:,i)=empty;
         A={'trial too short'}
      elseif length(data{5})>(n_frames)
          data_in(:,i)=[data{5}(1:n_frames)];
          length(data{5})
     elseif length(data{5})<(n_frames) && length(data{5})>(n_frames-4)
         data_NaN=NaN(n_frames-length(data{5}),1);
         data_in(:,i)=[data{5};data_NaN];
     else
        data_in(:,i)=data{5};
      end
    % smooth data use 10% of data points
    if isnan(data_in(:,i))
        S(:,i)=NaN(n_frames,1);
    else
        S(:,i)=smooth(data_in(:,i),0.05,'rloess');
     end
   end
     clear data
end 


for j=1:n_trials_hab     
   B=binning(S(:,j),[bin_factor,1],'mean');
   if any(baseline_trials==j)
        baseline_hab(:,end+1)=B(1:N_values);
        elseif any(minus_trials==j)
          CS_minus_hab(:,end+1)=B(1:N_values);
        elseif any(plus_trials==j)
        CS_plus_hab(:,end+1)=B(1:N_values);
   end
   clear B  


end      
 


%% retrieval
 baseline_ret=[];
 CS_minus_ret=[];
 CS_plus_ret=[];
 S_ret=[];
for i=1:n_trials_ret
    filename=char(strcat('SOMATA/pupil_dilation/',mouse_ID,'/24h/',num2str(i),'.txt')); %intensitiy values rows = ROIs columns=frames
    if exist(filename) == 0;
        empty= NaN(N_values_ret,1);
        data{5}=empty;
    else
        load(filename) 
        fRead = fopen(filename,'r');
        data = textscan(fRead, '%n %n %n %n %n','delimiter', '\t'); % read the whole file in the data variable and removes header (check help textscan)
        fclose(fRead);
     % check if too noisy -> not more than a certain number of frames in a
     if length(data{5})<(n_frames_ret-4)
         empty= NaN(n_frames_ret,1);
         data_in_ret(:,i)=empty;
     elseif length(data{5})>(n_frames_ret)
          data_in_ret(:,i)=[data{5}(1:n_frames_ret)];
          length(data{5})
     elseif length(data{5})<(n_frames_ret) && length(data{5})>(n_frames_ret-4)
         data_NaN=NaN(n_frames_ret-length(data{5}),1);
         data_in_ret(:,i)=[data{5};data_NaN];
     else
        data_in_ret(:,i)=data{5};
     end
    % smooth data use 10% of data points
        if isnan(data_in_ret(:,i))
            S_ret(:,i)=NaN(n_frames_ret,1);
        else
            S_ret(:,i)=smooth(data_in_ret(:,i),0.05,'rloess');
        end 
    end 

end


for j=1:n_trials_ret    
    B=binning(S_ret(:,j),[bin_factor,1],'mean');
    if any(baseline_trials_ret==j)
        baseline_ret(:,end+1)=B(1:N_values_ret);
    elseif any(minus_trials_ret==j)
          CS_minus_ret(:,end+1)=B(1:N_values_ret);
    elseif any(plus_trials_ret==j)
        CS_plus_ret(:,end+1)=B(1:N_values_ret);
    end  
    clear B  
end
% Calculate approx. timing
T=0;
f=(1/frame_rate);
for n= 1:(n_frames/bin_factor)
     t= n*f;
     T(n,:)= t;
end


%% calculate change in pupil dilation during CS presentation
%habituation
dd_d_minus=df_f_fun(CS_minus_hab,d_base);
dd_d_plus=df_f_fun(CS_plus_hab,d_base);
dd_d_baseline=df_f_fun(baseline_hab,[1,length(baseline_hab)]);

%retrieval
dd_d_minus_ret=df_f_fun(CS_minus_ret,d_base_recall);
dd_d_plus_ret=df_f_fun(CS_plus_ret,d_base_recall);
dd_d_baseline_ret=df_f_fun(baseline_ret,[1,length(baseline_ret)]);


%% calculate z-score in pupil dilation during CS presentation
%habituation
[mz_minus, z_minus]=z_score_fun(CS_minus_hab,d_base);
[mz_plus, z_plus]=z_score_fun(CS_plus_hab,d_base);
[mz_baseline, z_baseline]=z_score_fun(baseline_hab,[1,length(baseline_hab)]);

%retrieval
[mz_minus_ret,z_minus_ret]=z_score_fun(CS_minus_ret,d_base_recall);
[mz_plus_ret,z_plus_ret]=z_score_fun(CS_plus_ret,d_base_recall);
[mz_baseline_ret,z_baseline_ret]=z_score_fun(baseline_ret,[1,length(baseline_ret)]);
%% AUC
AUC_plus=AUC_fun2({dd_d_plus},d_base(2)+1,window_AUC_end, frame_rate);
AUC_minus=AUC_fun2({dd_d_minus},d_base(2)+1,window_AUC_end, frame_rate);

AUC_plus_ret=AUC_fun2({dd_d_plus_ret},d_base_recall(2)+1,window_AUC_end_ret, frame_rate_ret);
AUC_minus_ret=AUC_fun2({dd_d_minus_ret},d_base_recall(2)+1,window_AUC_end_ret, frame_rate_ret);

% AUC cell array cell 1=minus hab, 2=plus hab, 3=minus ret, 4=plus ret
AUC_pupil={{AUC_plus} {AUC_minus}  {} {} {AUC_plus_ret } {AUC_minus_ret} };

%% save data
file_name=char(strcat(mouse_ID,'_pupil_dilation_smoothed.mat'));
save(file_name);