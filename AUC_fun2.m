%% calculation of response integral using trapezoid function and getting rid of normalisation by multiplying with 1/frame_rate


% input:
% df_f...cell array each cell containing df/f values of one ROI
% tone_start...stimulus start [frame number]
% tone_stop...stimulus end [frame number]
% frame_rate...for unnormalise it
%
% output:
% AUC...array with AUC values of all trials (rows) of all ROIs (columns)
% 




function AUC=AUC_fun2(df_f,tone_start, tone_stop,frame_rate)
       if isempty(df_f)
           AUC=NaN;
       else
      for i=1:length(df_f);
        B= df_f{:,i};
        B_tone=B(tone_start:tone_stop,:);
        Q_ROI= trapz(B_tone);
        AUC(:,i)=(Q_ROI');
      end
     AUC=AUC*(1/frame_rate);
       end
end
 