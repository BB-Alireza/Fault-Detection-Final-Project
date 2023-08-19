%% class isolation
clc;
clear;
load('data.mat' );
data(:,34) = [];

data(:,1:end-1) = zscore(data(:,1:end-1));
N  = [] ;
F1 = [] ;
F2 = [] ;
F3 = [] ;
F4 = [] ;
F5 = [] ;



for i=1:length(data)
    switch data(i,end)
        case 1
          F1 = [F1 ; data(i,1:end-1)];  %class 1
        case 2
          F2 = [F2 ; data(i,1:end-1)];  %class 2  
        case 3
          F3 = [F3 ; data(i,1:end-1)];  %class 3
        case 4
          F4 = [F4 ; data(i,1:end-1)];  %class 4
        case 5
          F5 = [F5 ; data(i,1:end-1)];  %class 5
        case 6
          N  = [N  ; data(i,1:end-1)];  %class 6
    end
end

save('data_after_classisolation.mat','F1','F2','F3','F4','F5','N');
