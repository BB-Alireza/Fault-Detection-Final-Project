clc;
clear;
close all;

%% Import the scada_data
[~, ~, scada_data] = xlsread('C:\Users\Javad Hasanpour\Desktop\FAULT project\wind turbinr fault detection\kagel\code\scada.xlsx','scada_data','A2:BN49028');
scada_data(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),scada_data)) = {''};

%% Import the scada_fault
[~, ~, scada_fault] = xlsread('C:\Users\Javad Hasanpour\Desktop\FAULT project\wind turbinr fault detection\kagel\code\scada.xlsx','fault_data','A2:D554');
scada_fault(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),scada_fault)) = {''};

%%
[n_data ,p] = size(scada_data);
[n_fault,~] = size(scada_fault);
scada_data(1:end,67) = mat2cell(6,1,1) ;

for i=1:n_data
    k = 0;                                      % bazy dataha 2 ya 3 fault hamzaman daran
    for j=1:n_fault
        if strcmp(scada_data(i,1),scada_fault(j,1)) & (cell2mat(scada_data(i,2))==cell2mat(scada_fault(j,2)))
            if k==0
                 scada_data(i,67) = scada_fault(j,4);
                 k = k+1;
            else
                 scada_data = [scada_data;scada_data(i,:)];     %repeat data that it has 2 or more fault
                 [new_n_data,~] = size(scada_data);                
                 scada_data(new_n_data,67) = scada_fault(j,4); 
                 k = k+1;
            end
        end
    end
end
            
%% check
[n_data ,p] = size(scada_data);
h=0;
for i=1:n_data
    if cell2mat(scada_data(i,67))==4
        h=h+1;
    end
end

%% to .mat
data = scada_data(1:end,4:end);
data = cell2mat(data);

save('data.mat','data');




