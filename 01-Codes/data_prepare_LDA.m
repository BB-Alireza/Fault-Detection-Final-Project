
%% ---------------------------------- Load data -------------------------------------------
clc;
clear;
load('97.mat' );
load('105.mat');
load('118.mat');
load('130.mat');
N  = X097_DE_time ;     %class 1
F1 = X105_DE_time ;     %class 2  
F2 = X130_DE_time ;     %class 3
F3 = X118_DE_time ;     %class 4

hold on;
plot(X097_DE_time(1:200))     %200 bode mishe
plot(X097_DE_time(201:400))
plot(X097_DE_time(401:600))

%% --------------------------------- create data -----------------------------------------------
Normal=[];
Fault1=[];
Fault2=[];
Fault3=[];
for i=1:100
    Normal = [Normal;N(1:200)'] ;
    N(1:200)=[];
end
for i=1:100
    Fault1 = [Fault1;F1(1:200)'] ;
    F1(1:200)=[];
end
for i=1:100
    Fault2 = [Fault2;F2(1:200)'] ;
    F2(1:200)=[];
end
for i=1:100
    Fault3 = [Fault3;F3(1:200)'] ;
    F3(1:200)=[];
end
%Normal = Normal' ;
%Fault1 = Fault1' ;
%Fault2 = Fault2' ;
%Fault3 = Fault3' ;


%% ------------------------------ Feature Extraction -------------------------------------------
Normal = [Normal,Feature_extraction(Normal)];
Fault1 = [Fault1,Feature_extraction(Fault1)];
Fault2 = [Fault2,Feature_extraction(Fault2)];
Fault3 = [Fault3,Feature_extraction(Fault3)];


%% ---------------------------------- START LDA ----------------------------------------------

xf1 = Fault1' ;
xf2 = Fault2' ;
xf3 = Fault3' ;
xn  = Normal' ;

[p,lenN ] = size(xn) ;
[~,lenF1] = size(xf1);
[~,lenF2] = size(xf2);
[~,lenF3] = size(xf3);
x = [xf1 , xf2 , xf3 , xn] ;
num_class = 4 ;                     % number of class 

%% ------------------------------------ LDA -----------------------------------------------------

mean_n  = mean(xn ,2);                 % mean normal
mean_f1 = mean(xf1,2);                 % mean fault1
mean_f2 = mean(xf2,2);                 % mean fault2
mean_f3 = mean(xf3,2);                 % mean fault3


%S1 = zeros(p,p);                         %bedone estefade az cov(x)         
%S2 = zeros(p,p);
%S3 = zeros(p,p);
%for i=1:num_n
%    S1 = S1 + (xn(:,i)-mean_n)*(xn(:,i)-mean_n)' ;
%end
%for i=1:num_f1
%    S2 = S2 + (xf1(:,i)-mean_f1)*(xf1(:,i)-mean_f1)' ;
%end
%for i=1:num_f2
%    S3 = S3 + (xf2(:,i)-mean_f2)*(xf2(:,i)-mean_f2)' ;
%end

S1 = cov(xn' );
S2 = cov(xf1');
S3 = cov(xf2');
S4 = cov(xf3');
Sw = S1 + S2 + S3 + S4 ;


Mean = mean(x,2);
SB1 = lenN  * (mean_n  - Mean)*(mean_n  - Mean)' ;
SB2 = lenF1 * (mean_f1 - Mean)*(mean_f1 - Mean)' ;
SB3 = lenF2 * (mean_f2 - Mean)*(mean_f2 - Mean)' ;
SB4 = lenF3 * (mean_f3 - Mean)*(mean_f3 - Mean)' ;
SB = SB1 + SB2 + SB3 + SB4 ;


invSw_SB = inv(Sw) * SB ;
[V,D] = eig (invSw_SB) ;
Landa = diag(D);
[Landa,Landa_label]=sort(Landa,'descend');
 
 W=[];
 for i=1:(num_class-1)
    W(:,i)=V(:,Landa_label(i));
 end
y = W' * x ;


%% ------------------------------------ output normalize ----------------------------------------
y=y';
Meany = mean(y) ;
inVary = diag(1 ./ sqrt(var(y)));
p=zeros(size(y));
for i=1:size(y,1)
   p(i,:) = inVary * (y(i,:)' - Meany') ; 
end
y=p;

FAULT1 = y(1:lenF1,:) ;
FAULT2 = y(lenF1+1:lenF1+lenF2,:) ;
FAULT3 = y(lenF1+lenF2+1:lenF1+lenF2+lenF3,:) ;
NORMAL = y(lenF1+lenF2+lenF3+1:lenF1+lenF2+lenF3+lenN ,:) ; 

% 3D plot
figure; 
plot3(FAULT1(:,1),FAULT1(:,2),FAULT1(:,3),'bo','MarkerFaceColor','b'); hold on; grid on;
plot3(FAULT2(:,1),FAULT2(:,2),FAULT2(:,3),'ro','MarkerFaceColor','r');
plot3(FAULT3(:,1),FAULT3(:,2),FAULT3(:,3),'yo','MarkerFaceColor','y');
plot3(NORMAL(:,1),NORMAL(:,2),NORMAL(:,3),'go','MarkerFaceColor','g');
title('LDA'); legend('FAULT1','FAULT2','FAULT3','NORMAL');


%% -------------------------- tarkib & taghsim test & train ------------------------------------
n=75 ;                                    %darsade train
[F1_tr,F1_ts] = test_train_DataGenerator (FAULT1,n);
[F2_tr,F2_ts] = test_train_DataGenerator (FAULT2,n);
[F3_tr,F3_ts] = test_train_DataGenerator (FAULT3,n);
[N_tr ,N_ts ] = test_train_DataGenerator (NORMAL,n);

%% ------------------------------ save data mat format -------------------------------------------
save('data.mat','F1_tr','F1_ts','F2_tr','F2_ts','F3_tr','F3_ts','N_tr','N_ts');


    
    
    
    