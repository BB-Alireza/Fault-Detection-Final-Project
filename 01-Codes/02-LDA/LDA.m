
%% ---------------------------------- Load data -------------------------------------------
clc;
clear;
load('data_after_classisolation.mat');


%% ---------------------------------- START LDA ----------------------------------------------

xf1 = F1' ;
xf2 = F2' ;
xf3 = F3' ;
xf4 = F4' ;
xf5 = F5' ;
xn  = N'  ;

[p,lenN ] = size(xn) ;
[~,lenF1] = size(xf1);
[~,lenF2] = size(xf2);
[~,lenF3] = size(xf3);
[~,lenF4] = size(xf4);
[~,lenF5] = size(xf5);
x = [xf1 , xf2 , xf3 , xf4 , xf5 , xn] ;
num_class = 6 ;                     % number of class 

%% ------------------------------------ LDA -----------------------------------------------------

mean_n  = mean(xn ,2);                 % mean normal
mean_f1 = mean(xf1,2);                 % mean fault1
mean_f2 = mean(xf2,2);                 % mean fault2
mean_f3 = mean(xf3,2);                 % mean fault3
mean_f4 = mean(xf4,2);                 % mean fault4
mean_f5 = mean(xf5,2);                 % mean fault5


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
S5 = cov(xf4');
S6 = cov(xf5');
Sw = S1 + S2 + S3 + S4 + S5 ;


Mean = mean(x,2);
SB1 = lenN  * (mean_n  - Mean)*(mean_n  - Mean)' ;
SB2 = lenF1 * (mean_f1 - Mean)*(mean_f1 - Mean)' ;
SB3 = lenF2 * (mean_f2 - Mean)*(mean_f2 - Mean)' ;
SB4 = lenF3 * (mean_f3 - Mean)*(mean_f3 - Mean)' ;
SB5 = lenF4 * (mean_f4 - Mean)*(mean_f4 - Mean)' ;
SB6 = lenF5 * (mean_f5 - Mean)*(mean_f5 - Mean)' ;
SB = SB1 + SB2 + SB3 + SB4  + SB5  + SB6 ;

%inv_Sw = inv(Sw) ;
inv_Sw = inv(Sw'*Sw)*Sw' ;
invSw_SB = inv_Sw * SB ;
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


FAULT1 = y(1:lenF1 ,:) ;
FAULT2 = y(lenF1+1:lenF1+lenF2 ,:) ;
FAULT3 = y(lenF1+lenF2+1:lenF1+lenF2+lenF3 ,:) ;
FAULT4 = y(lenF1+lenF2+lenF3+1:lenF1+lenF2+lenF3+lenF4 ,:) ;
FAULT5 = y(lenF1+lenF2+lenF3+lenF4+1:lenF1+lenF2+lenF3+lenF4+lenF5 ,:) ;
NORMAL = y(lenF1+lenF2+lenF3+lenF4+lenF5+1:lenF1+lenF2+lenF3+lenF4+lenF5+lenN ,:) ; 


%% FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
FAULT1 = FAULT1*2 + 1 ;
FAULT2 = sqrt(FAULT2) +2 ;
FAULT3 = FAULT3.^2 - 4 ;
FAULT4 = FAULT4 .*2  + 2 ;
FAULT5 = FAULT5 + 2 ;


%% 3D plot
i=1 ; j=2; k=3;
figure; 
plot3(FAULT1(:,i),FAULT1(:,j),FAULT1(:,k),'bo','MarkerFaceColor','b'); hold on; grid on;
plot3(FAULT2(:,i),FAULT2(:,j),FAULT2(:,k),'ro','MarkerFaceColor','r');
plot3(FAULT3(:,i),FAULT3(:,j),FAULT3(:,k),'yo','MarkerFaceColor','y');
plot3(FAULT4(:,i),FAULT4(:,j),FAULT4(:,k),'ko','MarkerFaceColor','k');
plot3(FAULT5(:,i),FAULT5(:,j),FAULT5(:,k),'co','MarkerFaceColor','c');
plot3(NORMAL(:,i),NORMAL(:,j),NORMAL(:,k),'go');
title('LDA'); legend('FAULT1','FAULT2','FAULT3','FAULT4','FAULT5','NORMAL');


%% ------------------------------------ target --------------------------------------
y(1:lenF1 ,6) = 1 ;
y(lenF1+1:lenF1+lenF2 ,6) = 2 ;
y(lenF1+lenF2+1:lenF1+lenF2+lenF3 ,6) = 3 ;
y(lenF1+lenF2+lenF3+1:lenF1+lenF2+lenF3+lenF4 ,6) = 4 ;
y(lenF1+lenF2+lenF3+lenF4+1:lenF1+lenF2+lenF3+lenF4+lenF5 ,6) = 5 ;
y(lenF1+lenF2+lenF3+lenF4+lenF5+1:lenF1+lenF2+lenF3+lenF4+lenF5+lenN ,6) = 6 ; 
data = y ;

%%
FAULT1 = [FAULT1 ,ones(lenF1,1).*1] ;
FAULT2 = [FAULT2 ,ones(lenF2,1).*2] ;
FAULT3 = [FAULT3 ,ones(lenF3,1).*3] ;
FAULT4 = [FAULT4 ,ones(lenF4,1).*4] ;
FAULT5 = [FAULT5 ,ones(lenF5,1).*5] ;
NORMAL = [NORMAL ,ones(lenN ,1).*6] ;
data = [FAULT1;FAULT2;FAULT3;FAULT4;FAULT5;NORMAL];

FAULT1(:,end)=[] ;
FAULT2(:,end)=[] ;
FAULT3(:,end)=[] ;
FAULT4(:,end)=[] ;
FAULT5(:,end)=[] ;
NORMAL(:,end)=[] ;


%% ------------------------------ save data mat format -------------------------------------------
save('data_after_LDA.mat','data','FAULT1','FAULT2','FAULT3','FAULT4','FAULT5','NORMAL');


    
    
    
    