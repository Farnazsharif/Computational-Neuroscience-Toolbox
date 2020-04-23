% filename1='PowerT3_n_15';
% load([filename1 '.mat']);
% close all
figure
plot((j1))
hold on
y1=200000;


xtr1=find(trial(:,4)==1, 1, 'last')
%   s2=find(I>xtr1, 1, 'first')-1
 s2=find(pow(:,2)>xtr1, 1, 'first')-1
xi=s2*ones(y1,1);


y=1:1:y1;
plot(xi,y,'--k','LineWidth',2)
hold on
xtr2=find(trial(:,4)==2, 1,'last');
% s2=find(I>xtr2, 1, 'first')-1
s2=find(pow(:,2)>xtr2, 1, 'first')-1
xj=s2*ones(y1,1);
plot(xj,y,'--k','LineWidth',2)

%% Making corellaton bars
clf
L=[];
% W=[ 1 6;7 41;42 48; 49 90;91 97;98 124;125 140;1 140 ]; % for rest
%  W=[ 1 6;7 41;42 48; 49 90;91 97;98 124;125 137;1 137 ];
W=[ 1 6;7 33;34 41;42 48; 49 82;83 90;91 97;98 117;118 124;125 137;1 137];
L=lc(pow(:,2));

for i=1:length(W)
    
L1=L(W(i,1):W(i,2));
Po1=j1(W(i,1):W(i,2));
[R1,P1] = corrcoef(L1,Po1);
RF(i)=R1(1,2);
PF(i)=P1(1,2);
end 

b=diag(RF);
hold on
for k=1:length(W)
    h(k)=bar(1:length(RF),b(:,k));
end


% Setnam={ 'B set1'; 'set1'; 'B set2';'set2'; ; 'B set3';'set3';'End'; 'all'}
% Setnam={ 'B set1'; 'set1'; 'E set1'; 'B set2';'set2'; 'E set2'; 'B set3';'set3'; 'E set3';'End'; 'all'};
% set(gca,'xticklabel',Setnam)
ylabel('Corelation of licking and thtea power at rest','fontsize',14)
% set(h(1),'facecolor','r') ;set(h(2),'facecolor','r')
% set(h(3),'facecolor','b') ;set(h(4),'facecolor','b') 
% set(h(5),'facecolor','k') ;set(h(6),'facecolor','k') 
% set(h(7),'facecolor','cyan') ;set(h(8),'facecolor','magenta') 
% 
set(h(1),'facecolor','r') ;set(h(2),'facecolor','r') ;set(h(3),'facecolor','r') 
set(h(4),'facecolor','b') ;set(h(5),'facecolor','b') ;set(h(6),'facecolor','b') 
set(h(7),'facecolor','k') ;set(h(8),'facecolor','k') ;set(h(9),'facecolor','k')
set(h(10),'facecolor','cyan') ;set(h(11),'facecolor','magenta') 
set(gcf,'color','w');
 xlim([0.5 11.5])


%% Makin lick bar
L=[];
W=[ 1 6;7 41;42 48; 49 90;91 97;98 124];
% W=[ 1 6;7 33;34 41;42 48; 49 82;83 90;91 97;98 117;118 124];
L=lc(pow(:,2));
for i = 1:length(W)
    Lic(i)=mean(L(W(i,1):W(i,2)),1);
end

 figure
% bar(Lic)
b=diag(Lic);
hold on
for k=1:length(W)
    h(k)=bar(1:length(Lic),b(:,k));
end
% Setnam={ 'B set1'; 'set1'; 'B set2';'set2'; ; 'B set3';'set3';};
% Setnam={ 'B set1'; 'set1'; 'E set1'; 'B set2';'set2'; 'E set2'; 'B set3';'set3'; 'E set3'};
% set(gca,'xticklabel',Setnam)

set(h(1),'facecolor','r') ;set(h(2),'facecolor','r')
set(h(3),'facecolor','b') ;set(h(4),'facecolor','b') 
set(h(5),'facecolor','k') ;set(h(6),'facecolor','k') 

% set(h(1),'facecolor','r') ;set(h(2),'facecolor','r') ;set(h(3),'facecolor','r') 
% set(h(4),'facecolor','b') ;set(h(5),'facecolor','b') ;set(h(6),'facecolor','b') 
% set(h(7),'facecolor','k') ;set(h(8),'facecolor','k') ;set(h(9),'facecolor','k') 
 xlim([0.5 6.5])
 ylabel('Mean Theta Power','fontsize',14)
 set(gcf,'color','w');
 
%% making window for plot J matrix power
clc
W=[ 1 6;7 41;42 48; 49 90;91 97;98 124];
% W=[ 1 6;7 33;34 41;42 48; 49 82;83 90;91 97;98 117;118 124];
J=[];
for i = 1:length(W)
    J(i)=mean(j1(W(i,1):W(i,2)),1);
end

for i=1:length(W)
    Width(i)= length(W(i,1):W(i,2));
end

% Width=Width*0.04;
Width=Width*0.03;
b=diag(J);
figure
hold on
for k=1:length(W)
    h(k)=bar(1:length(J),b(:,k),Width(k),'FaceColor', [k/20 (.05*k)/2 (k+50)/100]);
end
Setnam={ 'B set1'; 'set1'; 'B set2';'set2'; ; 'B set3';'set3';};
% Setnam={ 'B set1'; 'set1'; 'E set1'; 'B set2';'set2'; 'E set2'; 'B set3';'set3'; 'E set3'};
set(gca,'xticklabel',Setnam)

set(h(1),'facecolor','r') ;set(h(2),'facecolor','r')
set(h(3),'facecolor','b') ;set(h(4),'facecolor','b') 
set(h(5),'facecolor','k') ;set(h(6),'facecolor','k') 
% set(h(1),'facecolor','r') ;set(h(2),'facecolor','r') ;set(h(3),'facecolor','r') 
% set(h(4),'facecolor','b') ;set(h(5),'facecolor','b') ;set(h(6),'facecolor','b') 
% set(h(7),'facecolor','k') ;set(h(8),'facecolor','k') ;set(h(9),'facecolor','k') 
 xlim([0.5 6.5])
 ylabel('Mean Theta Power','fontsize',14)
 set(gcf,'color','w');
%% for power matrix before making mean
W=[1 6;7 45;46 52;53 94;95 101;102 128;129 I(end)];
ndx=[];
for i=1: length(W)
   a= find(unique(I)>=W(i,1), 1, 'first');
   b= find(unique(I)>W(i,2), 1, 'first')-1;
   if isempty(b)==1
      b=length(unique(I));
   end
   ndx(i,1)=a;
   ndx(i,2)=b;
end
%%
for i=1:length(ndx);
PO(i)= mean(PowerT3(ndx(i,1):ndx(i,2),theta),1)
end

figure(2)
bar(PO)
Setnam={'Begin of set1'; 'rest of set1'; 'Begin of set2'; 'rest of set2';'Begin of set3'; 'rest of set3' };
set(gca,'xticklabel',Setnam)
ylabel('Mean Theta Power','fontsize',14)
set(gcf,'color','w');
%% makin windo for corelation with j

W=[1 6;7 41;42 46;47 90;91 94;95 124;125 137];
% (initial window = 7:45);
% (initial window = 46:50);
%(initial window = 51:94);
%(initial window = 95:98);
% (initial window = 99:128);
% (initial window = 129:142);
