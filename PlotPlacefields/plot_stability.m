% S1
% T=3.599613200000000e+03;
% t_tr=  2.6829e+03;
% S2
% T=  2.8030e+03+973.5183;
% t_o= 973.5183;
% S3

 T= 924.3243+3.0969e+03;
t_tr=3.0969e+03;
%% plot pca and std
sample_win=10;
time_win=[2]
% filename1='DJ52_S3C';
filename1='DJ52_S3_1T';
load([filename1 '.mat'])
recording='c';
figure
% close all
for i=1

[mean_coeaff,std_coeaff]=plot_pca (filename1,int_celln(i),sample_win,time_win,recording);

end


%% make noise and displacement
sample_win=10;
time_win=[2]
filename1='DJ25_2';
load([filename1 '.mat'])
recording='a';
% [mean_coeaff,std_coeaff]=plot_pca (filename1,int_celln(13),sample_win,time_win,recording);

% load (filename1)
sesion=3;
for time_win=[2]
    
stability(filename1, sample_win,time_win,recording,sesion);


% save(['chr_noise_' num2str(time_win) '_' num2str(sesion) ],'noise')
% save(['chr_displacement_' num2str(time_win) '_' num2str(sesion)],'displacement')
% % 
% save(['A_noise_' num2str(time_win) '_' num2str(sesion) ],'noise')
% save(['A_displacement_' num2str(time_win) '_' num2str(sesion)], 'displacement')
end
%%

clear
clc
% celln=7;
% a=16;
% sesion=2;
sn1=31;
sn2=16;
mean_coeaff1=[];
mean_coeaff2=[];
std_coeaff1=[];
std_coeaff2=[];

for sesion=1:3
    cn=[16 7 10];
    for k=1:cn(sesion)
        
        if sesion==2
            i=2;
            load(['mean_coeaff_' num2str(k) '_' num2str(i)  '_'  num2str(sesion)],'mean_coeaff')
            load(['std_coeaff_' num2str(k) '_' num2str(i)  '_'  num2str(sesion)],'std_coeaff')
%             mean_coeaff=mean_coeaff/(max(mean_coeaff));
            mean_coeaff=flip(mean_coeaff);
            T1=1:length(mean_coeaff)/(sn1):length(mean_coeaff);
            
            mean_coeaff1=[mean_coeaff1; interp1(1:length(mean_coeaff),mean_coeaff,T1)];
            std_coeaff1=[std_coeaff1; interp1(1:length(std_coeaff),std_coeaff,T1)];
            
            i=1;
            load(['mean_coeaff_' num2str(k) '_' num2str(i)  '_'  num2str(sesion)],'mean_coeaff')
            load(['std_coeaff_' num2str(k) '_' num2str(i)  '_'  num2str(sesion)],'std_coeaff')
%             mean_coeaff=mean_coeaff/(max(mean_coeaff));
            mean_coeaff=flip(mean_coeaff);
            T2=1:length(mean_coeaff)/(sn2+mod(length(mean_coeaff),2)):length(mean_coeaff);
            
            mean_coeaff2=[mean_coeaff2; interp1(1:length(mean_coeaff),mean_coeaff,T2)];
            std_coeaff2=[std_coeaff2; interp1(1:length(std_coeaff),std_coeaff,T2)];
        else
            
            i=1;
            load(['mean_coeaff_' num2str(k) '_' num2str(i)  '_'  num2str(sesion)],'mean_coeaff')
            load(['std_coeaff_' num2str(k) '_' num2str(i)  '_'  num2str(sesion)],'std_coeaff')
%             mean_coeaff=mean_coeaff/(max(mean_coeaff));
            T1=1:length(mean_coeaff)/(sn1):length(mean_coeaff);
            
            mean_coeaff1=[mean_coeaff1; interp1(1:length(mean_coeaff),mean_coeaff,T1)];
           std_coeaff1=[std_coeaff1; interp1(1:length(std_coeaff),std_coeaff,T1)];
            
            i=2;
            load(['mean_coeaff_' num2str(k) '_' num2str(i)  '_'  num2str(sesion)],'mean_coeaff')
            load(['std_coeaff_' num2str(k) '_' num2str(i)  '_'  num2str(sesion)],'std_coeaff')
%             mean_coeaff=mean_coeaff/(max(mean_coeaff));
            T2=1:length(mean_coeaff)/(sn2+mod(length(mean_coeaff),2)):length(mean_coeaff);
            
            mean_coeaff2=[mean_coeaff2; interp1(1:length(mean_coeaff),mean_coeaff,T2)];
            std_coeaff2=[std_coeaff2; interp1(1:length(std_coeaff),std_coeaff,T2)];
        end
    end
end


%%
sn=20;
  T=a(1):length(a)/sn:a(end)';
  X=interp1(a,a,T);
%%
% clear
% clc
% celln=10
% chr_noise=[];
% chr_displacement=[];
% chr_noise_o=[];
% chr_displacement_o=[];
sesion=1;

% for sesion=1:3
  cn=[16 7 10];
for k=1%:cn(sesion)
   
i=1;

load(['mean_coeaff_' num2str(k) '_' num2str(i) '_'  num2str(sesion) ],'mean_coeaff')
load(['std_coeaff_' num2str(k) '_' num2str(i) '_'  num2str(sesion) ],'std_coeaff')
% mean_coeaff1=mean_coeaff;
% std_coeaff1=std_coeaff;
mean_coeaff1=mean_coeaff/(max(mean_coeaff));
std_coeaff1=std_coeaff/(max(mean_coeaff));

chr_noise=[chr_noise mean(std_coeaff1)];
chr_displacement=[chr_displacement mean(abs(diff( mean_coeaff1))) ];

i=2;

load(['mean_coeaff_' num2str(k) '_' num2str(i) '_'  num2str(sesion) ],'mean_coeaff')
load(['std_coeaff_' num2str(k) '_' num2str(i) '_'  num2str(sesion) ],'std_coeaff')
% mean_coeaff2=mean_coeaff;
% std_coeaff2=std_coeaff;
mean_coeaff2=mean_coeaff/(max(mean_coeaff));
std_coeaff2=std_coeaff/(max(mean_coeaff));

chr_noise_o=[chr_noise_o mean(std_coeaff2)];
chr_displacement_o=[chr_displacement_o mean(abs(diff(mean_coeaff2))) ];

% mean_coeaffT(k,:)=[mean_coeaff1 mean_coeaff2];
% std_coeaffT(k,:)=[std_coeaff1 std_coeaff2];


end

%%
save('chr_noise_o' ,'chr_noise_o')
save('chr_displacement_o'  ,'chr_displacement_o')
save('chr_noise' ,'chr_noise')
save('chr_displacement' ,'chr_displacement')

%% %acute


clear
clc

A_noise=[];
A_displacement=[];


for sesion=1:3
    cn=[35 33 14];
    
    for k=1:cn(sesion)
        
        
        
        load(['mean_coeaff_' num2str(k) '_'   num2str(sesion) ],'mean_coeaff')
        load(['std_coeaff_' num2str(k) '_'   num2str(sesion) ],'std_coeaff')
       
        mean_coeaff1=mean_coeaff/(max(mean_coeaff));
        std_coeaff1=std_coeaff/(max(mean_coeaff));
        
        A_noise=[A_noise mean(std_coeaff1)];
        A_displacement=[A_displacement mean(abs(diff( mean_coeaff1))) ];
        
    end
end


%%

save('A_noise' ,'A_noise')
save('A_displacement' ,'A_displacement')


%%

load(['mean_coeaff1'])
load(['std_coeaff1'])
load(['mean_coeaff2'])
load(['std_coeaff2']) 
mean_coeaffT=[mean_coeaff1 mean_coeaff2];
std_coeaffT=[std_coeaff1 std_coeaff2];

%%
%s1
% tr=22;
% tr=44;
% S2
% tr=8;
% S3
tr=30
for k=1:33
% plot(DX(k,:),1:57,'-.')
% errorbar(1:length(mean_coeaffT),mean_coeaffT(k,:),std_coeaffT(k,:),'--')
plot(abs(diff(mean_coeaffT(k,:))),'.','Markersize',1)
hold on 

end
hold on
%%
errorbar(1:length(mean_coeaffT)-1,mean(abs(diff(mean_coeaffT,1,2))),std(abs(diff(mean_coeaffT,1,2))),'k.')
hold on
% plot(mean(abs(diff(mean_coeaffT,1,2))),'k','Linewidth',2)
hold on
%  plot([tr tr],[-0.05 0.25],'k','Linewidth',2);
plot(1:length(mean_coeaffT)-1,mean(abs(diff(mean_coeaffT,1,2))),'k.','markersize',15)

xlim([24 36])
ylim([-0.05 0.25])

set(gca,'CameraUpVector',[1,0,0],'YDir','reverse','XAxisLocation','top')
%  set(gca,'box','off', 'XTick',[],'YTick',[],'TickDir','out')
%%
% tr=23;
%  mean_coeaff=mean_coeaff/(max(mean_coeaff));

figure
    plot(abs(diff(mean_coeaff)),1:21,'k.','Markersize',15)
%     hold on
%     plot([-0.05 0.1],[tr tr],'r','Linewidth',1);
 xlim([-0.05 0.05])
%  set(gca,'box','off', 'XTick',[],'YTick',[],'TickDir','out')
%  ylim([0 30])


%%

hold on
plot([tr tr],[0 1],'r');
set(gca,'CameraUpVector',[1,0,0],'YDir','reverse','XAxisLocation','top')

%% loading files
close all

clear 

time_win =2;
chr_noise=[];
chr_noise_o=[];
chr_displacement=[];
chr_displacement_o=[];
A_noise=[];
A_displacement=[];

 for sesion=1:3
i=sesion;
filename=(['chr_noise_' num2str(time_win) '_' num2str(sesion)]);
load([filename '.mat'])

chr_noise=[chr_noise noise(1,:)];
chr_noise_o=[chr_noise_o noise(2,:)];

filename=(['chr_displacement_' num2str(time_win) '_' num2str(sesion)]);
load([filename '.mat'])

chr_displacement=[chr_displacement displacement(1,:)];
chr_displacement_o=[chr_displacement_o displacement(2,:)];

filename=(['A_noise_' num2str(time_win) '_' num2str(sesion)]);
load([filename '.mat'])
A_noise=[A_noise noise];

filename=(['A_displacement_' num2str(time_win) '_' num2str(sesion)]);
load([filename '.mat'])
A_displacement=[A_displacement displacement]; 
 end
 x=A_displacement;
 y=chr_displacement;
 z=chr_displacement_o;

%% plotting noise and displace ment
cd N2min_A
load('A_displacement.mat')
load('A_noise.mat')
cd ..
cd N2min_chr
load('chr_displacement.mat')
load('chr_displacement_o.mat')
load('chr_noise.mat')
load('chr_noise_o.mat')
cd ..
A_noise(79)=[];
A_displacement(79)=[];
%% scatter  plot
figure

scatter(A_displacement,A_noise,'LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b')

hold on

scatter(chr_displacement,chr_noise,'LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r')

hold on

scatter(chr_displacement_o,chr_noise_o,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k')



h=lsline;
set(h(1),'color','k')
set(h(2),'color','r')
set(h(3),'color','b')

%%  bar graph
  figure    
subplot(1,2,1)
bar(1,mean(A_noise),'FaceColor',[0 0 .7])
hold on
b1=errorbar(1,mean(A_noise),std(A_noise)/sqrt(length(A_noise)))
set(b1,'linewidth',3,'color','k')
hold on
bar(2,mean(chr_noise),'FaceColor',[ 1 0.7 0.4 ])
hold on
b2=errorbar(2,mean(chr_noise),std(chr_noise)/sqrt(length(chr_noise)))
set(b2,'linewidth',3,'color','k')
hold on
bar(3,mean(chr_noise_o),'FaceColor',[ .9 0 0 ])
hold on
b3=errorbar(3,mean(chr_noise_o),std(chr_noise_o)/sqrt(length(chr_noise_o)))
set(b3,'linewidth',3,'color','k')
set(gca,'box','off','TickDir','out','XTick',[])
subplot(1,2,2)

bar(1,mean(A_displacement),'FaceColor',[0 0 .7])
hold on
b4=errorbar(1,mean(A_displacement),std(A_displacement)/sqrt(length(A_displacement)))
set(b4,'linewidth',3,'color','k')
hold on
bar(2,mean(chr_displacement),'FaceColor',[ 1 0.7 0.4 ])
hold on
b5=errorbar(2,mean(chr_displacement),std(chr_displacement)/sqrt(length(chr_displacement)))
set(b5,'linewidth',3,'color','k')
hold on
bar(3,mean(chr_displacement_o),'FaceColor',[ .9 0 0 ])
hold on
b6=errorbar(3,mean(chr_displacement_o),std(chr_displacement_o)/sqrt(length(chr_displacement_o)))
set(b6,'linewidth',3,'color','k')
% legend('Location', 'northeastoutside')
set(gca,'box','off','TickDir','out','XTick',[])


 %% statistics
name1=repmat({'gp1'},1,length(x));
name2=repmat({'gp2'},1,length(y));
name3=repmat({'gp3'},1,length(z));

figure
groups=[y,z];
names=[name2,name3];
[p,anovatab,stats] = anova1(groups, names);

%%
w=0:0.03:0.5;
figure
n1=hist(A_noise,w);
b1=bar(n1/length(A_noise));
b1.FaceAlpha=0.5;
figure
n2=hist(chr_noise_o,w);
b2=bar(n2/length(chr_noise_o));
b2.FaceAlpha=0.5;
figure
n3=hist(chr_noise,w);
b3=bar(n3/length(chr_noise));
b3.FaceAlpha=0.5;


%%
figure
subplot(3,1,1)
b1=bar(n1/length(A_noise),'BarWidth', 1,'FaceColor',[0 0 .7],'EdgeColor',[0 0 .7]);
subplot(3,1,3)
b2=bar(n2/length(chr_noise_o),'BarWidth', 1,'FaceColor',[ .9 0 0 ],'EdgeColor',[ .9 0 0 ]);
subplot(3,1,2)
b3=bar(n3/length(chr_noise),'BarWidth', 1,'FaceColor',[ 1 0.7 0.4 ],'EdgeColor',[ 1 0.7 0.4 ]);
% b1.FaceAlpha=0.2;
% b2.FaceAlpha=0;
% b3.FaceAlpha=0.;






