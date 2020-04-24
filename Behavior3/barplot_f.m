function  barplot_f(CA3,CA1,x_lable)
THR=0.05;

% CA3
V2=CA3.C_1;
V3=CA3.P_1;
V5=CA3.C_2;
V6=CA3.P_2;


V1=[V2; V3];
V4=[V5; V6];
x=V1;
% CA1
V2a=CA1.C_1;
V3a=CA1.P_1;
V5a=CA1.C_2;
V6a=CA1.P_2;

V1a=[V2a; V3a];
V4a=[V5a; V6a];

C={'k','r'};
% close all
x_min=0;
x_max=1;
if min(x)<0 
    x_min=min(x)-0.2;
elseif max(x)>1
    x_max=max(x)+2;
end

fx=900;
fh=300;

figure1=figure('Position',[100 100 fx fh]);

%  TRD2 & TRD2 & TRD2
g_an1=0;
g_an2=1;
g_an3=.16;
g_an4=0.07;

g_1=0.04; 
g_2= 0.12;
g_3= 0.13;
g_4=0.55;
Pa=[];Pk=[];
[Pa(1,1),Pk(1,1)] = Stest( V1a',V1' );[Pa(1,2),Pk(1,2)] = Stest( V4a',V4');
[Pa(2,1),Pk(2,1)] = Stest( V1a',V4a' );[Pa(2,2),Pk(2,2)] = Stest( V1',V4');
Pa_1=Pa;Pk_1=Pk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Annotation


% dim1 = [g_an1 g_an2-(0.08)*1 g_an3 g_an4];
% dim2 = [g_an1 g_an2-(0.08)*2 g_an3 g_an4];
% dim3 = [g_an1 g_an2-(0.08)*3 g_an3 g_an4];
% dim4 = [g_an1 g_an2-(0.08)*4 g_an3 g_an4];
% 
% str1=(['CA1-TRD1, CA3-TRD1  ' 'Pa=' num2str(round(Pa(1,1),4)) '  '  'Pk=' num2str(round(Pk(1,1),4))]);
% str2=(['CA1-TRD2, CA3-TRD2  ' 'Pa=' num2str(round(Pa(1,2),4)) '  '  'Pk=' num2str(round(Pk(1,2),4))]);
% str3=(['CA1-TRD1, CA1-TRD2  ' 'Pa=' num2str(round(Pa(2,1),4)) '  '  'Pk=' num2str(round(Pk(2,1),4))]);
% str4=(['CA3-TRD1, CA3-TRD2  ' 'Pa=' num2str(round(Pa(2,2),4)) '  '  'Pk=' num2str(round(Pk(2,2),4))]);
% 
% if Pa(1,1) < THR || Pk(1,1) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim1,'String',str1,'FitBoxToText','on','LineStyle','none','color',C{i});
% 
% if Pa(1,2) < THR || Pk(1,2) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim2,'String',str2,'FitBoxToText','on','LineStyle','none','color',C{i});
% 
% if Pa(2,1) < THR || Pk(2,1) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim3,'String',str3,'FitBoxToText','on','LineStyle','none','color',C{i});
% 
% if Pa(2,2) < THR || Pk(2,2) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim4,'String',str4,'FitBoxToText','on','LineStyle','none','color',C{i});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axes('Position',[g_1 g_2 g_3 g_4])
% plot(sort(V4a),(1:length(V4a))/length(V4a),'Linewidth',3,'color','r')
% hold on
% plot(sort(V4),(1:length(V4))/length(V4), 'Linewidth',3,'color','r','LineStyle',':')
% hold on
% plot(sort(V1a),(1:length(V1a))/length(V1a),'Linewidth',3,'color','b')
% hold on
% plot(sort(V1),(1:length(V1))/length(V1), 'Linewidth',3,'color','b','LineStyle',':')
% box ('off')
%  
% xlim([x_min x_max])
% xlabel(x_lable,'FontWeight','bold')

% Cue cells

g_an1=g_an3*1.6;
g_1=g_an1; 
Pa=[];Pk=[];
[Pa(1,1),Pk(1,1)] = Stest( V2a',V2' );[Pa(1,2),Pk(1,2)] = Stest( V5a',V5');
[Pa(2,1),Pk(2,1)] = Stest( V3a',V3' );[Pa(2,2),Pk(2,2)] = Stest( V6a',V6');
Pa_2=Pa;Pk_2=Pk;

dim1 = [g_an1 g_an2-(0.08)*1 g_an3 g_an4];
dim2 = [g_an1 g_an2-(0.08)*2 g_an3 g_an4];
dim3 = [g_an1 g_an2-(0.08)*3 g_an3 g_an4];
dim4 = [g_an1 g_an2-(0.08)*4 g_an3 g_an4];

% str1=(['CA1-Cue1, CA3-Cue1  ' 'Pa=' num2str(round(Pa(1,1),4)) '  '  'Pk=' num2str(round(Pk(1,1),4))]);
% str2=(['CA1-Cue2, CA3-Cue2  ' 'Pa=' num2str(round(Pa(1,2),4)) '  '  'Pk=' num2str(round(Pk(1,2),4))]);
% str3=(['CA1-PI1, CA3-PI1  ' 'Pa=' num2str(round(Pa(2,1),4)) '  '  'Pk=' num2str(round(Pk(2,1),4))]);
% str4=(['CA1-PI2, CA3-PI2  ' 'Pa=' num2str(round(Pa(2,2),4)) '  '  'Pk=' num2str(round(Pk(2,2),4))]);
% 
% if Pa(1,1) < THR || Pk(1,1) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim1,'String',str1,'FitBoxToText','on','LineStyle','none','color',C{i});
% 
% if Pa(1,2) < THR || Pk(1,2) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim2,'String',str2,'FitBoxToText','on','LineStyle','none','color',C{i});
% 
% if Pa(2,1) < THR || Pk(2,1) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim3,'String',str3,'FitBoxToText','on','LineStyle','none','color',C{i});
% 
% if Pa(2,2) < THR || Pk(2,2) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim4,'String',str4,'FitBoxToText','on','LineStyle','none','color',C{i});
% axes('Position',[g_1 g_2 g_3 g_4])
% PI cells

g_an1=g_an3*3.2;
g_1=g_an1; 

Pa=[];Pk=[];
[Pa(1,1),Pk(1,1)] = Stest( V2a',V5a' );[Pa(1,2),Pk(1,2)] = Stest( V2',V5');
[Pa(2,1),Pk(2,1)] = Stest( V3a',V6a' );[Pa(2,2),Pk(2,2)] = Stest( V3',V6');
Pa_3=Pa;Pk_3=Pk;

dim1 = [g_an1 g_an2-(0.08)*1 g_an3 g_an4];
dim2 = [g_an1 g_an2-(0.08)*2 g_an3 g_an4];
dim3 = [g_an1 g_an2-(0.08)*3 g_an3 g_an4];
dim4 = [g_an1 g_an2-(0.08)*4 g_an3 g_an4];

% str1=(['CA1-Cue1, CA1-Cue2  ' 'Pa=' num2str(round(Pa(1,1),4)) '  '  'Pk=' num2str(round(Pk(1,1),4))]);
% str2=(['CA3-Cue1, CA3-Cue2  ' 'Pa=' num2str(round(Pa(1,2),4)) '  '  'Pk=' num2str(round(Pk(1,2),4))]);
% str3=(['CA1-PI1, CA1-PI2  ' 'Pa=' num2str(round(Pa(2,1),4)) '  '  'Pk=' num2str(round(Pk(2,1),4))]);
% str4=(['CA3-PI1, CA3-PI2  ' 'Pa=' num2str(round(Pa(2,2),4)) '  '  'Pk=' num2str(round(Pk(2,2),4))]);
% 
% if Pa(1,1) < THR || Pk(1,1) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim1,'String',str1,'FitBoxToText','on','LineStyle','none','color',C{i});
% 
% if Pa(1,2) < THR || Pk(1,2) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim2,'String',str2,'FitBoxToText','on','LineStyle','none','color',C{i});
% 
% if Pa(2,1) < THR || Pk(2,1) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim3,'String',str3,'FitBoxToText','on','LineStyle','none','color',C{i});
% 
% if Pa(2,2) < THR || Pk(2,2) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim4,'String',str4,'FitBoxToText','on','LineStyle','none','color',C{i});


g_an1=g_an3*4.8;


Pa=[];Pk=[];
[Pa(1,1),Pk(1,1)] = Stest( V2a',V3a' );[Pa(1,2),Pk(1,2)] = Stest( V5a',V6a');
[Pa(2,1),Pk(2,1)] = Stest( V2',V3'  );[Pa(2,2),Pk(2,2)] = Stest( V5',V6');
Pa_4=Pa;Pk_4=Pk;

dim1 = [g_an1 g_an2-(0.08)*1 g_an3 g_an4];
dim2 = [g_an1 g_an2-(0.08)*2 g_an3 g_an4];
dim3 = [g_an1 g_an2-(0.08)*3 g_an3 g_an4];
dim4 = [g_an1 g_an2-(0.08)*4 g_an3 g_an4];

% str1=(['CA1-Cue1, CA1-PI1  ' 'Pa=' num2str(round(Pa(1,1),4)) '  '  'Pk=' num2str(round(Pk(1,1),4))]);
% str2=(['CA1-Cue2, CA1-PI2  ' 'Pa=' num2str(round(Pa(1,2),4)) '  '  'Pk=' num2str(round(Pk(1,2),4))]);
% str3=(['CA3-Cue1, CA3-PI1  ' 'Pa=' num2str(round(Pa(2,1),4)) '  '  'Pk=' num2str(round(Pk(2,1),4))]);
% str4=(['CA3-Cue2, CA3-PI2  ' 'Pa=' num2str(round(Pa(2,2),4)) '  '  'Pk=' num2str(round(Pk(2,2),4))]);
% 
% if Pa(1,1) < THR || Pk(1,1) < THR
%     i=2;
% else
%     i=1;
% end
% %
% annotation('textbox',dim1,'String',str1,'FitBoxToText','on','LineStyle','none','color',C{i});
% %
% if Pa(1,2) < THR || Pk(1,2) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim2,'String',str2,'FitBoxToText','on','LineStyle','none','color',C{i});
% 
% if Pa(2,1) < THR || Pk(2,1) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim3,'String',str3,'FitBoxToText','on','LineStyle','none','color',C{i});
% 
% if Pa(2,2) < THR || Pk(2,2) < THR
%     i=2;
% else
%     i=1;
% end
% annotation('textbox',dim4,'String',str4,'FitBoxToText','on','LineStyle','none','color',C{i});



%% bar graph

% {bule,red,green,black,green,black}; 
C=[8/256 133/256 161/256;193/256 90/256 99/256;99/256 209/256 62/256;88/256 88/256 88/256;99/256 209/256 62/256;88/256 88/256 88/256];
ed_C={'b';'r';'b';'b';'r';'r'};
Vec={'V1a','V4a','V2a','V3a','V5a','V6a';'V1','V4','V2','V3','V5','V6'};
sp=1.25;
Xpos_b=[1 1+sp 7 7+sp 7+sp*2 7+sp*3;1+sp*2 1+sp*3 13 13+sp 13+sp*2 13+sp*3]; 
% Xpos_b=[1 1+sp 7 13 7+sp*2 13+sp*2;1+sp*2 1+sp*3 7+sp  13+sp 7+sp*3 13+sp*3]; % cue cue PI PI

axes('Position',[0.07 0.07 0.43 .55])

T=[];


for k=1:2
    
for l=1:6
    
eval(['T=' Vec{k,l} ';']);
%**************************************************************************
% [~,meanO]=meanphase(T); % circular mean
% bar_value=meanO;
bar_value=nanmean(T); % arethmatic mean 
%**************************************************************************


hp=bar(Xpos_b(k,l),bar_value,'FaceColor',C(l,:),'EdgeColor',ed_C{l},'LineWidth',2);

set(gca ,'XTickLabel',[])
hold on
b1=errorbar(Xpos_b(k,l),bar_value,nanstd(T)/sqrt(length(T)));
set(b1,'linewidth',2,'color','k')

% if k==2
%   hatchfill2(hp,'single','HatchAngle',30,'HatchSpacing',20,'HatchLineStyle','-.');  
% end

end

end

set(gca,'FontSize',20)
set(gca,'xtick',[])
set(gca,'xlim',([0.25 17.5]))
b1 = blanks(10);b2=blanks(25);b3=blanks(35);b4=blanks(10);
xlabel(['CA1' b1 'CA3' b2 'CA1' b3 'CA3' b4],'FontWeight','bold','FontSize',12)
title(x_lable,'FontSize',12)

% ylim([0 .9])
%% adding sigbars
%1

groups_t={[Xpos_b(1,1) Xpos_b(1,2)],[Xpos_b(2,1) Xpos_b(2,2)],[Xpos_b(1,1) Xpos_b(2,1)],[Xpos_b(1,2) Xpos_b(2,2)]};
test=[Pa_1(2,1),Pa_1(2,2),Pa_1(1,1),Pa_1(1,2)];
ndx=(find(test<THR));
b1=sigstar(groups_t(ndx),test(ndx),0,THR);

% 2
groups_t={[Xpos_b(1,3) Xpos_b(1,4)],[Xpos_b(1,5) Xpos_b(1,6)],[Xpos_b(2,3) Xpos_b(2,4)],[Xpos_b(2,5) Xpos_b(2,6)]};
test=[Pa_4(1,1),Pa_4(1,2),Pa_4(2,1),Pa_4(2,2)];
ndx=(find(test<THR));
b1=sigstar(groups_t(ndx),test(ndx),0,THR);

% 3
groups_t={[Xpos_b(1,3) Xpos_b(1,5)],[Xpos_b(2,3) Xpos_b(2,5)],[Xpos_b(1,4) Xpos_b(1,6)],[Xpos_b(2,4) Xpos_b(2,6)]};
test=[Pa_3(1,1),Pa_3(1,2),Pa_3(2,1),Pa_3(2,2)];
ndx=(find(test<THR));
b1=sigstar(groups_t(ndx),test(ndx),0,THR);

% 4
groups_t={[Xpos_b(1,3) Xpos_b(2,3)],[Xpos_b(1,5) Xpos_b(2,5)],[Xpos_b(1,4) Xpos_b(2,4)],[Xpos_b(1,6) Xpos_b(2,6)]};
test=[Pa_2(1,1),Pa_2(1,2),Pa_2(2,1),Pa_2(2,2)];
ndx=(find(test<THR));
b1=sigstar(groups_t(ndx),test(ndx),0,THR);
end

