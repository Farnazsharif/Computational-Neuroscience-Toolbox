function  Boxplot_f(CA3,CA1,x_lable)
%%

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


%%
C={'k','r'};
% close all
x_min=0;
x_max=1;
if min(x)<0 
    x_min=min(x)-0.2;
elseif max(x)>1
    x_max=max(x)+2;
end



% %  TRD2 & TRD2 & TRD2
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

Pa=[];Pk=[];
[Pa(1,1),Pk(1,1)] = Stest( V2a',V2' );[Pa(1,2),Pk(1,2)] = Stest( V5a',V5');
[Pa(2,1),Pk(2,1)] = Stest( V3a',V3' );[Pa(2,2),Pk(2,2)] = Stest( V6a',V6');
Pa_2=Pa;Pk_2=Pk;

Pa=[];Pk=[];
[Pa(1,1),Pk(1,1)] = Stest( V2a',V5a' );[Pa(1,2),Pk(1,2)] = Stest( V2',V5');
[Pa(2,1),Pk(2,1)] = Stest( V3a',V6a' );[Pa(2,2),Pk(2,2)] = Stest( V3',V6');
Pa_3=Pa;Pk_3=Pk;

Pa=[];Pk=[];
[Pa(1,1),Pk(1,1)] = Stest( V2a',V3a' );[Pa(1,2),Pk(1,2)] = Stest( V5a',V6a');
[Pa(2,1),Pk(2,1)] = Stest( V2',V3'  );[Pa(2,2),Pk(2,2)] = Stest( V5',V6');
Pa_4=Pa;Pk_4=Pk;



%% box plot
C=[8/256 133/256 161/256;193/256 90/256 99/256;99/256 209/256 62/256;88/256 88/256 88/256;99/256 209/256 62/256;88/256 88/256 88/256];
ed_C={'b';'r';'b';'b';'r';'r'};


figure('Position',[10 -500 500 800]);

hp=[];
data=[];



ORD=[6 5 4 3 2 1];
Xpos =[1 1.25  2.25 2.5 2.75 3;1.5 1.75  3.5 3.75 4 4.25];
data_1=[V1a;V4a;V2a;V3a;V5a;V6a];
data_2=[V1;V4;V2;V3;V5;V6];

for kk=1:2
group_1= [repmat(Xpos(kk,1),length(V1a),1);repmat(Xpos(kk,2),length(V4a),1);repmat(Xpos(kk,3),length(V2a),1);repmat(Xpos(kk,4),length(V3a),1);repmat(Xpos(kk,5),length(V5a),1);repmat(Xpos(kk,6),length(V6a),1)];%; 
group_2= [repmat(Xpos(kk,1),length(V1),1);repmat(Xpos(kk,2),length(V4),1);repmat(Xpos(kk,3),length(V2),1);repmat(Xpos(kk,4),length(V3),1);repmat(Xpos(kk,5),length(V5),1);repmat(Xpos(kk,6),length(V6),1)];%; 

eval(['data=' 'data_' num2str(kk) ';']);
eval(['group=' 'group_' num2str(kk) ';']);
hpb=boxplot(data,group,'positions', Xpos(kk,:),'Notch','on');%,'Labels',{'all','Cue','PI','all','Cue','PI'}

h = findobj(gca,'Tag','Box');
hold on


for u=6:-1:1

    hp=patch(get(h(u),'XData'),get(h(u),'YData'),'y','FaceAlpha',0.9,'FaceColor',C(ORD(u),:),'EdgeColor',ed_C{ORD(u)},'LineWidth',1.5);
    if kk==2
    hatchfill2(hp,'single','HatchAngle',30,'HatchSpacing',10,'HatchLineStyle','-.');
    end

end 
end
xlim([0.75 4.5])
set(gca,'XTickLabel',[])
set(gca,'xtick',[])
set(gca,'FontSize',20)
title(x_lable,'FontSize',12)

% % M_data=[mean(V1a) mean(V4a) mean(V2a) mean(V3a) mean(V5a) mean(V6a) mean(V1) mean(V4) mean(V2) mean(V3) mean(V5) mean(V6)];
% % plot([Xpos(1,:) Xpos(2,:)],M_data, 'dk','LineWidth',1.5)


%% adding sigbars
groups_t={[Xpos(1,1) Xpos(1,2)],[Xpos(2,1) Xpos(2,2)],[Xpos(1,1) Xpos(2,1)],[Xpos(1,2) Xpos(2,2)]};
test=[Pk_1(2,1),Pk_1(2,2),Pk_1(1,1),Pk_1(1,2)];
ndx=(find(test<THR));
hp=sigstar(groups_t(ndx),test(ndx),0,THR);

% 2
groups_t={[Xpos(1,3) Xpos(1,4)],[Xpos(1,5) Xpos(1,6)],[Xpos(2,3) Xpos(2,4)],[Xpos(2,5) Xpos(2,6)]};
test=[Pk_4(1,1),Pk_4(1,2),Pk_4(2,1),Pk_4(2,2)];
ndx=(find(test<THR));
hp=sigstar(groups_t(ndx),test(ndx),0,THR);
% 3
groups_t={[Xpos(1,3) Xpos(1,5)],[Xpos(2,3) Xpos(2,5)],[Xpos(1,4) Xpos(1,6)],[Xpos(2,4) Xpos(2,6)]};
test=[Pk_3(1,1),Pk_3(1,2),Pk_3(2,1),Pk_3(2,2)];
ndx=(find(test<THR));
hp=sigstar(groups_t(ndx),test(ndx),0,THR);
% 4
groups_t={[Xpos(1,3) Xpos(2,3)],[Xpos(1,5) Xpos(2,5)],[Xpos(1,4) Xpos(2,4)],[Xpos(1,6) Xpos(2,6)]};
test=[Pk_2(1,1),Pk_2(1,2),Pk_2(2,1),Pk_2(2,2)];
ndx=(find(test<THR));
hp=sigstar(groups_t(ndx),test(ndx),0,THR);

% ylim([-200 200])


