function ReorderPower(filename,shank2,Frq) 
% ReorderPower('po_',8,8) 
% shank2=8; sh 1:8 = 0;  sh 8:16 = 8;
% filename='po_';
% Frq=9;
% chorder=[5 4 6 3 7 2 8 1];
chorder=[1 8 2 7 3 6 4 5];
shank_sapce=100;
nsite=8;
nshank=8;
reorder=[];
for ii = 1:16
    reorder = [reorder chorder+(ii-1)*nsite];
end

ch=reorder;

load ([filename,num2str(65)]);

P=PP(:,1);
matC=[];
smoothC=smooth1D(P,10,1);
xbinNumber=100;
matC2=reshape(smoothC,xbinNumber,length(P)/xbinNumber)';
[a,b]=size(matC2);
%%

Wpo=zeros(a*8,b*8*2+7*shank_sapce);


y=1:a:(a*nsite)+1;
shankx=[0:(b*2+shank_sapce):(b*2+shank_sapce)*7];
% CFCMAT=[];
for k=1:8
    
x=[1+shankx(k) (b+1)+shankx(k) (1)+shankx(k) (b+1)+shankx(k) (1)+shankx(k) (b+1)+shankx(k) (1)+shankx(k) (b/2)+shankx(k)];
shank=ch(((k+shank2)*8)-7:(k+shank2)*8);

for i= 1:8
    
load ([filename,num2str(shank(i))]);

P=PP(:,Frq);
matC=[];
smoothC=smooth1D(P,10,1);
xbinNumber=100;
matC2=reshape(smoothC,xbinNumber,length(P)/xbinNumber)';

% Wpo(y(i):(y(i+1)-1),x(i):x(i)+b-1)=abs(matnorm(matC2,1));
Wpo(y(i):(y(i+1)-1),x(i):x(i)+b-1)=abs(matnorm(matC2,2));

end

end
save  (['Wpo',num2str(Frq)],'Wpo')
% imagesc(abs(matnorm(CFCMAT,2)))
figure
imagesc(Wpo)
title(['Freq = ' num2str(Frq)  ] ,'fontsize',18)
%% plot rectangulars
hold on
for R=1:7
    
xx=[(2*b)*R+(shank_sapce)*(R-1)-1+0.5];
yy=1;
h=rectangle('Position',[xx,yy,shank_sapce,length(Wpo)],'FaceColor',[0 0 0.5]);
set(h, 'edgecolor','none')

end
%% plot
set(gca,'XTickLabel',[])
set(gca,'XTickLabel',[1 2 3  4 5  6 7 8 ]);
set(gca,'Xtick',[b (b)*3+shank_sapce (b)*5+shank_sapce*2 (b)*7+shank_sapce*3 (b)*9+shank_sapce*4 (b)*11+shank_sapce*5 (b)*13+shank_sapce*6 (b)*15+shank_sapce*7 ])

set(gca,'YTickLabel',[])
set(gca,'YTickLabel',[1 8 2 7 3 6 4 5]);
set(gca,'Ytick',[a/2 a/2*3 a/2*5 a/2*7 a/2*9 a/2*11 a/2*13 a/2*15])
 set(gcf,'color','w');