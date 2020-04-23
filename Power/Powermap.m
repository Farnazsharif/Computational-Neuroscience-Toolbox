function  Powermap( filename1,filename2, i,ch)
% [ xtsP ] = Powermap( filename, i) Farnaz
% close all
% i=1;
% filename1='FM05_1';
% filename2='PP_65';
% load('PP.mat');
load([filename2,num2str(ch) '.mat'])
P=PP(:,i);
matC=[];
%%

% close all
% matC3=[];
% B=[];
% scale=10;
% B=imresize(P,scale,'lanczos3');
% matC3=reshape(B(:,1),xbinNumber,length(B)/xbinNumber)';
% imagesc(abs(matnorm(matC3,2)));
%%
smoothC=smooth1D(P,10,1);
% smoothC=smooth(P,'moving');
xbinNumber=100;
% matC1=reshape(P,100,length(P)/xbinNumber)';
matC2=reshape(smoothC,xbinNumber,length(P)/xbinNumber)';
figure ('Position',[700 700 400 500])
gu1=.07; gu2=0.08; gu3=0.9; gu4=0.8;

subplot('Position',[gu1 gu2 gu3 gu4])
% imagesc(matnorm(matC1,2))
% figure
imagesc(abs(matnorm(matC2,2)))
% title('PowerT3 map');
xlabel('bin number (belt lenght = 100 cm)','fontweight','bold');
ylabel('trials','fontweight','bold');  

hold on
%%
load([filename1 '.mat'])
%% remove negative times
M=behav.txlrts;
trial=5;
setNumber=6;
Reward=4;
negndx=find(M(:,1)<0);
length(M);
M(negndx,:)=[];
length(M);
%%
L=zeros(length(M),1);
txrtsl=M;    
%5) making bin
    % 5-1) find position zero
    ndxzero=find([1;diff(txrtsl(:,2))]<-10);
    % 5-2) for each belt rotation, compute rate in each space bin
    xbin=[];
    trialnb=[];
    setnb=[];
    tmp=[];
    % for each cycle 
for ii = 2:length(ndxzero)
    %compute xbin edges
    xstep=txrtsl(ndxzero(ii)-1,2)/xbinNumber;
    binedges=0:xstep:txrtsl(ndxzero(ii)-1,2);
    xbin=[xbin;(1:xbinNumber)'];
    %select cycles from 0 to end
    tmp=txrtsl(ndxzero(ii-1):ndxzero(ii)-1,:);
    for jj=1:length(binedges)-1
    %trial # and set#
        rndx=find(tmp(:,2)>binedges(jj)&tmp(:,2)<=binedges(jj+1));
        if ~isempty(rndx)
            trialnb=[trialnb;tmp(rndx(end),trial)];
            setnb=[setnb;tmp(rndx(end),setNumber)];
        else
            trialnb=[trialnb;trialnb(end)];
            setnb=[setnb;setnb(end)];
        end
    end
end 

xtsP=[xbin  trialnb setnb P ];

hold on

% 2) plot Reward position 
incr=15;
sets=unique(xtsP(:,3));
setyxz=nan(length(sets),4);
for ii = 1:length(sets)
    ndx1=find(xtsP(:,3)==sets(ii));
    setyxz(ii,1)=floor(ndx1(1)/xbinNumber);
    setyxz(ii,2)=floor(ndx1(end)/xbinNumber);
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&    
    %Reward
    ndx1=find(M(:,setNumber)==sets(ii));
    tmp21=M(ndx1,:);    
    setyxz(ii,3)=mode(tmp21(find(tmp21(:,Reward)==1),2)); 
    
    %Vacuum
    setyxz(ii,4)=setyxz(ii,3)+incr;
    
end

TS=unique(M(:,trial:setNumber),'rows');
maxpos=[];
for ii=1:length(TS(:,1))
    ndx=find(M(:,trial)==TS(ii,1)& M(:,setNumber)==TS(ii,2));
    maxpos=[maxpos;max(M(ndx,2))];
end
maxpos=mode(maxpos);
setyxz(:,3:4)=setyxz(:,3:4)/maxpos*xbinNumber;

 for jj=1:length(setyxz(:,1))
        plot([setyxz(jj,3) setyxz(jj,3)],[setyxz(jj,1) setyxz(jj,2)],'w');
        plot([setyxz(jj,4) setyxz(jj,4)],[setyxz(jj,1) setyxz(jj,2)],'w');
        plot([setyxz(jj,3) setyxz(jj,4)],[setyxz(jj,1) setyxz(jj,1)],'w');
        plot([setyxz(jj,3) setyxz(jj,4)],[setyxz(jj,2) setyxz(jj,2)],'w');
 end

hold on

%3) plot Belt
% [0.126 0.075 0.5 0.8]
 gu2=0.89; gu3=(Beltinfo.object_end(end))./230-gu1; gu4=0.05;
subplot('Position', [gu1, gu2, gu3, gu4]);
  
    w=Beltinfo.object_end-Beltinfo.object_bgn;
    y = zeros(length(Beltinfo.object_bgn),1);
    dy = ones(length(Beltinfo.object_bgn),1).*10;
    hold on

    for jj=1:length(Beltinfo.object_bgn) 
        u=rectangle('position',[Beltinfo.object_bgn(jj) y(jj) w(jj) dy(jj)],'FaceColor',Beltinfo.objectC{jj});
    end

    set(gca,'YTicklabel',[])
    set(gca,'xticklabel',[])

    title(['Freq = ' num2str(i)  '  ch = ' num2str(ch)] ,'fontsize',18)
    
    hold off
    set(gcf,'color','w');
    
    
end

