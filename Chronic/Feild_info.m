function [Field_edge,rows] = Feild_info( filename,CellVector,smooth,delta,hw)
% delta=0.2; % place_field width
% hw=35; % place_field ranges
% load('SFA4_S3_TRD1.mat')
% filename='SFA4_S3_TRD1';
% smooth=10;
% delta=0.25;
% hw=40;
% CellVector=Cell_group.C(:,3);


load([filename '.mat']);
load([filename 'PlaceField.mat'])
nPFbin=100;
[a1,a2]=size(CellVector);
for i=1:a1

    Celln=CellVector(i);
    smoothT=smooth1D(xttsc(:,2),smooth,1);
    smoothC=smooth1D(xttsc(:,Celln+4),smooth,1);
    ncell=length(smoothC(1,:));
    rate=smoothC./repmat(smoothT,1,ncell);
    xtsr=[xttsc(:,[1 3 4]) rate];
    rr=xtsr(:,4);
    matC=reshape(rr,nPFbin,length(rr)/nPFbin)';
    matC=matnorm(matC,2);
    rows(i,:)=mean(matC,1);
    
end

%%
[s1,~]=size(CellVector);

Field_edge=[];
for h=1:s1;
    x1=[];
    x2=[];
    PF_C2=[];
    
[PF_C1,PF_C2]=max(rows(h,:));
Field_edge(h,2)=PF_C2;

hw2=[PF_C2:PF_C2+hw];
hw1=[PF_C2-hw:PF_C2];

if PF_C2 > 100-hw
hw3=[PF_C2:100 1:(hw-length(PF_C2:100))];
x1 = hw1(find(rows(h,hw1) < delta*PF_C1, 1, 'last'))-1;
x2 =hw3(find(rows(h,hw3) < delta*PF_C1, 1, 'first'));

elseif PF_C2 < hw
hw4=[PF_C2:-1:1 100:-1:(100-[hw-length(1:PF_C2)])];
x1 = hw4(find(rows(h,hw4) < delta*PF_C1, 1, 'first'))-1;
x2 = hw2(find(rows(h,hw2) < delta*PF_C1, 1, 'first'));

else
x1 = hw1(find(rows(h,hw1) < delta*PF_C1, 1, 'last'))-1;
x2 = hw2(find(rows(h,hw2) < delta*PF_C1, 1, 'first'));

end
Field_edge(h,1)=x1;
Field_edge(h,3)=x2;
end


%%
for i=1:s1;
    if Field_edge(i,2) < Field_edge(i,1) 
        
       Field_edge(i,4)= Field_edge(i,1);
       Field_edge(i,5)= Field_edge(i,2)+100;
       Field_edge(i,6)= Field_edge(i,3)+100;
    
    elseif Field_edge(i,2) > Field_edge(i,3)
    
       Field_edge(i,4)= Field_edge(i,1);
       Field_edge(i,5)= Field_edge(i,2);
       Field_edge(i,6)= Field_edge(i,3)+100 ;
       
    else
       Field_edge(i,4)= Field_edge(i,1);
       Field_edge(i,5)= Field_edge(i,2);
       Field_edge(i,6)= Field_edge(i,3); 
    end  
        
        
        
end


end

