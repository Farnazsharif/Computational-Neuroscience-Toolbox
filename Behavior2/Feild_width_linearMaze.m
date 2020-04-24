function [Field_edge,rows,PF] = Feild_width_linearMaze(CellVector,smooth,delta,Matrix,tr,b)
% Matrix=Rate_Matrix;
% CellVector=1;
% smooth=1;
% tr=TR_all;

trials_port=tr;
nd_int=[];
for i=1:length(trials_port)
    nd_int=[nd_int;find(Matrix(:,4)==trials_port(i))];
end

Matrix=Matrix(nd_int,:);

nPFbin=max(Matrix(:,1));
[a1,~]=size(CellVector);

% for i=1:a1
   
%     Celln=CellVector(i);
    Celln=CellVector;
    smoothT=smooth1D(Matrix(:,2),smooth,1);
    smoothC=smooth1D(Matrix(:,Celln+b),smooth,1);
    ncell=length(smoothC(1,:));
    rate=smoothC./repmat(smoothT,1,ncell);
    xtsr=[Matrix(:,[1 3 4]) rate];
    rr=xtsr(:,4);
    matC=reshape(rr,nPFbin,length(rr)/nPFbin)';
    MatC=matC;
    
%   PF(i,:,:)=MatC;  
%   rows(i,:)=mean(MatC,1);
  
  PF=MatC;  
  rows=mean(MatC,1);
    
% end
imresiz_magnitude=nPFbin;
L= nPFbin*imresiz_magnitude;
C=imresize(rows,[a1 L],'lanczos3');
u=[C,C,C];

%%
[s1,~]=size(CellVector);
MS1=[];
MS2=[];
Field_edge=[];
for h=1:s1
    M=u(h,:);
    [I,J]=sort(M,'descend');
    Field_edge(h,2)=J(2);
    [~,b1] = find(M(1:J(2)) < delta*I(1),1,'last');
    if isempty(b1)==1;
        b1=NaN;
        MS1(h,:)=NaN(1,50);
        MS2(h,:)=NaN(1,50);
        Field_edge(h,1)=b1;
        Field_edge(h,3)=b1;
    else
        Field_edge(h,1)=b1;
%         MS1(h,:)=imresize( [M(Field_edge(h,1):J(2))]  , [1 50],'lanczos3') ;
        [~,b2] = find(M(J(2):end) < delta*I(1),1,'first');
        Field_edge(h,3)=b2-1+J(2);
%         MS2(h,:)= imresize( [M(J(2):Field_edge(h,3))] , [1 50],'lanczos3');
    end
end
% J-100
Field_edge=(Field_edge./imresiz_magnitude);
Field_edge2=Field_edge-nPFbin;

[a,~]=size(CellVector);

%%
x1=[];
x2=[];
x3=[];
for j=1:a;
    
    if isnan(Field_edge(j,1))==0;
        if Field_edge2(j,1)<0;
            x1(j,:)=[nPFbin+Field_edge2(j,1) nPFbin]; x2(j,:)=[0 Field_edge2(j,2)]; x3(j,:)=[Field_edge2(j,2) Field_edge2(j,3)];
        elseif Field_edge2(j,3)>nPFbin;
            x1(j,:)=[Field_edge2(j,1) Field_edge2(j,2)]; x2(j,:)=[Field_edge2(j,2) nPFbin]; x3(j,:)=[0 Field_edge2(j,3)-nPFbin];
        else isnan(Field_edge2(j,1))==0;
            x1(j,:)=[Field_edge2(j,1) Field_edge2(j,2)]; x2(j,:)=[Field_edge2(j,2) Field_edge2(j,3)]; x3(j,:)=[nan nan];
        end
    else
        x1(j,:)=[nan nan]; x2(j,:)=[nan nan]; x3(j,:)=[nan nan];
    end
    
end

Field_edge=[Field_edge x1 x2 x3];

end

