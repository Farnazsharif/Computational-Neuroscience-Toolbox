function [slope_All,slope_set1,slope_set2,columns ]=Evolution (PF,tr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


[a1,~,~]=size(PF);
for i=6
    
    p=[];X1=[];Y1=[];Mat=[];
    
    Mat(:,:)=PF(i,:,:); 
    columns(:,i)=mean(Mat,2);

    
    X1=(1:tr(3))';
    Y1=columns(X1,i);
    P = polyfit(X1,Y1,1);
    slope_All(i) = P(1);
    
    p=[];X1=[];Y1=[];
    X1=(1:tr(2))';
    Y1=columns(X1,i);
    P=polyfit(X1,Y1,1); 
    slope_set1(i) = P(1);
    
    p=[];X1=[];Y1=[];
    X1=(tr(2)+1:tr(3))';
    Y1=columns(X1,i);
    P=polyfit(X1,Y1,1); 
    slope_set2(i) = P(1);
    
%     M=[Mat,Mat,Mat]; 
%     In_PF=M(:,Field_edge(i,1):Field_edge(i,3));
%     columns(:,i)=length(In_PF);
 
end

end

