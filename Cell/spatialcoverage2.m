function [SpatialC]=spatialcoverage2(filename1,filename3)

% load('FM05_1.mat')
% filename3='FM05_1PlaceField';
load([filename1 '.mat'])
load([filename3 '.mat'])
binN=100;
ndxzero=find(xttsc(:,1)==1);
sets=unique(xttsc(:,4));
for i=1:length(sets)
s=find(xttsc(ndxzero,4)==sets(i));
Set(i,1)=s(1);
Set(i,2)=s(end);
end

%%
scale=10;

for j=1:length(G)
 
SpatialC(j,1)=G(j);
cc=reshape(xttsc(:,4+find(G==G(j))),binN,142);

for i =1 : length(sets)
    
e=cc(:,Set(i,1):Set(i,2));
e1=matnorm(imresize(sum(e,2)',[1 binN*scale],'lanczos3'),2);
SpatialC(j,i+1)=SpatialCoverage(e1');

end

end

