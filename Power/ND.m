function [ndx]=ND(I,w)

ndx(1,1)=find(I>=w(1), 1, 'first');
ndx(1,2)=find(I>w(end), 1, 'first')-1;

end