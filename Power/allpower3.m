function allpower3(filename,Field_edge,frq1,frq2,EEGsamplerate,eeg,cellN,Folder,CN)

%% Define the PF Ndx
[T_PF] = PF_Ndx( filename,Field_edge,CN);
% [T_PF2] = PF_Ndx2( filename,Field_edge,CN,smooth,Matrix2); more precise
Ndx_q1=round([T_PF(:,1) (T_PF(:,1)+(T_PF(:,2)-T_PF(:,1))./2)].*EEGsamplerate); % Ndx_q1
Ndx_q2=round([ (T_PF(:,1)+(T_PF(:,2)-T_PF(:,1))./2) T_PF(:,2)].*EEGsamplerate); % Ndx_q2
Ndx_q3=round([T_PF(:,2) (T_PF(:,2)+(T_PF(:,3)-T_PF(:,2))./2)].*EEGsamplerate); % Ndx_q3
Ndx_q4=round([(T_PF(:,2)+(T_PF(:,3)-T_PF(:,2))./2) T_PF(:,3) ].*EEGsamplerate); % Ndx_q4
Ndx_q5=round([T_PF(:,1) T_PF(:,2)].*EEGsamplerate); % Ndx_hf1
Ndx_q6=round([T_PF(:,2) T_PF(:,3)].*EEGsamplerate); % Ndx_hf2
Ndx_q7=round([T_PF(:,1) T_PF(:,3)].*EEGsamplerate); % Ndx_F
%% 
% u1=round(T_PF(1,1).*EEGsamplerate);
u2=round(T_PF(end,3).*EEGsamplerate)+1000;


freq = frq1:frq2;
scale = frq2scal(freq,EEGsamplerate);


% for jj=ch1:ch2
tic
jj=1; % for one cell
tic
PowerT3=[];
S=[];
eegch=eeg(:,jj);
eegch=eegch(1:u2+1000);

S = cwt(eegch,scale,'morl');
PowerT3 = (envelop(S.*S))';
Power=PowerT3;
%%
 cd (Folder)
 for q=1:7
 PP=[];   
 Ndx=eval(['Ndx_q' num2str(q)]);
 %%
for i=1:length(Ndx)
PP(i,:)=mean(Power(Ndx(i,1):Ndx(i,2),:),1); %eslah shavad, mean pp baraye ha cell save shavad
end
%% saving
if q < 5
    save (['PP' num2str(cellN) '_'  num2str(q)], 'PP');
    
elseif  q == 5
    save (['P_hf' num2str(cellN) '_'  num2str(1)], 'PP');
elseif  q == 6
    save (['P_hf' num2str(cellN) '_'  num2str(2)], 'PP');
elseif q==7
    save (['P_F' num2str(cellN) ], 'PP');
end 

 end
%% b)-2 For cancatted windows

 numrows=1000;
 
 
 for q=1:7
 Ndx=eval(['Ndx_q' num2str(q)]);
%%
C=[];
 ca=[];
 as=[];
 B=[];
for xx=1:length(Ndx)
 
 as=Power(Ndx(xx,1):Ndx(xx,2),:);
 B = imresize(as,[numrows frq2-frq1+1]);
 C=cat(1,B,C);

end
%% saving
if q < 5
    save (['CP' num2str(cellN) '_'  num2str(q)], 'C');
    
elseif  q == 5
    save (['CP_hf' num2str(cellN) '_'  num2str(1)], 'C');
elseif  q == 6
    save (['CP_hf' num2str(cellN) '_'  num2str(2)], 'C');
elseif q==7
    save (['CP_F' num2str(cellN) ], 'C');
end 



end
cd ..
% end
%%





end

