
clc
clear
% dir = 'A:\OptoMECLEC\OML19\day2\';
dirData = 'Y:\OptoMECLEC\';
animal  = {'OML18','OML18','OML18','OML18','OML19','OML19'};
sessions  = {'day1','day2','day4','day5','day2','day3'};
ses=6;
dirses = [dirData animal{ses} '\'  sessions{ses}];
cd(dirses);
filename=[animal{ses} sessions{ses}]
load([filename '.mat'])
dirsave = 'Metrics\';
load([ sessions{ses} '_Pos.mat']);
load('thetaP.mat')

cd PCsAnt
load('PCs.mat')

tic
% speed thrs
for i = 1:length(Pos.ON)
    if Pos.ON(i,3) < 0.2
       Pos.ON(i,2) = NaN;
    end
end
for i = 1:length(Pos.Off)
    if Pos.Off(i,3) < 0.2
       Pos.Off(i,2) = NaN;
    end
end   
figure; plot(Pos.ON(:,1),Pos.ON(:,2),'.b');hold on;  plot(Pos.Off(:,1),Pos.Off(:,2),'.r');hold on;
posTrials{1} = Pos.ON(:,1:2);  
posTrials{2} = Pos.Off(:,1:2);  

for i=1:size(curve,2)
    for j=1:2
        for k=1:length(stats{i}{j}.peak) 
            if stats{i}{j}.peak(k) ~= 0
                boundaries{i}{j}(k,1)= curve{i}{j}.x(stats{i}{j}.fieldX(k,1));
                boundaries{i}{j}(k,2)= curve{i}{j}.x(stats{i}{j}.fieldX(k,2));
            else
                boundaries{i}{j}(k,1)= NaN;
                boundaries{i}{j}(k,2)= NaN;
            end
        end   
    end
end

 % per each PF individually
for i=1:size(curve,2)
    
   for j=1
          for k=1:length(stats{i}{j}.x) %size(boundaries{i}{j},1) Number of the fields
               if ~isnan (stats{i}{j}.x(k))%(boundaries{i}{j}(k,1)) % for each PF not removed 
                   i
                [dataPPon{k,i},statsPPon{k,i}] = PhasePrecession(posTrials{j},spikes.Cellinfo(i).times,thetaP,'boundaries',boundaries{i}{j}(k,:));
              else
               dataPPon{k,i}=nan;statsPPon{k,i}=nan;
              end
          end
   end
end

for i=1:size(curve,2)
    
   for j=2
          for k=1:length(stats{i}{j}.x) %size(boundaries{i}{j},1)
              if ~isnan (stats{i}{j}.x(k))%(boundaries{i}{j}(k,1)) % for each PF not removed
                  i
                [dataPPoff{k,i},statsPPoff{k,i}] = PhasePrecession(posTrials{j},spikes.Cellinfo(i).times,thetaP,'boundaries',boundaries{i}{j}(k,:));
               else
               dataPPoff{k,i}=nan;statsPPoff{k,i}=nan;
              end
          end
   end   
end
cd ..

save([dirsave 'PhasePrec.mat'],'dataPPon','statsPPon','dataPPoff','statsPPoff');
toc
