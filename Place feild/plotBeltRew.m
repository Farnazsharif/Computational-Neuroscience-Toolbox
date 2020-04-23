function plotBeltRew(filename,yspan)
load([filename '.mat']);
hold on
% yspan=12; %number of triala
    for jj=1:length(Beltinfo.object_bgn)
        if Beltinfo.object_bgn(jj) > Beltinfo.object_end(jj)
           plot([0 0 Beltinfo.object_end(jj) Beltinfo.object_end(jj) 0],[yspan+1 -1 -1 yspan+1 yspan+1],Beltinfo.objectC{jj}) 
           hold on
           plot([Beltinfo.object_bgn(jj) Beltinfo.object_bgn(jj) Beltinfo.Length  Beltinfo.Length Beltinfo.object_bgn(jj)],[yspan+1 -1 -1 yspan+1 yspan+1],Beltinfo.objectC{jj},'LineStyle',':','LineWidth',0.5)
           hold on
        else
%    
        plot([Beltinfo.object_bgn(jj) Beltinfo.object_bgn(jj) Beltinfo.object_end(jj) Beltinfo.object_end(jj) Beltinfo.object_bgn(jj)],[yspan+1 -1 -1 yspan+1 yspan+1],Beltinfo.objectC{jj},'LineStyle',':','LineWidth',0.5)
      end  
    end  
    
   Rew =Reward(filename);
   
   hold on 
   for kk=1:length(Rew)
       
    plot([Rew(kk,1) Rew(kk,1) Rew(kk,2) Rew(kk,2) Rew(kk,2)],[yspan+1 -1 -1 yspan+1 yspan+1],'Color',[1 .5 0],'LineStyle','--')
   
   end
   
   hold off 
       ylim([0 yspan])
       
       set(gca,'YTickLabel',[])
%              set(gca,'XTickLabel',[])