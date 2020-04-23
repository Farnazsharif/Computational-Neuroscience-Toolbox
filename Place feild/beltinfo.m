 
function beltinfo(filename)

%%save belt info
% Beltinfo.object_bgn=[0  25  50  85 110 140];
% Beltinfo.object_end=[5  30  60  90 115 150];
% Beltinfo.objectC={'y','b','w','b','y','w'};
% Beltinfo.object_ID={'Yellw Flag','Blue tube', 'sharp edged plastic', 'Blue tube', 'Yellw Flag','sharp edged plastic'};
% Beltinfo.length=230;
% Beltinfo.whole=41; 
% save('TRD1_Belt', 'Beltinfo' )

% Beltinfo.object_bgn=[0  25  60  80 115 145];
% Beltinfo.object_end=[5  35  65  85 125 150];
% Beltinfo.objectC={'g','w','g','y','w','y'};
% Beltinfo.object_ID={'Glues','white Velcros', 'Glues', 'Yellow tubes', 'white Velcros','Yellow tubes'};
% Beltinfo.length=230;
% Beltinfo.whole=41; 
% save('TRD2_Belt', 'Beltinfo' )

% Beltinfo.object_bgn=[0  35  65  90 120 145];
% Beltinfo.object_end=[10  40  70  100 125 150];
% Beltinfo.objectC={'g','m','w','g','w','m'};
% Beltinfo.object_ID={'white spine','colord wires', 'plastic', 'white spine', 'plastic','colord wires'};
% Beltinfo.length=230;
% Beltinfo.whole=41; 
% save('TRD3_Belt', 'Beltinfo' )


%%
%display belt info
fb = axes('Position', [0.13, .91, .775, .06]);

w=Beltinfo.object_end-Beltinfo.object_bgn;
y = zeros(length(Beltinfo.object_bgn),1);
dy = ones(length(Beltinfo.object_bgn),1).*10;

hold on

for ii=1:length(Beltinfo.object_bgn)  
    
    u=rectangle('position',[Beltinfo.object_bgn(ii) y(ii) w(ii) dy(ii)],'FaceColor',Beltinfo.objectC{ii});
end

set(gca,'YTickLabel',[])
set(gca,'Color',[0.8,0.8,0.8])
set(gca,'XTickLabel',[])
% set(gca,'XTickLabel',[0 16 29  48 63  81.5 93 115 130  149  162  181 193 200])
% set(gca,'Xtick',[0 16 29  48 63  81.5 93 115 130  149  162  179 190 200])

hold off

% pink code [1,0.4,0.6]