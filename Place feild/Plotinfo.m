
clear
clc


fl = axes('Position', [0.13 0.175 0.775 0.755]);

open('hoho.fig');
h_hoho=get(gca,'Children');
x_hoho=get(h_hoho,'XData');
y_hoho=get(h_hoho,'YData');

open('lolo.fig');
h_lolo=get(gca,'Children');
x_lolo=get(h_lolo,'XData');
y_lolo=get(h_lolo,'YData');

%normal gca: [0.13 0.11 0.775 0.815]

figure (3)
hold on
bar(x_hoho,y_hoho,'r')
cellfun(@plot,x_lolo,y_lolo);

  set(gca,'YAxisLocation','left');
  set(gca,'XTickLabel',[]);
  ylabel('speed');
  ylim([0 180]);
  xlim([0 200]);
  
  
  haxes1 = gca;
  haxes1_pos = get(haxes1,'Position');
  haxes2 = axes('Position',haxes1_pos,...
              'YAxisLocation','right',...
              'Color','none');


  ylabel('lickrate');
  set(gca,'XTickLabel',[]);
  ylim([0 180]);
  xlim([0 200]);
  
  close(figure(2))
  close(figure(1))