
%% Prevoius Steps
function Movei(filename1,filename2,ch1,ch2)
% clear
%  filename1='FM05_1';
%  filename2='po_';
% [PowerT3,PP ] = allpower2( filename, frq1,frq2,ch);
% [T,Ndx]=binT(filename,xbinNumber)
load([filename1 '.mat'])
load([filename2,num2str(ch1) '.mat'])
xbinNumber=100;
[a,b]=size(PP);
for ii=ch1:ch2
ii
%% Set up the movie
% axis ij;
writerObj = VideoWriter([filename2,num2str(ii) '.avi']); 
writerObj.FrameRate = 2;
open(writerObj); 
fid = figure;
b=130;
for i=1:60   
    pause(0.5);
    figure(fid );
%     Powermap( filename, i);
    Powermap( filename1,filename2, i,ii);
    open(writerObj); 
    %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
    %end
 close
end
hold off
close(writerObj); 
end
%% Movie Test
% %% Set up some function. 
% % Sine between -2*pi and 2*pi.
% x = (10*-pi:0.1:10*pi)'; % Note the transpose.
% y = sin(x);
% fid = figure;
% hold on
% % The final plot.
% plot(x,y, '*')
% % Set up the movie.
% writerObj = VideoWriter('out.avi'); % Name it.
% writerObj.FrameRate = 1; % How many frames per second.
% open(writerObj); 
%  
% for i=1:size(y)      
%     % We just use pause but pretend you have some really complicated thing here...
%     pause(0.1);
%     figure(fid ); % Makes sure you use your desired frame.
%     plot(x(i),y(i),'or');
%  
%     %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
%         frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%         writeVideo(writerObj, frame);
%     %end
%  
% end
% hold off
% close(writerObj); % Saves the movie.
%%
