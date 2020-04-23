%% Find cell numbers at each shank
ds=[];
Cell_sh=[];
name='Place_field';
dirinfo=dir(name);
[ds,~]=size(dirinfo); %%check the number of folders

fln=0;
for hh=1:ds
    if isempty(strfind(dirinfo(hh).name,'h'))==1 %% detect folder name 'shank'
        fln=fln+1;
    end
end

cd (name)

Cell_sh=[];
for sh=1:ds-fln;
    dirinfo_sh=[];
    ds_sh=[];
dirinfo_sh=dir(dirinfo(sh+fln).name);
[ds_sh,~]=size(dirinfo_sh); %%check the number of folders

 fln_sh=0;
for hh=1:ds_sh
    if isempty(strfind(dirinfo_sh(hh).name,'h'))==1 ||  isempty(strfind(dirinfo_sh(hh).name,'._'))==0 %% detect folder name 'shank'
        fln_sh=fln_sh+1;
    end
end

Cell_sh=[Cell_sh, ds_sh-fln_sh];
end

cd ..
%%
ndx=[1 9 10 11 12 13 14 15 16 2 3 4 5 6 7 8 ];
dirinfo(ndx+fln).name
Cell_N=Cell_sh(ndx);
Fullpath='C:\Users\Behnamty\Documents\Pics\FM01_2'; 
%%  Making ppt file for each shank
% https://www.mathworks.com/matlabcentral/answers/99150-is-there-an-example-of-using-matlab-to-create-powerpoint-slides
gu1_a=500;gu2_a=500;
gu1_b=350;gu2_b=280;
gu1_c=350;gu2_c=250;

Gf1_a=0;Gf2_a=20;
Gf1_b=gu1_a+100;Gf2_b=0;
Gf1_c=gu1_a+100;Gf2_c=gu2_b+10;


for k=Cell_N:-1:1; 
    k
h = actxserver('PowerPoint.Application')
% h.Visible = 1;
Presentation = h.Presentation.Add

% Select the favorit slide type : Item(7)

 blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
% Generate your favorit number of slides and  insert your figure at your favorit size in your desired position

for i=Cell_N(k):-1:1
Slide(i) = Presentation.Slides.AddSlide(1,blankSlide);

Image(i) = Slide(i).Shapes.AddPicture(([Fullpath '\Placefield_images\Shank' num2str(k) '\' 'shank' num2str(k) 'cell' num2str(i) '.eps']),'msoFalse','msoTrue',Gf1_a,Gf2_a,gu1_a,gu2_a)
Image(i+1) = Slide(i).Shapes.AddPicture(([Fullpath '\Phase_images\Shank' num2str(k) '\' 'shank' num2str(k) 'cell' num2str(i) '.eps']),'msoFalse','msoTrue',Gf1_b,Gf2_b,gu1_b,gu2_b)
Image(i+2) = Slide(i).Shapes.AddPicture(([Fullpath '\Ripples_images\Shank' num2str(k) '\' 'shank' num2str(k) 'cell' num2str(i) '.jpg']),'msoFalse','msoTrue',Gf1_c,Gf2_c,gu1_c,gu2_c)
 end
Presentation.SaveAs([Fullpath '\Shank' num2str(k) '.ppt'])
invoke(Presentation,'Close');
end

h.Quit;
h.delete;

%%  Making 1 ppt file for all  shank

gu1_a=500;gu2_a=500;
gu1_b=350;gu2_b=280;
gu1_c=350;gu2_c=250;

Gf1_a=0;Gf2_a=20;
Gf1_b=gu1_a+100;Gf2_b=0;
Gf1_c=gu1_a+100;Gf2_c=gu2_b+10;

h = actxserver('PowerPoint.Application')
h.Visible = 1;
Presentation = h.Presentation.Add;
% Select the favorit slide type : Item(7)
blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);
% Generate your favorit number of slides and  insert your figure at your favorit size in your desired position


for k=length(Cell_N):-1:1; 

for i=Cell_N(k):-1:1
   
Slide(i) = Presentation.Slides.AddSlide(1,blankSlide);

Image(i) = Slide(i).Shapes.AddPicture(([Fullpath '\Placefield_images\Shank' num2str(k) '\' 'shank' num2str(k) 'cell' num2str(i) '.eps']),'msoFalse','msoTrue',Gf1_a,Gf2_a,gu1_a,gu2_a)
Image(i+1) = Slide(i).Shapes.AddPicture(([Fullpath '\Phase_images\Shank' num2str(k) '\' 'shank' num2str(k) 'cell' num2str(i) '.eps']),'msoFalse','msoTrue',Gf1_b,Gf2_b,gu1_b,gu2_b)
Image(i+2) = Slide(i).Shapes.AddPicture(([Fullpath '\Ripples_images\Shank' num2str(k) '\' 'shank' num2str(k) 'cell' num2str(i) '.jpg']),'msoFalse','msoTrue',Gf1_c,Gf2_c,gu1_c,gu2_c)
 end

end

Presentation.SaveAs([Fullpath '\Shanks_all.ppt'])
h.Quit;
h.delete;



%% 
% if there is matlab .fig, this command convert and save it to ppt with
% high resolution

open('s.fig')
print('-dmeta','s')
Image = Slide.Shapes.AddPicture((['C:\Users\Behnamty\Documents\Pics\FM01_1\s.emf']),'msoFalse','msoTrue',Gf1_a,Gf2_a,gu1_a,gu2_a);
close all

%% This command is for changing the format of the existing eps files

for k=1:16;
    k
    s1=[];
   [s1,~]=size(dir(['Shank' num2str(k)]));
    cd (['Shank' num2str(k)])
for i=1:(s1-2)
   
eps2xxx((['shank' num2str(k) 'cell' num2str(i) '.eps']),{'jpeg'},'C:\Program Files\gs\gs9.20\bin\gswin64c.exe');

end
cd ..
end





