function appendBeltinfo_DJ(filename, ver)

if ver==7 % New object added task DJ 17-19
    
    BeltInfo.Zero = 166;
    BeltInfo.Length = 201;
    BeltInfo.object_bgn = [1 25.5 51 75.5 100.5 125.5 150.5 176];
    BeltInfo.object_end = [2 26.5 52 76.5 101.5 126.5 151.5 177];
    BeltInfo.Newobject_bgn = 84 ;
    BeltInfo.Newobject_end = 94;
    
    BeltInfo.object_ID = {'RBTube_L' 'WTube_C' 'YTube_A' 'RTube_R' 'BTube_R' '2GTube_C' 'YGTube_L' 'BTube_A'};
    %BeltInfo.Newobject_ID = {'Spines' 'Valcro' 'GTube_A'};
    switch (filename)
        case 'DJ17_1'
            BeltInfo.Newobject_ID = {'Spines'};
            BeltInfo.origin = 2:71;
            BeltInfo.New = 74:163;
            BeltInfo.New2 = -1;
            BeltInfo.NewobjectC = {'w'};
        case 'DJ17_2'
            BeltInfo.Newobject_ID = {'Valcro'};
            BeltInfo.origin = 3:71;
            BeltInfo.New = 72:129;
            BeltInfo.New2 = -1;
            BeltInfo.NewobjectC = {'w:'};
        case 'DJ18_1'
            BeltInfo.Newobject_ID = {'Spines'};
            BeltInfo.origin = 2:70;
            BeltInfo.New = 72:144;
            BeltInfo.New2 = -1;
            BeltInfo.NewobjectC = {'w'};
        case 'DJ18_2'
            BeltInfo.Newobject_ID = {'Valcro' 'GTube_A'};
            BeltInfo.origin = 3:74;
            BeltInfo.New = 76:175;
            BeltInfo.New2 = 184:215;
            BeltInfo.NewobjectC = {'w:'};
        case 'DJ19_1'
            BeltInfo.Newobject_ID = {'Spines'};
            BeltInfo.origin = 2:73;
            BeltInfo.New = 77:156;
            BeltInfo.New2 = -1;
            BeltInfo.NewobjectC = {'w'};
        case 'DJ19_2'
            BeltInfo.Newobject_ID = {'Valcro'};
            BeltInfo.origin = 2:46;
            BeltInfo.New = 48:123;
            BeltInfo.New2 = -1;
            BeltInfo.NewobjectC = {'w:'};
    end
    
    
    BeltInfo.objectC = {'r' 'b' 'g' 'c' 'c' 'b' 'r' 'g'};
    
    
elseif ver==8 % periodic  and aperiodic  DJ 21-25
    
    BeltInfo.Zero = 173;
    BeltInfo.Length = 201;
    BeltInfo.object_bgn = [10 30 50 70 90 100 119 130 170 179] ;
    BeltInfo.object_end = [11 31 51 71 91 101 120 131 171 180] ;
    
    
    BeltInfo.object_ID = {'WTube_L' 'Spine_R' 'WTube_L' 'Spine_R' 'WTube_L' 'Spine_R' 'Spine_R' 'WTube_L' 'WTube_L' 'Spine_R'};
    BeltInfo.objectC = {'r' 'g' 'r' 'g' 'r' 'g' 'g' 'r' 'r', 'g'};
    
elseif (ver == 9 || ver ==10) % reward position changing task 9 for cue based 10 for random position.
    
    BeltInfo.Zero = 173.8;
    BeltInfo.Length = 200.5;
    BeltInfo.object_bgn = [9.5 29.5 50 82 99.5 122 150 179.5] ;
    BeltInfo.object_end = [10.5 30.5 59 88 100.5 128 159 180.5] ;
    
    
    BeltInfo.object_ID = {'BTube_L' 'GYTube_R' 'Spine' 'Vertical Tubes' 'GYTube_R' 'Vertical Tubes' 'Spine' 'BTube_L'};
    BeltInfo.objectC = {'b:' 'c:' 'r' 'g' 'c:' 'g' 'r' 'b:' };
    
elseif ver==11 % Dajung's 1st task with laser
    clear BeltInfo
    BeltInfo.Zero = 160;
    BeltInfo.Length = 200;
    BeltInfo.object_bgn = [173 193 40 63 106 130];
    BeltInfo.object_end = [174 10 41 82 107 150];
    BeltInfo.object_ID = {'BigTube1' 'tubes' 'BigTube2' 'spines' 'BigTube3' 'Velcro'};
    BeltInfo.objectC = {'g:' 'm' 'g:' 'c' 'g:' 'r'};
    
    
elseif ver==12 % Dajung's 2nd task with laser and adding object
    clear BeltInfo
    BeltInfo.Zero = 65;
    BeltInfo.Length = 200;
    BeltInfo.object_bgn = [74,105,113,143,153,183,23,33,62];
    BeltInfo.object_end = [94,106,133,144,173,184,24,53,63];
    BeltInfo.Newobject_bgn = 193;
    BeltInfo.Newobject_end =  13;
    
    BeltInfo.object_ID = {'spines' 'BigTube1' 'tubes' 'BigTube2' 'spines' 'BigTube3' 'BigTube4' 'tubes' 'BigTube5' };
    BeltInfo.objectC = {'c' 'g:' 'm' 'g:'  'c'  'g:' 'g:' 'm' 'g:'};
    BeltInfo.Newobject_ID = {'spines+Tube'};
    BeltInfo.NewobjectC = {'r'};
    switch (filename)
        case 'DJ_m13_3'
           BeltInfo.New = 31; 
        case 'DJ_m14_1'
           BeltInfo.New = 74;  
       case 'DJ_m16_2'
           BeltInfo.New = 100; 
            
    end
    
    
end
if ~exist('task_version')
task_version = ver;
save(filename,'-append','task_version')
end
save(filename,'-append','BeltInfo')