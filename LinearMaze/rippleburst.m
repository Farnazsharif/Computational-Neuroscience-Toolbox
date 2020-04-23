
% This will locate ripple burst and quantify and isolate them
% USAGE:
%     ripplesISO(filename,cluster_group)
%
% INPUTS:
%     filename: base name of the .kwik and .kwx file to be read
%     cluster_group: vector with the "name" of the different
%     cluster_group, it can either start with 0 or 1,
%     %example: filename = 'AB1'
%               cluster_group = 1:8
% 
% Aza 2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AB4 01032013
filebase = ['D:\aza\analysis\data' '\AB4' '\01032013' ];
ratname='\AB4_030113';
cluster_group=1:7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load cells organize per shank (cluster group)
[spiket, spikeind, numclus, iEleClu] = ReadCluRes([filebase '\clures' '\jdl512_1_m_01032013'],cluster_group);
%Reading neurons classification and spikes
whowhat = readwhowhat([filebase '\general' ratname '.xlsx']);
clc

%PC analysis for... ca1
ca1 = whowhat.id(find(whowhat.region==1));
ca1pyrlayer = ca1(find(whowhat.layer(ca1)==1));
ca1pyrlayerpyrcell = ca1pyrlayer(find(whowhat.type(ca1pyrlayer)==1));
ca1pyrlayerintcell = ca1pyrlayer(find(whowhat.type(ca1pyrlayer)==2));

%PC analysis for... ca2
ca2 = whowhat.id(find(whowhat.region==2));
ca2pyrlayer = ca2;
ca2pyrlayerpyrcell = ca2pyrlayer(find(whowhat.type(ca2pyrlayer)==1));
ca2pyrlayerintcell = ca2pyrlayer(find(whowhat.type(ca2pyrlayer)==2));

%PC analysis for... ca3
ca3 = whowhat.id(find(whowhat.region==3));
ca3pyrlayer = ca3; %todo int and pyr
ca3pyrlayerpyrcell = ca3pyrlayer(find(whowhat.type(ca3pyrlayer)==1));
ca3pyrlayerintcell = ca3pyrlayer(find(whowhat.type(ca3pyrlayer)==2));

%% load ripples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SR = 20000; Fs = 1250;
%Load ripples for diff region...
load('Z:\AYA\data\AB4\jdl512_1_m\01032013\ripples3SD\AB4_030113ripsALL_3SD.mat');
%ca1
ca1rips(:,1) = floor(ripsCA13.trigs);ca1rips(:,2:3) = floor(ripsCA13.edges);
%ca2
ca2rips(:,1) = floor(ripsCA2.trigs);ca2rips(:,2:3) = floor(ripsCA2.edges);
%ca3
ca3rips(:,1) = floor(ripsCA35.trigs);ca3rips(:,2:3) = floor(ripsCA35.edges);
clear ripsCA11 ripsCA12 ripsCA13 ripsCA14 ripsCA15 ripsCA2 ripsCA34 ripsCA35 ripsCA36 ripsCA37
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Timing Description...
t1 = 167*60+54; %first = -> sleep
t2 = t1+8*60+19; %second = 
t3 = t2+26*60+37; %third =
t4 = t3+8*60+38; %fourth = 

load('AB4_030113.STATES.mat');

%% Separation in diff states...

select = {ca1rips;ca2rips;ca3rips};

for ii=1:length(select)
    selectrips = select{ii};
    ripsa = selectrips(selectrips(:,1)<t1*Fs,1);ripsb = selectrips(selectrips(:,1)>t4*Fs,1);
    ripspeakSLEEP = cat(1,ripsa,ripsb);clear ripsa ripsb %peak
    ripsc = selectrips(selectrips(:,2)<t1*Fs,2);ripsd = selectrips(selectrips(:,2)>t4*Fs,2);
    ripsstartSLEEP = cat(1,ripsc,ripsd);clear ripsc ripsd %start
    ripsz = selectrips((selectrips(:,1)>t1*Fs)&(selectrips(:,1)<t4*Fs),1); %peak
    ripsy = selectrips((selectrips(:,2)>t1*Fs)&(selectrips(:,2)<t4*Fs),2); %start
    ripspeakAWAKE = ripsz; ripsstartAWAKE = ripsy;
    ripspeak = cat(1,ripspeakSLEEP,ripspeakAWAKE);
    ripsstart = cat(1,ripsstartSLEEP,ripsstartAWAKE);
    RIPSpeak.REM=[];RIPSpeak.SWS=[];RIPSpeak.QUIET=[];RIPSpeak.RUN=[];
    RIPSstart.REM=[];RIPSstart.SWS=[];RIPSstart.QUIET=[];RIPSstart.RUN=[];

    count=0;
    for i = 1:length(ripspeak)
        if ~isempty(find(ripspeak(i)==STATES.REM))
           RIPSpeak.REM = cat(1,RIPSpeak.REM,ripspeak(i)); 
        elseif ~isempty(find(ripspeak(i)==STATES.SWS))
           RIPSpeak.SWS = cat(1,RIPSpeak.SWS,ripspeak(i));
        elseif ~isempty(find(ripspeak(i)==STATES.RUN))
           RIPSpeak.RUN = cat(1,RIPSpeak.RUN,ripspeak(i));
        elseif ~isempty(find(ripspeak(i)==STATES.QUIET))
           RIPSpeak.QUIET = cat(1,RIPSpeak.QUIET,ripspeak(i));
        else
            count=count+1;
        end
    end
    if count~=0
        disp(['Missing ' num2str(count)]);
    end
    clear i

    count=0;
    for i = 1:length(ripsstart)
        if ~isempty(find(ripsstart(i)==STATES.REM))
           RIPSstart.REM = cat(1,RIPSstart.REM,ripsstart(i)); 
        elseif ~isempty(find(ripsstart(i)==STATES.SWS))
           RIPSstart.SWS = cat(1,RIPSstart.SWS,ripsstart(i));
        elseif ~isempty(find(ripspeak(i)==STATES.RUN))
           RIPSstart.RUN = cat(1,RIPSstart.RUN,ripsstart(i));
        elseif ~isempty(find(ripspeak(i)==STATES.QUIET))
           RIPSstart.QUIET = cat(1,RIPSstart.QUIET,ripspeak(i));
        else
            count=count+1;
        end
    end
    if count~=0
        disp(['Missing ' num2str(count)]);
    end
    clear i
    
    if ii==1
        CA1RIPSpeak = RIPSpeak; CA1RIPSstart = RIPSstart;
    elseif ii==2
        CA2RIPSpeak = RIPSpeak; CA2RIPSstart = RIPSstart;
    elseif ii==3
        CA3RIPSpeak = RIPSpeak; CA3RIPSstart = RIPSstart;
    end
    
end
clear ii select
clear RIPSpeak RIPSstart ripspeakSLEEP ripsstartSLEEP ripspeakAWAKE ripspeakAWAKE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Locate ripples together and delete

window = floor(300*1e-3*Fs); %tiempo inter
CA1RIPSpeakISO.SWS = CA1RIPSpeak.SWS;CA1RIPSpeakISO.QUIET = CA1RIPSpeak.QUIET;
CA2RIPSpeakISO.SWS = CA2RIPSpeak.SWS;CA2RIPSpeakISO.QUIET = CA2RIPSpeak.QUIET;
CA3RIPSpeakISO.SWS = CA3RIPSpeak.SWS;CA3RIPSpeakISO.QUIET = CA3RIPSpeak.QUIET;
%ca1rips
for i = 1:length(CA1RIPSpeakISO.SWS)-1
    t0 = CA1RIPSpeakISO.SWS(i); t1 =  CA1RIPSpeakISO.SWS(i+1);
    distance = t1 - t0;
    %if consecutive ripples too close, kill next one
    if distance < window
        CA1RIPSpeakISO.SWS(i+1) = -1;
    end
end
clear i
CA1RIPSpeakISO.SWS = CA1RIPSpeakISO.SWS(CA1RIPSpeakISO.SWS~=-1);
for i = 1:length(CA1RIPSpeakISO.QUIET)-1
    t0 = CA1RIPSpeakISO.QUIET(i); t1 =  CA1RIPSpeakISO.QUIET(i+1);
    distance = t1 - t0;
    %if consecutive ripples too close, kill next one
    if distance < window
        CA1RIPSpeakISO.QUIET(i+1) = -1;
    end
end
clear i
CA1RIPSpeakISO.QUIET = CA1RIPSpeakISO.QUIET(CA1RIPSpeakISO.QUIET~=-1);

%ca2rips
for i = 1:length(CA2RIPSpeakISO.SWS)-1
    t0 = CA2RIPSpeakISO.SWS(i); t1 =  CA2RIPSpeakISO.SWS(i+1);
    distance = t1 - t0;
    %if consecutive ripples too close, kill next one
    if distance < window
        CA2RIPSpeakISO.SWS(i+1) = -1;
    end
end
clear i
CA2RIPSpeakISO.SWS = CA2RIPSpeakISO.SWS(CA2RIPSpeakISO.SWS~=-1);
for i = 1:length(CA2RIPSpeakISO.QUIET)-1
    t0 = CA2RIPSpeakISO.QUIET(i); t1 =  CA2RIPSpeakISO.QUIET(i+1);
    distance = t1 - t0;
    %if consecutive ripples too close, kill next one
    if distance < window
        CA2RIPSpeakISO.QUIET(i+1) = -1;
    end
end
clear i
CA2RIPSpeakISO.SWS = CA2RIPSpeakISO.SWS(CA2RIPSpeakISO.SWS~=-1);

%ca3rips
for i = 1:length(CA3RIPSpeakISO.SWS)-1
    t0 = CA3RIPSpeakISO.SWS(i); t1 =  CA3RIPSpeakISO.SWS(i+1);
    distance = t1 - t0;
    %if consecutive ripples too close, kill next one
    if distance < window
        CA3RIPSpeakISO.SWS(i+1) = -1;
    end
end
clear i
CA3RIPSpeakISO.SWS = CA3RIPSpeakISO.SWS(CA3RIPSpeakISO.SWS~=-1);
for i = 1:length(CA3RIPSpeakISO.QUIET)-1
    t0 = CA3RIPSpeakISO.QUIET(i); t1 =  CA3RIPSpeakISO.QUIET(i+1);
    distance = t1 - t0;
    %if consecutive ripples too close, kill next one
    if distance < window
        CA3RIPSpeakISO.QUIET(i+1) = -1;
    end
end
clear i
CA3RIPSpeakISO.SWS = CA3RIPSpeakISO.SWS(CA3RIPSpeakISO.SWS~=-1);

clear window t0 t1
%% Count single, double and triplets

CA1RIPSpeakBURST.SWS = CA1RIPSpeak.SWS; CA1RIPSpeakBURST.QUIET = CA1RIPSpeak.QUIET;
CA2RIPSpeakBURST.SWS = CA2RIPSpeak.SWS; CA2RIPSpeakBURST.QUIET = CA2RIPSpeak.QUIET;
CA3RIPSpeakBURST.SWS = CA3RIPSpeak.SWS; CA3RIPSpeakBURST.QUIET = CA3RIPSpeak.QUIET;
ca1swsrippleburst = zeros(5,1);ca1quietrippleburst = zeros(5,1);
ca2swsrippleburst = zeros(5,1);ca2quietrippleburst = zeros(5,1);
ca3swsrippleburst = zeros(5,1);ca3quietrippleburst = zeros(5,1);

window = floor(200*1e-3*Fs); %tiempo inter

%ca1 SWS
ripselect = CA1RIPSpeakBURST.SWS;burstselect = zeros(5,1);
for i = 1:length(ripselect)-3;
    if ripselect(i) ~= -1
    t0 = ripselect(i); t1 = ripselect(i+1);
    distance = t1-t0;
    count=0; %count of burst n
    if distance < window
        count=count+1;
        distance2 = ripselect(i+2) - ripselect(i+1);
        if distance2 < window
            count = count+1;
            distance3 = ripselect(i+3) - ripselect(i+2);
            if distance3 < window
                count = count+1;
                distance4 = ripselect(i+4) - ripselect(i+3);
                if distance4 < window
                    count = count+1;
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+4) = -1;
                else
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+3) = -1;
                end
            else
                burstselect(count,1) = burstselect(count,1) + 1;
                ripselect(i:i+2) = -1;
            end
        else
            burstselect(count,1) = burstselect(count,1) + 1;
            ripselect(i:i+1) = -1;
        end
    end
    end
end
ca1swsrippleburst = burstselect;
clear distance distance2 distance3 distance4 ripselect burstselect
%ca1 QUIET
ripselect = CA1RIPSpeakBURST.QUIET;burstselect = zeros(5,1);
for i = 1:length(ripselect)-3;
    if ripselect(i) ~= -1
    t0 = ripselect(i); t1 = ripselect(i+1);
    distance = t1-t0;
    count=0; %count of burst n
    if distance < window
        count=count+1;
        distance2 = ripselect(i+2) - ripselect(i+1);
        if distance2 < window
            count = count+1;
            distance3 = ripselect(i+3) - ripselect(i+2);
            if distance3 < window
                count = count+1;
                distance4 = ripselect(i+4) - ripselect(i+3);
                if distance4 < window
                    count = count+1;
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+4) = -1;
                else
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+3) = -1;
                end
            else
                burstselect(count,1) = burstselect(count,1) + 1;
                ripselect(i:i+2) = -1;
            end
        else
            burstselect(count,1) = burstselect(count,1) + 1;
            ripselect(i:i+1) = -1;
        end
    end
    end
end
ca1quietrippleburst = burstselect;
clear distance distance2 distance3 distance4 ripselect burstselect

%ca2 SWS
ripselect = CA2RIPSpeakBURST.SWS;burstselect = zeros(5,1);
for i = 1:length(ripselect)-3;
    if ripselect(i) ~= -1
    t0 = ripselect(i); t1 = ripselect(i+1);
    distance = t1-t0;
    count=0; %count of burst n
    if distance < window
        count=count+1;
        distance2 = ripselect(i+2) - ripselect(i+1);
        if distance2 < window
            count = count+1;
            distance3 = ripselect(i+3) - ripselect(i+2);
            if distance3 < window
                count = count+1;
                distance4 = ripselect(i+4) - ripselect(i+3);
                if distance4 < window
                    count = count+1;
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+4) = -1;
                else
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+3) = -1;
                end
            else
                burstselect(count,1) = burstselect(count,1) + 1;
                ripselect(i:i+2) = -1;
            end
        else
            burstselect(count,1) = burstselect(count,1) + 1;
            ripselect(i:i+1) = -1;
        end
    end
    end
end
ca2swsrippleburst = burstselect;
clear distance distance2 distance3 distance4 ripselect burstselect
%ca2 QUIET
ripselect = CA2RIPSpeakBURST.QUIET;burstselect = zeros(5,1);
for i = 1:length(ripselect)-3;
    if ripselect(i) ~= -1
    t0 = ripselect(i); t1 = ripselect(i+1);
    distance = t1-t0;
    count=0; %count of burst n
    if distance < window
        count=count+1;
        distance2 = ripselect(i+2) - ripselect(i+1);
        if distance2 < window
            count = count+1;
            distance3 = ripselect(i+3) - ripselect(i+2);
            if distance3 < window
                count = count+1;
                distance4 = ripselect(i+4) - ripselect(i+3);
                if distance4 < window
                    count = count+1;
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+4) = -1;
                else
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+3) = -1;
                end
            else
                burstselect(count,1) = burstselect(count,1) + 1;
                ripselect(i:i+2) = -1;
            end
        else
            burstselect(count,1) = burstselect(count,1) + 1;
            ripselect(i:i+1) = -1;
        end
    end
    end
end
ca2quietrippleburst = burstselect;
clear distance distance2 distance3 distance4 ripselect burstselect

%ca3 SWS
ripselect = CA3RIPSpeakBURST.SWS;burstselect = zeros(5,1);
for i = 1:length(ripselect)-3;
    if ripselect(i) ~= -1
    t0 = ripselect(i); t1 = ripselect(i+1);
    distance = t1-t0;
    count=0; %count of burst n
    if distance < window
        count=count+1;
        distance2 = ripselect(i+2) - ripselect(i+1);
        if distance2 < window
            count = count+1;
            distance3 = ripselect(i+3) - ripselect(i+2);
            if distance3 < window
                count = count+1;
                distance4 = ripselect(i+4) - ripselect(i+3);
                if distance4 < window
                    count = count+1;
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+4) = -1;
                else
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+3) = -1;
                end
            else
                burstselect(count,1) = burstselect(count,1) + 1;
                ripselect(i:i+2) = -1;
            end
        else
            burstselect(count,1) = burstselect(count,1) + 1;
            ripselect(i:i+1) = -1;
        end
    end
    end
end
ca3swsrippleburst = burstselect;
clear distance distance2 distance3 distance4 ripselect burstselect
%ca3 QUIET
ripselect = CA3RIPSpeakBURST.QUIET;burstselect = zeros(5,1);
for i = 1:length(ripselect)-3;
    if ripselect(i) ~= -1
    t0 = ripselect(i); t1 = ripselect(i+1);
    distance = t1-t0;
    count=0; %count of burst n
    if distance < window
        count=count+1;
        distance2 = ripselect(i+2) - ripselect(i+1);
        if distance2 < window
            count = count+1;
            distance3 = ripselect(i+3) - ripselect(i+2);
            if distance3 < window
                count = count+1;
                distance4 = ripselect(i+4) - ripselect(i+3);
                if distance4 < window
                    count = count+1;
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+4) = -1;
                else
                    burstselect(count,1) = burstselect(count,1) + 1;
                    ripselect(i:i+3) = -1;
                end
            else
                burstselect(count,1) = burstselect(count,1) + 1;
                ripselect(i:i+2) = -1;
            end
        else
            burstselect(count,1) = burstselect(count,1) + 1;
            ripselect(i:i+1) = -1;
        end
    end
    end
end
ca3quietrippleburst = burstselect;
clear distance distance2 distance3 distance4 ripselect burstselect

%%
clear ca1swsrippleburst ca1quietrippleburst ca2swsrippleburst ca2quietrippleburst ca3swsrippleburst ca3quietrippleburst
