%%
function positionDecodingMaxCorr=Dav_positionDecodingMaxCorr(phase_trains,spk_trains,position,smoothingRange)
%%
clc
smoothingRange=10;
position{1,1}=position_C{1,5};
phase_trains=[];
spk_trains=[];8
CN=1;
for i=1:28
phase_trains{1,1}{1,i}=phase_trains_C{1,5}{1,i}(CN,:);
spk_trains{1,1}{1,i}=spk_trains_C{1,5}{1,i}(CN,:);
end
size(phase_trains)
size(spk_trains)
size(position)

%%

positionDecodingMaxCorr=[];
for cond =1:size(position,2)
    
    positionDecodingMaxCorr{cond}.results = table;
    r = randperm(length(phase_trains{cond}));
    for iter = 1:length(r)
        r = circshift(r,1); % displacing the trials
%         for wind = smoothingRange
            wind = smoothingRange; 
            phase_trains_smooth_train=[];
            rates_trains_smooth_train=[];
            phase_trains_smooth_test=[];
            rates_trains_smooth_test=[];
            position_test=[];
            
            nCells=size(spk_trains{1,cond}{1,1},1)
        for cell = 1:nCells
            phase_trains_smooth = [];
            rates_trains_smooth = [];
            position_train = [];


            for t = 1:length(phase_trains{cond})
                phase_trains_smooth = [phase_trains_smooth; circ_smoothTS(phase_trains{cond}{r(t)}(cell,:),wind,'method','mean','exclude',0)];
                rates_trains_smooth = [rates_trains_smooth; smooth(spk_trains{cond}{r(t)}(cell,:),wind)*wind];
                position_train = [position_train; position{cond}{r(t)}'];% Finding the bin numbers in the displaced trials
            end
            
            phase_trains_smooth_train(cell,:) = phase_trains_smooth;
            rates_trains_smooth_train(cell,:) = rates_trains_smooth;          
            phase_trains_smooth_test(cell,:)= [circ_smoothTS(phase_trains{cond}{r(end)}(cell,:),wind,'method','mean','exclude',0)];
            rates_trains_smooth_test(cell,:) = [smooth(spk_trains{cond}{r(end)}(cell,:),wind)*wind];
            
        end
        position_test = position{cond}{r(end)};

        %% rate coding model
        cl = max_correlation_coefficient_CL;
        cl = train(cl,rates_trains_smooth_train,round(position_train));
        
        yfit_rate=[];
        for ts = 1:length(position_test)
           yfit_rate(ts) = test(cl,round(rates_trains_smooth_test(:,ts))); 
        end
        struct.mse_rate = mean((yfit_rate-position_test).^2);
        
        % chance rate
        rr = randperm(length(position_train));
        rrr = randperm(length(position_test));
        cl = max_correlation_coefficient_CL;
        cl = train(cl,rates_trains_smooth_train(:,rr),round(position_train));
        
        yfit_chance_rate=[];
        for ts = 1:length(position_test)
           yfit_chance_rate(ts) = test(cl,round(rates_trains_smooth_test(:,rrr(ts)))); 
        end
        struct.mse_chance_rate = mean((yfit_chance_rate-position_test).^2);
        
        %% phase coding model
        
        % discretize phase_trains here...
        phase_trains_smooth_train_cos = cos(phase_trains_smooth_train);
        phase_trains_smooth_train_sin = sin(phase_trains_smooth_train);
        phase_trains_smooth_test_cos = cos(phase_trains_smooth_test);
        phase_trains_smooth_test_sin = sin(phase_trains_smooth_test);
        
        phase_trains_smooth_train(phase_trains_smooth_train==0)=nan;
        phase_trains_smooth_test(phase_trains_smooth_test==0)=nan;
        phase_trains_smooth_train = discretize(phase_trains_smooth_train,-pi:.1:pi);
        phase_trains_smooth_test = discretize(phase_trains_smooth_test,-pi:.1:pi);
        phase_trains_smooth_train(isnan(phase_trains_smooth_train))=0;
        phase_trains_smooth_test(isnan(phase_trains_smooth_test))=0;
        
        phase_trains_smooth_train_cos(phase_trains_smooth_train_cos==0)=nan;
        phase_trains_smooth_test_cos(phase_trains_smooth_test_cos==0)=nan;
        phase_trains_smooth_train_cos = discretize(phase_trains_smooth_train_cos,-1:.1:1);
        phase_trains_smooth_test_cos = discretize(phase_trains_smooth_test_cos,-1:.1:1);
        phase_trains_smooth_train_cos(isnan(phase_trains_smooth_train_cos))=0;
        phase_trains_smooth_test_cos(isnan(phase_trains_smooth_test_cos))=0;
        
        phase_trains_smooth_train_sin(phase_trains_smooth_train_sin==0)=nan;
        phase_trains_smooth_test_sin(phase_trains_smooth_test_sin==0)=nan;
        phase_trains_smooth_train_sin = discretize(phase_trains_smooth_train_sin,-1:.1:1);
        phase_trains_smooth_test_sin = discretize(phase_trains_smooth_test_sin,-1:.1:1);
        phase_trains_smooth_train_sin(isnan(phase_trains_smooth_train_sin))=0;
        phase_trains_smooth_test_sin(isnan(phase_trains_smooth_test_sin))=0;
        
        % non-transformed circular decoding
        cl = max_correlation_coefficient_CL;
        cl = train(cl,[phase_trains_smooth_train],round(position_train));
        
        yfit_circ=[];
        for ts = 1:length(position_test)
           yfit_circ(ts) = test(cl,phase_trains_smooth_test(:,ts)); 
        end
        struct.mse_phase = mean((yfit_circ-position_test).^2);
        
        % cos and sin transformed circular decoding
        cl = max_correlation_coefficient_CL;
        cl = train(cl,phase_trains_smooth_train_cos,round(position_train));
        
        yfit_cos=[];
        for ts = 1:length(position_test)
           yfit_cos(ts) = test(cl,phase_trains_smooth_test_cos(:,ts)); 
        end
        struct.mse_phase_cos = mean((yfit_cos-position_test).^2);
        
        cl = max_correlation_coefficient_CL;
        cl = train(cl,phase_trains_smooth_train_sin,round(position_train));
        
        yfit_sin=[];
        for ts = 1:length(position_test)
           yfit_sin(ts) = test(cl,phase_trains_smooth_test_sin(:,ts)); 
        end
        struct.mse_phase_sin = mean((yfit_sin-position_test).^2);
        
        % all phase models in one...
        cl = max_correlation_coefficient_CL;
        all_phase_train = [phase_trains_smooth_train_cos;phase_trains_smooth_train_sin];
        all_phase_test =  [phase_trains_smooth_test_cos;phase_trains_smooth_test_sin];
        cl = train(cl,all_phase_train,round(position_train));
        
        yfit_circ_all=[];
        for ts = 1:length(position_test)
           yfit_circ_all(ts) = test(cl,all_phase_test(:,ts)); 
        end
        struct.mse_phase_all = mean((yfit_circ_all-position_test).^2);
        
        % chance phase
        cl = max_correlation_coefficient_CL;
        all_phase_train = [phase_trains_smooth_train_cos;phase_trains_smooth_train_sin];
        all_phase_test =  [phase_trains_smooth_test_cos;phase_trains_smooth_test_sin];
        cl = train(cl,all_phase_train(:,rr),round(position_train));
        
        yfit_chance=[];
        for ts = 1:length(position_test)
           yfit_chance(ts) = test(cl,all_phase_test(:,rrr(ts))); 
        end
        struct.mse_chance = mean((yfit_chance-position_test).^2);
               %% put data into struct/table
        struct.tau = wind;
        struct.condition = cond;
        struct.iter = iter;
        struct.trialOrder = r;
        positionDecodingMaxCorr{cond}.results = [positionDecodingMaxCorr{cond}.results;struct2table(struct)];
        
    end
end