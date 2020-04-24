%%
function positionDecodingMaxCorr=Dav_positionDecodingMaxCorr2(spk_trains,position,smoothingRange)
%%

positionDecodingMaxCorr=[];
for cond =1:size(position,2)
    
    positionDecodingMaxCorr{cond}.results = table;
    r = randperm(length(spk_trains{cond}));
    for iter = 1:length(r)
        r = circshift(r,1);
%         for wind = smoothingRange
            wind = smoothingRange; 
            rates_trains_smooth_train=[];
            rates_trains_smooth_test=[];
            position_test=[];
            
            nCells=size(spk_trains{1,cond}{1,1},1);
        for cell = 1:nCells
            rates_trains_smooth = [];
            position_train = [];


            for t = 1:length(spk_trains{cond})
                rates_trains_smooth = [rates_trains_smooth; smooth(spk_trains{cond}{r(t)}(cell,:),wind)*wind];
                position_train = [position_train; position{cond}{r(t)}'];
            end
            
            rates_trains_smooth_train(cell,:) = rates_trains_smooth;          
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
        
               %% put data into struct/table
        struct.tau = wind;
        struct.condition = cond;
        struct.iter = iter;
        struct.trialOrder = r;
        positionDecodingMaxCorr{cond}.results = [positionDecodingMaxCorr{cond}.results;struct2table(struct)];
        
    end
    
end


