function [trans,prob,prob_uncor,pred,hiBound,loBound] = GetTransProb_Roman(rawCCG,n,binSize,alpha,varargin)
% Extract the baseline corrected CCG + spike trans probe from raw CCG

% rawCCG = spike count between reference and target spike train
% n = number of reference spikes
% bin size = the binning of the CCG
% (optional input) = intwin = time bins in which synapse should inject
% excess synchrony

% alpha = .001;
% Define integration window
if ~isempty(varargin)
    intwin = varargin{1};
else
    % When binSize is 0.0008, integration window is between 0.8 to 2.8 ms
    intwin = round(length(rawCCG)/2) + round([.0008:binSize:.0028]/binSize);
end

% Define slow, network comodulation time scale in 
% seconds and in number of binsized chunks
conv_win_sec = .015; 
conv_win_binnum = round(conv_win_sec/binSize);

% Get CCG normalized by number of reference spikes
prob_uncor = rawCCG/n;

% Calculate low freq. network comodulation
[ ~, pred ] = cch_conv_Roman( rawCCG, conv_win_binnum);

hiBound=poissinv(1-alpha,pred);
loBound=poissinv(alpha, pred);

% Baseline subtract
prob = (rawCCG(:) - pred)/n;

% prob(prob<0) = 0;

% Integrate
trans = nansum(prob(intwin));
end