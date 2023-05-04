function MMSE_capacity=Precoding_MMSE(ChannelFullPolar,noise,SNR_value)
 
%% %MMSE for Perfect channels
P=noise*10^(0.1*SNR_value); 
 

[U,S,V] = svd(ChannelFullPolar*ChannelFullPolar'+ noise * eye(size(ChannelFullPolar,1))); 
T=S;
T(find(S~=0)) = 1./S(find(S~=0));

BFMatrix =ChannelFullPolar'*  V * T' * U';
weights = ones(size(ChannelFullPolar,1),1);
 

BFMatrixDiagNorm=diag(sqrt(BFMatrix'*BFMatrix))';
BFMatrix = BFMatrix./repmat(BFMatrixDiagNorm,size(BFMatrix,1),1);
 
rhos = diag(abs(ChannelFullPolar*BFMatrix).^2)';
powerAllocationwMMSE = functionHeuristicPowerAllocation(rhos,P,weights);
 
BFMatrix = kron(sqrt(powerAllocationwMMSE),ones(size(BFMatrix,1),1)).*BFMatrix;


H_G_MMSE=ChannelFullPolar*BFMatrix;
MMSE_channelGains = abs(H_G_MMSE).^2;
MMSE_signalGains = diag(MMSE_channelGains);
MMSE_interferenceGains = sum(MMSE_channelGains,2)-MMSE_signalGains;
MMSE_rates = log2(1+MMSE_signalGains./(MMSE_interferenceGains+noise));
% MMSE_rates = log2(1+powerAllocationwMMSE./(var_noise));

MMSE_capacity=sum(MMSE_rates);
end