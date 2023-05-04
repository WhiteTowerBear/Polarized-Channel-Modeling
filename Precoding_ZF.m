function ZF_test_capacity=Precoding_ZF(ChannelFullPolar,noise,SNR_value)

%% ZF for Perfect channels
P=noise*10^(0.1*SNR_value); 

[U,S,V] = svd(ChannelFullPolar*ChannelFullPolar'); 
T=S;
T(find(S~=0)) = 1./S(find(S~=0));
 
ChannelInv =ChannelFullPolar'* V * T' * U'; 
 

weights = ones(size(ChannelFullPolar,1),1);
[n_r,~]=size(ChannelFullPolar);

for i=1:n_r
    ChannelColNorm(1,i)=norm(ChannelInv(:,i));
    BFMatrix(:,i)=sqrt(1/n_r)*ChannelInv(:,i)/ChannelColNorm(1,i);
end

% BFMatrixDiagNorm=diag(sqrt(ChannelInv'*ChannelInv))';
% BFMatrix = sqrt(1/n_r)*ChannelInv./repmat(BFMatrixDiagNorm,size(ChannelInv,1),1);
 
rhos = diag(abs(ChannelFullPolar*BFMatrix).^2)';
powerAllocationwZFBF = functionHeuristicPowerAllocation(rhos,P,weights);
 
% ZF_rate =  log2(1+powerAllocationwZFBF./( var_noise ));

ZF_rate =  log2(1+powerAllocationwZFBF./(n_r*noise*ChannelColNorm.^2));
ZF_capacity=sum(ZF_rate);

ZF_test=log2(1+P./(n_r*noise*ChannelColNorm.^2));
ZF_test_capacity=sum(ZF_test);

end