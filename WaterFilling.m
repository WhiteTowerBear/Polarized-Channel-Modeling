function [PowerQ,RatePA1]=WaterFilling(ChannelOri,trans_power)
% Water-filling method

[~,ChanLambda,~]=svd(ChannelOri*ChannelOri');
ChanLambda(ChanLambda==0)=[]; 
 
test1	= [1:length(ChanLambda)]./ChanLambda;
test2	= cumsum(1./ChanLambda) + trans_power;
q	= max(find(test1<test2));
if(isempty(q))
    RatePA1	= 0;
    PowerQ	= zeros(size(ChanLambda));
else
    X	= test2(q)/q;
    PowerQ	= max(X - 1./ChanLambda,0);
    RatePA1	= sum(log2(1+PowerQ.*ChanLambda));
end

if(nargout>1 & q<length(ChanLambda))
    PowerQ(length(ChanLambda))	= 0;
end

end