function [RatePA4_sum,RatePA4_AllUsers]=PowerAllEqual(SNR_value,noise,users,receive,ChannelBD)
% PA4: all users are allocated with equal power
trans_power=noise*10^(0.1*SNR_value); 
PolarPower=1/3; 
trans_power=PolarPower*trans_power;
for index_user=1:users.num
    ChannelOriX=ChannelBD.XX((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:);
    [~,ChanLambdaX,~]=svd(ChannelOriX*ChannelOriX'); 
    ChanLambdaX=diag(ChanLambdaX);
    PowerQ=trans_power/length(ChanLambdaX);

    ChannelOriY=ChannelBD.YY((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:);
    [~,ChanLambdaY,~]=svd(ChannelOriY*ChannelOriY'); 
    ChanLambdaY=diag(ChanLambdaY);

    ChannelOriZ=ChannelBD.ZZ((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:);
    [~,ChanLambdaZ,~]=svd(ChannelOriZ*ChannelOriZ'); 
    ChanLambdaZ=diag(ChanLambdaZ);

    RatePA4_AllUsers(index_user,1) = sum(log2(1+PowerQ.*ChanLambdaX))+...
        sum(log2(1+PowerQ.*ChanLambdaY))+sum(log2(1+PowerQ.*ChanLambdaZ));
end

RatePA4_sum=sum(RatePA4_AllUsers); 

end