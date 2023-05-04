function [RatePA4_sum,RatePA4_AllUsers]=PowerEqualAllocation(SNR_value,noise,users,receive,ChannelBD)
% PA4: The equal power allocation is assigned to the three polarizations,
% then the water-filling is performed in all polarized channels
trans_power=noise*10^(0.1*SNR_value); 
PolarPower=1/3; 
trans_power=PolarPower*trans_power;
for index_user=1:users.num
    [PowerG_AllUsers_XX(:,index_user),RatePA1_AllUsers_XX(index_user,1)]...
            =WaterFilling(ChannelBD.XX((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:),trans_power); 
    [PowerG_AllUsers_YY(:,index_user),RatePA1_AllUsers_YY(index_user,1)]...
            =WaterFilling(ChannelBD.YY((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:),trans_power); 
    [PowerG_AllUsers_ZZ(:,index_user),RatePA1_AllUsers_ZZ(index_user,1)]...
            =WaterFilling(ChannelBD.ZZ((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:),trans_power);
    PowerG_AllUsers=[PowerG_AllUsers_XX,PowerG_AllUsers_YY,PowerG_AllUsers_ZZ];
    RatePA4_AllUsers(index_user,1)=(RatePA1_AllUsers_XX(index_user,1)...
        +RatePA1_AllUsers_YY(index_user,1)+RatePA1_AllUsers_ZZ(index_user,1));
end

RatePA4_sum=sum(RatePA4_AllUsers); 

end