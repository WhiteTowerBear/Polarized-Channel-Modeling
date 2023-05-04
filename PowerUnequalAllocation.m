function [RatePA3_sum,RatePA3_AllUsers]=PowerUnequalAllocation(SNR_value,noise,users,receive,ChannelBD)
% PA3: The unequal power allocation is assigned to the three polarizations,
% then the water-filling is performed in all polarized channels

trans_power=noise*10^(0.1*SNR_value); 
for index_user=1:users.num 
    temp_powerX=norm(ChannelBD.XX((index_user-1)*receive.totalNum...
        +1:index_user*receive.totalNum,:),'fro');
    temp_powerY=norm(ChannelBD.YY((index_user-1)*receive.totalNum...
        +1:index_user*receive.totalNum,:),'fro');
    temp_powerZ=norm(ChannelBD.ZZ((index_user-1)*receive.totalNum...
        +1:index_user*receive.totalNum,:),'fro');
    ChannPower=diag([temp_powerX,temp_powerY,temp_powerZ]);  
    [PowerQ,~]=WaterFilling(ChannPower,1);
    trans_PowerFirst=PowerQ*trans_power;

    [PowerG_AllUsers_XX(:,index_user),RatePA3_AllUsers_XX(index_user,1)]...
            =WaterFilling(ChannelBD.XX((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:),trans_PowerFirst(1)); 
    [PowerG_AllUsers_YY(:,index_user),RatePA3_AllUsers_YY(index_user,1)]...
            =WaterFilling(ChannelBD.YY((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:),trans_PowerFirst(2)); 
    [PowerG_AllUsers_ZZ(:,index_user),RatePA3_AllUsers_ZZ(index_user,1)]...
            =WaterFilling(ChannelBD.ZZ((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:),trans_PowerFirst(3));
    PowerG_AllUsers=[PowerG_AllUsers_XX,PowerG_AllUsers_YY,PowerG_AllUsers_ZZ];
    RatePA3_AllUsers(index_user,1)=(RatePA3_AllUsers_XX(index_user,1)...
        +RatePA3_AllUsers_YY(index_user,1)+RatePA3_AllUsers_ZZ(index_user,1));
     
end

RatePA3_sum=sum(RatePA3_AllUsers); 

end