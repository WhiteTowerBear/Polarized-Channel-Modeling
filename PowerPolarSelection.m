function [RatePA1_sum,RatePA1_AllUsers]=PowerPolarSelection(SNR_value,noise,users,receive,ChannelBD)
% PA1: The water-filling power allocation is performed in the best
% polarized channel
trans_power=noise*10^(0.1*SNR_value);
PolarSelection=zeros(users.num*3,1);
for index_user=1:users.num
    temp_PolSel=...
        [norm(ChannelBD.XX((index_user-1)*receive.totalNum...
        +1:index_user*receive.totalNum,:),'fro')^2,norm(ChannelBD.YY...
        ((index_user-1)*receive.totalNum+1:index_user*receive.totalNum,:)...
        ,'fro')^2,norm(ChannelBD.ZZ((index_user-1)*receive.totalNum...
        +1:index_user*receive.totalNum,:),'fro')^2]; 
    temp_PolSel(temp_PolSel~=max(temp_PolSel))=0;
    temp_PolSel(temp_PolSel==max(temp_PolSel))=1;
    % each three terms belong to one user, e.g., [0,1,0] to one user in Y
    % polarization
    PolarSelection((index_user-1)*3+1:index_user*3)=temp_PolSel; 
    

    if isequal(temp_PolSel,[1,0,0]) 
        [PowerG_AllUsers(:,index_user),RatePA1_AllUsers(index_user,1)]...
            =WaterFilling(ChannelBD.XX((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:),trans_power);
    elseif isequal(temp_PolSel,[0,1,0]) 
        [PowerG_AllUsers(:,index_user),RatePA1_AllUsers(index_user,1)]...
            =WaterFilling(ChannelBD.YY((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:),trans_power);
    elseif isequal(temp_PolSel,[0,0,1])
        [PowerG_AllUsers(:,index_user),RatePA1_AllUsers(index_user,1)]...
            =WaterFilling(ChannelBD.ZZ((index_user-1)*receive.totalNum...
                    +1:index_user*receive.totalNum,:),trans_power);
    end
end

RatePA1_sum=sum(RatePA1_AllUsers);

end