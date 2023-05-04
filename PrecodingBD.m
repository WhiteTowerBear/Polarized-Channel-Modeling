function  [PrecodeMatrixF]=PrecodingBD(ChannelOri,receive,users)
% This function is the BD (block diagonalization) precoding

PrecodeMatrixF=[];
for index_user=1:users.num
    tempChannelPart=ChannelOri;
    tempChannelPart((index_user-1)*receive.totalNum+1:index_user...
        *receive.totalNum,:)=[];
    [~,~,barV]=svd(tempChannelPart);
    barVTNull=barV(:,rank(tempChannelPart)+1:end);
    dotChannel=ChannelOri((index_user-1)*receive.totalNum+1:index_user...
        *receive.totalNum,:)*barVTNull;
    [~,~,dotV]=svd(dotChannel); 
    dotVTSig=dotV(:,1:rank(dotChannel));
    PrecodeMatrixF =[PrecodeMatrixF,barVTNull*dotVTSig];

    % test result
    test=tempChannelPart*barVTNull;
    test2=ChannelOri((index_user-1)*receive.totalNum+1:index_user...
        *receive.totalNum,:)*barVTNull;
    test3=test2*dotVTSig;
end
end



