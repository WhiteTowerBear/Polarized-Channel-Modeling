%========================================================a
%  Author: Wei Li;
%  Date:   26/08/ 2022,  SUTD;
%  Version: V1.0 
%  Note: This is polarized HMIMOS channel modeling in near
%        field using Green's function, the precoding design
%        and power allocation. 

% The relevant paper is L. Wei et al., "Tri-Polarized 
% Holographic MIMO Surfaces for Near-Field Communications: 
% Channel Modeling and Precoding Design," in IEEE Transactions 
% on Wireless Communications, doi: 10.1109/TWC.2023.3266298.

% Please cite this paper properly if you use any part of this code
%%========================================================
clc
clear   
% fid=fopen('result.txt','a+');
%%========================================================
%  Set the parameters
% Environment setting
v_light=3*10^8;
waveFreq=1e3;  
% par.waveLambda=v_light/waveFreq;
par.waveLambda=1;
par.wavenumber=2*pi/par.waveLambda;
% system parameter setting
[transmit,receive,users,nearFieldRegion]=SystemConfig(par);

% The distance between each user and the transmitter
temp_UETX_dis.X=zeros(users.num*receive.totalNum,transmit.totalNum);
temp_UETX_dis.Y=zeros(users.num*receive.totalNum,transmit.totalNum);

for index_user=1:users.num 
    temp_UETX_dis.Z=users.start.z(index_user);
    temp_users_coordX=users.start.h(index_user):receive.spacing.h...
        :users.start.h(index_user)+(receive.num.h-1)*receive.spacing.h;
    users.coordinateX(index_user,:)=kron(ones(1,receive.num.v),...
        temp_users_coordX);
    temp_users_coordY=users.start.v(index_user):receive.spacing.v...
        :users.start.v(index_user)+(receive.num.v-1)*receive.spacing.v; 
    users.coordinateY(index_user,:)=kron(temp_users_coordY,...
        ones(1,receive.num.h));
    % The distance between each user and the transmitter
    temp_UETX_dis.X((index_user-1)*receive.totalNum+1:index_user...
        *receive.totalNum,:)=bsxfun(@minus,transmit.coordinateX,...
        users.coordinateX(index_user,:).');
    temp_UETX_dis.Y((index_user-1)*receive.totalNum+1:index_user...
        *receive.totalNum,:)=bsxfun(@minus,transmit.coordinateY,...
        users.coordinateY(index_user,:).');
    temp_UETX=sqrt(temp_UETX_dis.X((index_user-1)*receive.totalNum...
        +1:index_user*receive.totalNum,:).^2+temp_UETX_dis.Y((index_user-1)...
        *receive.totalNum+1:index_user*receive.totalNum,:).^2 ...
        +temp_UETX_dis.Z.^2);
    par.UETX.distance((index_user-1)*receive.totalNum+1:index_user...
        *receive.totalNum,:)=temp_UETX;
end

[ChanPol]=NF_ChannelGen(transmit,receive,par,temp_UETX_dis);
ChannelFullPolar=[ChanPol.XX,ChanPol.XY,ChanPol.XZ;...
    ChanPol.XY,ChanPol.YY,ChanPol.YZ;...
    ChanPol.XZ,ChanPol.YZ,ChanPol.ZZ];


%% Precoding: UE cluster precoding 
[ChannelUEPre]=Precoding_UECluster(transmit,receive,users,ChanPol);
UEPrecodedChan=blkdiag(ChannelUEPre.XX,ChannelUEPre.YY,ChannelUEPre.ZZ);


%% Precoding: Two-layer precoding
[ChannelBD]=Precoding_TwoLayer(transmit,receive,users,ChanPol);
ChannelBDMatrix=blkdiag(ChannelBD.XX,ChannelBD.YY,ChannelBD.ZZ);

%% Power allocation 
% PA1: Polarization selection first, then power allocation
noise=1;
SNR_value_range=-10:1:4;
iter=1;
for SNR_value=SNR_value_range
    % The Two-layer precoding +PA1 (polarization selection based PA)
    [RatePA1_sum_TL(iter),~]=PowerPolarSelection(SNR_value,noise,users,receive,ChannelBD); 
    % The UE +PA1
    [RatePA1_sum_UE(iter),~]=PowerPolarSelection(SNR_value,noise,users,receive,ChannelUEPre); 
    
    % PA2:Equal power allocation to all polarized channels, then perform
    % water-filling to each polarized channel
    % The Two-layer precoding +PA2
    [RatePA2_sum_TL(iter),~]=PowerAllEqual(SNR_value,noise,users,receive,ChannelBD);
    % The UE +PA2
    [RatePA2_sum_UE(iter),~]=PowerAllEqual(SNR_value,noise,users,receive,ChannelUEPre);
    
    % PA3: Unequal power allocation to the different polarization channel
    % The Two-layer precoding +PA3
    [RatePA3_sum_TL(iter),~]=PowerUnequalAllocation(SNR_value,noise,users,receive,ChannelBD);
    % The UE +PA3 is not available
    RatePA3_sum_UE=0;
    
%     % PA4:Equal power allocation to all polarized channels, then perform
%     % water-filling to each polarized channel
%     % The Two-layer precoding +PA4
%     [RatePA4_sum_TL(iter),~]=PowerEqualAllocation(SNR_value,noise,users,receive,ChannelBD);
%     % The UE +PA4 is not available
%     RatePA4_sum_UE=0;
    
    iter=iter+1;
end
% Conclusion: PA1<PA2<(PA3,PA4), PA3 is similar to PA4, i.e., there is
% little improvement in using water filling in polarized channel power
% matrix (compared with equal power allocation to different poallrizations)

figure;
hold on;   
plot(SNR_value_range,RatePA1_sum_TL/users.num,'ro-','LineWidth',2);
plot(SNR_value_range,RatePA2_sum_TL/users.num,'bx-','LineWidth',2);
plot(SNR_value_range,RatePA3_sum_TL/users.num,'ks-','LineWidth',2);
plot(SNR_value_range,RatePA1_sum_UE/users.num,'--','LineWidth',2,'Color',[1,0.38,0]);
plot(SNR_value_range,RatePA2_sum_UE/users.num,'--','LineWidth',2,'Color',[0.54,0.17,0.89]);
xlabel('SNR (dB)','Interpreter','latex');
ylabel('Spectral efficiency (bits/s/Hz)','Interpreter','latex');
grid on
legend('Two-layer precoding+PA1','Two-layer precoding+PA2',...
    'Two-layer precoding+PA3','UE precoding+PA1',...
    'UE precoding+PA2','Interpreter','Latex');
 
 

 
 
 
 