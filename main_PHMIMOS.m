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
par.waveLambda=20;
par.wavenumber=2*pi/par.waveLambda;
% system parameter setting
[transmit,receive,users,nearFieldRegion]=SystemConfig(par);

% The distance between each user and the transmitter
temp_UETX_dis.X=zeros(users.num*receive.totalNum,transmit.totalNum);
temp_UETX_dis.Y=zeros(users.num*receive.totalNum,transmit.totalNum);
temp_UETX_dis.Z=repmat(kron(users.start.z.',ones(receive.totalNum,1)),...
    1,transmit.totalNum);
for index_user=1:users.num
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
        +temp_UETX_dis.Z((index_user-1)*receive.totalNum...
        +1:index_user*receive.totalNum,:).^2);
    par.UETX.distance((index_user-1)*receive.totalNum+1:index_user...
        *receive.totalNum,:)=temp_UETX;
end

[ChanPol]=NF_ChannelGen(transmit,receive,par,temp_UETX_dis);
% Full polarized channel
ChannelFullPolar=[ChanPol.XX,ChanPol.XY,ChanPol.XZ;...
    ChanPol.XY,ChanPol.YY,ChanPol.YZ;...
    ChanPol.XZ,ChanPol.YZ,ChanPol.ZZ];

SNR=10;
var_noise=10^(-0.1*SNR);
noise=sqrt(var_noise/2)*(randn(users.num*receive.totalNum,1)...
    +1i*randn(users.num*receive.totalNum,1));

%% Precoding: UE cluster precoding 
[ChannelUEPre]=Precoding_UECluster(transmit,receive,users,ChanPol);
UEPrecodedChan=blkdiag(ChannelUEPre.XX,ChannelUEPre.YY,ChannelUEPre.ZZ);


%% Precoding: Two-layer precoding
[ChannelBD]=Precoding_TwoLayer(transmit,receive,users,ChanPol);
ChannelBDMatrix=blkdiag(ChannelBD.XX,ChannelBD.YY,ChannelBD.ZZ);

%% Power allocation 
% PA1: Polarization selection first, then power allocation
SNR_value=10;
noise=1;
% The Two-layer precoding +PA1
[RatePA1_sum_TL,RatePA1_AllUsers]=PowerPolarSelection(SNR_value,noise,users,receive,ChannelBD); 
% The UE +PA1
[RatePA1_sum_UE,~]=PowerPolarSelection(SNR_value,noise,users,receive,ChannelUEPre); 

% PA2:Equal power allocation to all polarized channels, then perform
% water-filling to each polarized channel
% The Two-layer precoding +PA2
[RatePA2_sum_TL,RatePA2_AllUsers]=PowerAllEqual(SNR_value,noise,users,receive,ChannelBD);
% The UE +PA2
[RatePA2_sum_UE,~]=PowerAllEqual(SNR_value,noise,users,receive,ChannelUEPre);


% PA3: Unequal power allocation to the different polarization channel
% The Two-layer precoding +PA3
[RatePA3_sum_TL,RatePA3_AllUsers]=PowerUnequalAllocation(SNR_value,noise,users,receive,ChannelBD);
% The UE +PA3 is not available
RatePA3_sum_UE=0;

% PA4:Equal power allocation to all polarized channels, then perform
% water-filling to each polarized channel
% The Two-layer precoding +PA4
[RatePA4_sum_TL,RatePA4_AllUsers]=PowerEqualAllocation(SNR_value,noise,users,receive,ChannelBD);
% The UE +PA4 is not available
RatePA4_sum_UE=0;

% Conclusion: PA1<PA2<(PA3,PA4), PA3 is similar to PA4, i.e., there is
% little improvement in using water filling in polarized channel power
% matrix (compared with equal power allocation to different poallrizations)


 
TwoLayerPA=[RatePA1_sum_TL,RatePA2_sum_TL,RatePA3_sum_TL,RatePA4_sum_TL]
UEPA=[RatePA1_sum_UE,RatePA2_sum_UE,RatePA3_sum_UE,RatePA4_sum_UE] 
 
% CorChan.X=ChannelBD.XX*ChannelBD.XX'; 
% DiversityX=trace(CorChan.X)^2/(trace(CorChan.X.^2))
% 
% CorChan.Y=ChannelBD.YY*ChannelBD.YY'; 
% DiversityY=trace(CorChan.Y)^2/(trace(CorChan.Y.^2))
% 
% CorChan.Z=ChannelBD.ZZ*ChannelBD.ZZ'; 
% DiversityZ=trace(CorChan.Z)^2/(trace(CorChan.Z.^2))
% 
% DiversityX+DiversityY+DiversityZ

CorChan.X=ChanPol.XX*ChanPol.XX'; 
DiversityX=trace(CorChan.X)^2/(trace(CorChan.X.^2))

CorChan.Y=ChanPol.YY*ChanPol.YY'; 
DiversityY=trace(CorChan.Y)^2/(trace(CorChan.Y.^2))

CorChan.Z=ChanPol.ZZ*ChanPol.ZZ'; 
DiversityZ=trace(CorChan.Z)^2/(trace(CorChan.Z.^2))

DiversityX+DiversityY+DiversityZ

 
%% Channel correlation computation 
CorDisTransmitX=repmat(transmit.coordinateX,1,transmit.totalNum)...
    -kron(transmit.coordinateX,ones(1,transmit.totalNum));
% CorDisTransmitX=reshape(CorDisTransmitX,transmit.totalNum,transmit.totalNum);
% CorDisTransmitX=CorDisTransmitX.';
CorDisTransmitY=repmat(transmit.coordinateY,1,transmit.totalNum)...
    -kron(transmit.coordinateY,ones(1,transmit.totalNum));
% CorDisTransmitY=reshape(CorDisTransmitY,transmit.totalNum,transmit.totalNum);
% CorDisTransmitY=CorDisTransmitY.';
CorDisTransmit=CorDisTransmitX.^2+CorDisTransmitY.^2; % Squared dis between two ant.

% corCom_fun=@(DisXY,DisZ)(8*DisZ^2-DisXY)./((4*DisZ^2+DisXY).^(5/2));
% diagElements=0;
% for index_user=1:users.num 
%     CorMatrixEachUE= corCom_fun(CorDisTransmit,users.start.z(index_user));
%     diagElements= 1/(4*users.start.z(index_user)^3);
%     CorMatrixEachUE(logical(eye(size(CorMatrixEachUE))))=...
%         diagElements*ones(transmit.totalNum,1);
%     CorMatrixEachUE=CorMatrixEachUE./diagElements;
%     delta=CorDisTransmit(1,:)/users.start.z(index_user);  
%      figure 
%     plot(delta(1,:), CorMatrixEachUE(1,:),'LineWidth',2,'Color','r')
%  
% end

CorDisTransmit=unique(sort(CorDisTransmit)); 
CorDisTransmit= (sort(CorDisTransmit)); 
delta=sqrt(CorDisTransmit)/users.start.z(1);

rest=(1-delta.^2/8)./(1+delta.^2/4).^(5/2);

figure 
plot(delta, rest,'LineWidth',2,'Color','r')
 







