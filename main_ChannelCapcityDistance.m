%========================================================a
%  Author: Wei Li;
%  Date:   25/10/ 2022,  SUTD;
%  Version: V1.0 
%  Note: This is polarized HMIMOS channel modeling in near
%        field using Green's function, the precoding design
%        and power allocation of the TP. 

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


 
%% Power allocation 
% PA1: Polarization selection first, then power allocation
noise=1; 
distance_range=0.5:1:15.5;
distance_range=distance_range*par.waveLambda; 
iter=1;
SNR_value=10;
for distance_value=distance_range
    users.start.z=distance_value;
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

[ChanPol]=DP_ChannelGen(transmit,receive,par,temp_UETX_dis);
ChannelFullPolar=[ChanPol.XX,ChanPol.XY,ChanPol.XZ;...
    ChanPol.XY,ChanPol.YY,ChanPol.YZ;...
    ChanPol.XZ,ChanPol.YZ,ChanPol.ZZ];
ChannelDP=[ChanPol.XX,ChanPol.XY;...
    ChanPol.XY,ChanPol.YY;...
    ChanPol.XZ,ChanPol.YZ];
ChannelSP=ChanPol.XX;

    trans_power=noise*10^(0.1*SNR_value);
    [~,lambda_FP,~]=svd(ChannelFullPolar*ChannelFullPolar');
    lambda_FP(lambda_FP==0)=[];
    SumRate_FP(1,iter)=sum(log2(1+trans_power*lambda_FP));
    [~,lambda_DP,~]=svd(ChannelDP*ChannelDP');
    lambda_DP(lambda_DP==0)=[];
    SumRate_DP(1,iter)=sum(log2(1+trans_power*lambda_DP));

    [~,lambda_SP,~]=svd(ChannelSP*ChannelSP');
    lambda_SP(lambda_SP==0)=[];
    SumRate_SP(1,iter)=sum(log2(1+trans_power*lambda_SP));
    
    iter=iter+1;
end
 

figure;
hold on;   
plot(distance_range,SumRate_FP/users.num,'ro-','LineWidth',2);
plot(distance_range,SumRate_DP/users.num,'bx-','LineWidth',2);
plot(distance_range,SumRate_SP/users.num,'ks-','LineWidth',2);
xlabel('distance $z (\mathrm{m})$ ','Interpreter','latex');
ylabel('Channel capacity (bits/s/Hz)','Interpreter','latex');
grid on
legend('TP HMIMOS','DP HMIMOS', 'Conventional HMIMOS',...
    'Interpreter','Latex');
 
 

 
 
 
 