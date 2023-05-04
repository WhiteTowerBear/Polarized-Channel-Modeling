%========================================================a
%  Author: Wei Li;
%  Date:   26/08/ 2022,  SUTD;
%  Version: V1.0 
%  Note: This tests the impact of imperfect CSI, using LS 
%        CE to recover signal X, i.e., X=(QHP)^{-1}Y, where
%        Q is allocated power, P is precoding, Y is received

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

IP_var=10e-5;
ChanPol_IP.XX=ChanPol.XX+IP_var; % Imperfect CSI
ChanPol_IP.YY=ChanPol.YY+IP_var; 
ChanPol_IP.ZZ=ChanPol.ZZ+IP_var; 
ChanPol_IP.XY=ChanPol.XY+IP_var;  
ChanPol_IP.YZ=ChanPol.YZ+IP_var; 
ChanPol_IP.XZ=ChanPol.XZ+IP_var; 
%% Precoding: UE cluster precoding 
[ChannelUEPre]=Precoding_UECluster(transmit,receive,users,ChanPol);
UEPrecodedChan=blkdiag(ChannelUEPre.XX,ChannelUEPre.YY,ChannelUEPre.ZZ);

[ChannelUEPre_IP]=Precoding_UECluster(transmit,receive,users,ChanPol_IP);

%% Precoding: Two-layer precoding
[ChannelBD]=Precoding_TwoLayer(transmit,receive,users,ChanPol);
ChannelBDMatrix=blkdiag(ChannelBD.XX,ChannelBD.YY,ChannelBD.ZZ);

[ChannelBD_IP]=Precoding_TwoLayer(transmit,receive,users,ChanPol_IP);
ChannelBDMatrix_IP=blkdiag(ChannelBD_IP.XX,ChannelBD_IP.YY,ChannelBD_IP.ZZ);
%% Power allocation 
% PA1: Polarization selection first, then power allocation
SNR_value_range=-5:5:30;
iter=1;
for SNR_value=SNR_value_range 
    var_noise=10^(-0.1*SNR_value); 
    %% received signal    
    [~,ChanLambda,~]=svd(ChannelBDMatrix*ChannelBDMatrix');
    ChanLambda=diag(ChanLambda); 

    trans_X=sqrt(1/2)*(randn(size(ChanLambda,2),1));
    noise=sqrt(var_noise)*(randn(size(ChanLambda,1),1));
    Rec_Y=ChanLambda*trans_X+noise;

    trans_X_est=pinv(ChanLambda)*Rec_Y;  

    BER_all(iter)=norm(trans_X- trans_X_est,'fro')^2/(norm(trans_X,'fro')^2);

    [~,ChanLambda_IP,~]=svd(ChannelBDMatrix_IP*ChannelBDMatrix_IP');
    if length(ChanLambda_IP)<length(ChanLambda)
        ChanLambda_IP=[ChanLambda_IP,zeros(length(ChanLambda)...
            -length(ChanLambda_IP))]; 
    elseif length(ChanLambda_IP)>length(ChanLambda)
        ChanLambda_IP(length(ChanLambda)+1:end,:)=[];
    end
    ChanLambda_IP=diag(ChanLambda_IP); 

    trans_X_est_IP=pinv(ChanLambda_IP)*Rec_Y;  

    BER_all_IP(iter)=norm(trans_X- trans_X_est_IP,'fro')^2/(norm(trans_X,'fro')^2);

    iter=iter+1;
end 

figure;
hold on;   
plot(SNR_value_range,BER_all,'ko-','LineWidth',2);
plot(SNR_value_range,BER_all_IP,'ro-','LineWidth',2);
xlabel('SNR (dB)','Interpreter','latex');
ylabel('Bit error rate','Interpreter','latex');
grid on
legend('Two-layer precoding with perfect CSI',...
    'Two-layer precoding with imperct CSI','Interpreter','Latex');
 
 

 
 
 
 