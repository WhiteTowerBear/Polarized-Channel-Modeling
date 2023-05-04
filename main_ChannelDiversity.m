%========================================================a
%  Author: Wei Li;
%  Date:   06/09/ 2022,  SUTD;
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

transmit.len.h=15*par.waveLambda;
receive.len.h=15*par.waveLambda;
DoF_res=[]; 
transmitInitiaRange=4:9:36;
for transAddAnt=transmitInitiaRange
% system parameter setting 
par.UETX.distance=[];
transmit.num.h=transAddAnt;
transmit.num.v=1;
transmit.spacing.h=transmit.len.h/transmit.num.h; % horizontal spacing
transmit.spacing.v=0.5*par.waveLambda; % vertical spacing
 
transmit.len.v=transmit.num.v*transmit.spacing.v;
transmit.aperture=sqrt(transmit.len.h^2+transmit.len.v^2);
transmit.patchArea=transmit.spacing.h*transmit.spacing.v;
transmit.totalNum=transmit.num.h*transmit.num.v;


receive.num.h=4;
receive.num.v=1;
receive.spacing.h=receive.len.h/receive.num.h; 
receive.spacing.v=0.5*par.waveLambda; 
receive.len.v=receive.num.v*receive.spacing.v;
receive.aperture=sqrt(receive.len.h^2+receive.len.v^2);
receive.patchArea=receive.spacing.h*receive.spacing.v;
receive.totalNum=receive.num.h*receive.num.v;


nearFieldRegion=2*(transmit.aperture+receive.aperture)^2/par.waveLambda; 

transmit.start.h=0;
transmit.start.v=0;
transmit.coordinateX=transmit.start.h:transmit.spacing.h...
        :transmit.start.h+(transmit.num.h-1)*transmit.spacing.h;
transmit.coordinateX=kron(ones(1,transmit.num.v),transmit.coordinateX);
transmit.coordinateY=transmit.start.v:transmit.spacing.v...
        :transmit.start.v+(transmit.num.v-1)*transmit.spacing.v;
transmit.coordinateY=kron(transmit.coordinateY,ones(1,transmit.num.h));

users.num=1; 
users.start.h=receive.aperture:receive.aperture:users.num*receive.aperture;
users.start.h=users.start.h+3*par.waveLambda; % location of each user
users.start.v=receive.aperture:receive.aperture:users.num*receive.aperture;
users.start.v=users.start.v+3*par.waveLambda; 
% users.start.z=receive.aperture:1*receive.aperture:1*users.num*receive.aperture;
% users.start.z=0.2*par.waveLambda:0.2*par.waveLambda:0.2*par.waveLambda*users.num;
% users.start.z=0.5*nearFieldRegion;
users.start.z=0.5*par.waveLambda;



% The distance between each user and the transmitter
users.coordinateX=[];
users.coordinateY=[];
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
 

 
%% Diversity gain
chanDiagCoPolar=blkdiag(ChanPol.XX,ChanPol.YY,ChanPol.ZZ); 
corChan=ChannelFullPolar*ChannelFullPolar';
corChanX=ChanPol.XX'*ChanPol.XX; 
corChanY=ChanPol.YY'*ChanPol.YY;
corChanZ=ChanPol.ZZ'*ChanPol.ZZ;
DiversityX=(trace(corChanX)/norm(corChanX,'fro'))^2;
DiversityY=(trace(corChanY)/norm(corChanY,'fro'))^2;
DiversityZ=(trace(corChanZ)/norm(corChanZ,'fro'))^2;



DoF_res=[DoF_res,trace(corChan)^2/norm(corChan,'fro')^2 ];
end

figure
plot(transmitInitiaRange,DoF_res,'ro-','linewidth',1);
grid on
hold on


