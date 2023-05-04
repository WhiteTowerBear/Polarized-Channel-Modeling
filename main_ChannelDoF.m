%========================================================a
%  Author: Wei Li;
%  Date:   09/09/ 2022,  SUTD;
%  Version: V1.0 
%  Note: This is polarized HMIMOS channel modeling in near
%        field using Green's function, the DoF.
% Result: The fixed surface area has limit. 

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
% HMIMOS setting 

transmit.len.h=8*par.waveLambda;
transmit.len.v=8*par.waveLambda;
transmit.aperture=sqrt(transmit.len.h^2+transmit.len.v^2);
receive.len.h=8*par.waveLambda;
receive.len.v=8*par.waveLambda;
receive.aperture=sqrt(receive.len.h^2+receive.len.v^2);

users.num=3; 
users.start.h=[0,0,0]; 
users.start.v=[0,0,0]; 
users.start.z=[5*par.waveLambda,7*par.waveLambda,9*par.waveLambda];

for index_user=1:users.num
    DOF=[];
    source_num_range=[4:12]; % One dimension
    for source_num=source_num_range
        transmit.num.h=source_num;
        transmit.num.v=source_num;
        transmit.spacing.h=transmit.len.h/transmit.num.h; % horizontal spacing
        transmit.spacing.v=transmit.len.v/transmit.num.v; % vertical spacing 
        transmit.patchArea=transmit.spacing.h*transmit.spacing.v;
        transmit.totalNum=transmit.num.h*transmit.num.v;
        
         
        receive.num.h=source_num;
        receive.num.v=source_num;
        receive.spacing.h=receive.len.h/receive.num.h; % horizontal spacing
        receive.spacing.v=receive.len.v/receive.num.v; % vertical spacing 
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
         
        % The distance between each user and the transmitter
        temp_UETX_dis.X=[];
        temp_UETX_dis.Y=[];
        temp_UETX_dis.Z=users.start.z(index_user);
        
         
        users.coordinateX=[];
        users.coordinateY=[];
        par.UETX.distance=[];
        
        temp_users_coordX=users.start.h(index_user):receive.spacing.h...
            :users.start.h(index_user)+(receive.num.h-1)*receive.spacing.h;
        users.coordinateX=kron(ones(1,receive.num.v),...
            temp_users_coordX);
        temp_users_coordY=users.start.v(index_user):receive.spacing.v...
            :users.start.v(index_user)+(receive.num.v-1)*receive.spacing.v; 
        users.coordinateY=kron(temp_users_coordY,...
            ones(1,receive.num.h));
        % The distance between each user and the transmitter
        temp_UETX_dis.X=bsxfun(@minus,transmit.coordinateX,...
            users.coordinateX.');
        temp_UETX_dis.Y=bsxfun(@minus,transmit.coordinateY,...
            users.coordinateY.');
        temp_UETX=sqrt(temp_UETX_dis.X.^2+temp_UETX_dis.Y.^2 ...
            +temp_UETX_dis.Z.^2);
        par.UETX.distance =temp_UETX;
        
        
        [ChanPol]=NF_ChannelGen(transmit,receive,par,temp_UETX_dis);
        % Full polarized channel
        ChannelFullPolar=[ChanPol.XX,ChanPol.XY,ChanPol.XZ;...
            ChanPol.XY,ChanPol.YY,ChanPol.YZ;...
            ChanPol.XZ,ChanPol.YZ,ChanPol.ZZ];
         
        
        %% Theoretical analysis 
        % DoF
        corChan=ChannelFullPolar'*ChannelFullPolar; 
        DOF=[DOF,(trace(corChan)/norm(corChan,'fro'))^2];
    end
    figure;
    hold on;   
    plot(source_num_range.^2,DOF,'ks-','LineWidth',1);
    xlabel('Number of transmit antennas','Interpreter','latex');
    ylabel('DOF','Interpreter','latex');
    grid on
end

