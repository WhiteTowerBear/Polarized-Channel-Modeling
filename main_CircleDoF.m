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

users.num=1; 
users.start.h=[0,0,0]; 
users.start.v=[0,0,0]; 
users.start.z=[5*par.waveLambda];

xaxix_range=0;
for index_user=1:users.num
    DOF=[];
    source_num_range=4:1:12; % One dimension 4:12
    xaxix_range=source_num_range;
    xaxix_range=transmit.len.h./xaxix_range;
    for source_num=source_num_range
         transmit.num.h=source_num;
        transmit.num.v=source_num;
        transmit.spacing.h=transmit.len.h/transmit.num.h; % horizontal spacing
        transmit.spacing.v=transmit.len.v/transmit.num.v; % vertical spacing 
        transmit.patchArea=transmit.spacing.h*transmit.spacing.v;
        transmit.totalNum=transmit.num.h*transmit.num.v;
        

        radius_range= sqrt(transmit.len.h*transmit.len.v/pi); 

         
        transmit.patchArea=transmit.spacing.h*transmit.spacing.v;
        transmit.totalNum=transmit.num.h*transmit.num.v;
        
         
        receive.num.h=source_num;
        receive.num.v=source_num;
        receive.spacing.h=2*receive.len.h/receive.num.h; % horizontal spacing
        receive.spacing.v=2*receive.len.v/receive.num.v; % vertical spacing 
        receive.patchArea=receive.spacing.h*receive.spacing.v;
        receive.totalNum=receive.num.h*receive.num.v;
        
        
        nearFieldRegion=2*(transmit.aperture+receive.aperture)^2/par.waveLambda; 
        
        transmit.start.h=-transmit.num.h/2;
        transmit.start.v=-transmit.num.v/2;
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
        
        users.coordinateX=transmit.coordinateX;
        users.coordinateY=transmit.coordinateY;
        % The distance between each user and the transmitter
        temp_UETX_dis.X=bsxfun(@minus,transmit.coordinateX,...
            users.coordinateX.');
        temp_UETX_dis.Y=bsxfun(@minus,transmit.coordinateY,...
            users.coordinateY.');
%         temp_dis_XY=sqrt(temp_UETX_dis.X.^2+temp_UETX_dis.Y.^2);
        
 
%         [row_nul,col_nul,~]=find((temp_UETX_dis.X-radius_range).^2 ...
%             +(temp_UETX_dis.Y-radius_range).^2>radius_range^2);
    
        temp_UETX_dis.X((transmit.coordinateX-radius_range).^2 ...
            +(transmit.coordinateY-radius_range).^2>radius_range^2)=0;
        temp_UETX_dis.Y((transmit.coordinateX-radius_range).^2 ...
            +(transmit.coordinateY-radius_range).^2>radius_range^2)=0;
        temp_dis_XY=sqrt(temp_UETX_dis.X.^2+temp_UETX_dis.Y.^2);
        
%         if length(temp_dis_XY)>square_numLimit
%             temp_UETX_dis.X(square_numLimit+1:end)=0;
%             temp_UETX_dis.Y(square_numLimit+1:end)=0;
%             temp_dis_XY(square_numLimit+1:end)=0; 
%         end

        temp_UETX=sqrt(temp_dis_XY.^2+temp_UETX_dis.Z.^2);
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
    plot(source_num_range.^2,DOF,'md-','LineWidth',1);
    xlabel('Number of patch antennas','Interpreter','latex');
    ylabel('DOF','Interpreter','latex');
    grid on
end

