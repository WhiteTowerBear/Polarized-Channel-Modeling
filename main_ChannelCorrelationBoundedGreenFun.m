%========================================================a
%  Author: Wei Li;
%  Date:   22/10/ 2022,  SUTD;
%  Version: V1.0 
%  Note: Channel correlation of dyadic Green function due
%        to original source and image source

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


 

%% Free-space dydadic Green function (original source)
CorDisTransmitX=transmit.coordinateX;
CorDisTransmitY=transmit.coordinateY;
CorDisTransmit_original=CorDisTransmitX.^2+CorDisTransmitY.^2; % Squared dis between two ant.
CorDisTransmit_original=sqrt(CorDisTransmit_original);
prodTerm=par.wavenumber*CorDisTransmit_original;
ImGreen_var=3*CorDisTransmit_original.^2.*sin(prodTerm)./(par.wavenumber^2*CorDisTransmit_original.^5)...
    -3*CorDisTransmit_original.^2.*cos(prodTerm)./(par.wavenumber*CorDisTransmit_original.^4)...
    -CorDisTransmit_original.^2.*sin(prodTerm)./(CorDisTransmit_original.^3);
ImGreen_ave=sin(prodTerm)./CorDisTransmit_original...
    +cos(prodTerm)./(par.wavenumber*CorDisTransmit_original.^2)...
    -sin(prodTerm)./(par.wavenumber^2*CorDisTransmit_original.^3)+ImGreen_var/3;

%    ImGreen_ave=ImGreen_ave/ImGreen_ave(1,2);
figure
%  plot(1:transmit.totalNum, ImGreen_ave,'b:','LineWidth',1) 
hold on
 %% Dydadic Green function (image source) 
CorDisTransmitImage=CorDisTransmitX.^2+CorDisTransmitY.^2+users.start.z.^2; % Squared dis between two ant.
CorDisTransmitImage=sqrt(CorDisTransmitImage);
prodTerm=par.wavenumber*CorDisTransmitImage;
ImGreen_var_image=3*CorDisTransmitImage.^2.*sin(prodTerm)./(par.wavenumber^2*CorDisTransmitImage.^5)...
    -3*CorDisTransmitImage.^2.*cos(prodTerm)./(par.wavenumber*CorDisTransmitImage.^4)...
    -CorDisTransmitImage.^2.*sin(prodTerm)./(CorDisTransmitImage.^3);
ImGreen_ave_image=-sin(prodTerm)./CorDisTransmitImage...
    -cos(prodTerm)./(par.wavenumber*CorDisTransmitImage.^2)...
    +sin(prodTerm)./(par.wavenumber^2*CorDisTransmitImage.^3)+ImGreen_var_image/3;

%    ImGreen_ave_image=ImGreen_ave_image/ImGreen_ave_image(1,2);
%  plot(1:transmit.totalNum, ImGreen_ave_image,'k--','LineWidth',1) 

grid on 

ImGreen_ave_all=ImGreen_ave+ImGreen_ave_image;
% ImGreen_ave_all=ImGreen_ave_all/ImGreen_ave_all(1,2);
%  plot(1:transmit.totalNum, ImGreen_ave_all,'b-','LineWidth',1) 
% 
%  legend('Original source', 'Image source', 'The sum of two sources')
% 
% xlabel('Number of transmit patch antennas','Interpreter','latex');
% ylabel('Spatial correlation','Interpreter','latex'); 
% legend('$\check{z}=0.4\lambda$','Interpreter','latex')

% figure
%%%%% ImGreen_ave_all=ImGreen_ave_all/ImGreen_ave_all(1,2);
plot(1:transmit.totalNum, ImGreen_ave_all,'k','LineWidth',1) 
xlabel('Number of transmit patch antennas','Interpreter','latex');
ylabel('spatial correlation','Interpreter','latex'); 
legend('$\Delta^s_x=\Delta^s_y=0.01\lambda,\check{z}=0.4\lambda$','Interpreter','latex')
grid on











