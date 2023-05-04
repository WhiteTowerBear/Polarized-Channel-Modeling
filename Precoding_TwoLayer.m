function [ChannelBD]=Precoding_TwoLayer(transmit,receive,users,ChanPol)
% function: 
% The two-layer precoding scheme, i.e., the first layer removes 
% cross-polarization interference using GE (Gaussian elimination) method, 
% the second layer removes inter-user interference using BD method. 
% The first-layer precoding 

ChannInterference=[zeros(users.num*receive.totalNum,transmit.totalNum),...
    ChanPol.XY,ChanPol.XZ;...
    ChanPol.XY,zeros(users.num*receive.totalNum,transmit.totalNum),ChanPol.YZ;...
    ChanPol.XZ,ChanPol.YZ,zeros(users.num*receive.totalNum,transmit.totalNum)];
ChannelSignal=blkdiag(ChanPol.XX,ChanPol.YY,ChanPol.ZZ); 

temp_invXY=pinv(ChanPol.XY);
temp_invYZ=pinv(ChanPol.YZ);
FirstBarA=temp_invXY*ChanPol.XZ;
FirstBarB=temp_invXY*ChanPol.YZ;
FirstBarC=-temp_invYZ*ChanPol.XZ*temp_invXY*ChanPol.YZ...
    -temp_invXY*ChanPol.XZ;
 
FirstBarCzero=eye(transmit.totalNum)-FirstBarC'*pinv(FirstBarC...
    *FirstBarC')*FirstBarC;
   
FirstPrecode.Z=ChanPol.ZZ'*pinv(ChanPol.ZZ*FirstBarCzero*...
    ChanPol.ZZ')*ChanPol.ZZ*FirstBarCzero;
FirstPrecode.Y=-FirstBarA*FirstPrecode.Z;
FirstPrecode.X=-FirstBarB*FirstPrecode.Z; 

FirstPrecodeMatrix=[FirstPrecode.X;FirstPrecode.Y;FirstPrecode.Z];
  
%  Precoding: Gaussian elimination to the whole channel
% temp_interference=eye(3*transmit.totalNum)-ChannInterference'*pinv(ChannInterference...
%     *ChannInterference')*ChannInterference;
% FirstPrecodeMatrix=ChannelSignal'*pinv(ChannelSignal*temp_interference*...
%     ChannelSignal')*ChannelSignal*temp_interference;
% ChannelFirstPrecode=ChannelFullPolar*FirstPrecodeMatrix; 

temp_signal=ChannelSignal*FirstPrecodeMatrix; % channel after perfect precoding
temp_noise=ChannInterference*FirstPrecodeMatrix;

% The second-layer precoding
% The channel components with the first precoding 
ChannelGE.XX=ChanPol.XX*FirstPrecode.X;
ChannelGE.YY=ChanPol.YY*FirstPrecode.Y;
ChannelGE.ZZ=ChanPol.ZZ*FirstPrecode.Z;

[SecondPrecode.X]=PrecodingBD(ChannelGE.XX,receive,users);
[SecondPrecode.Y]=PrecodingBD(ChannelGE.YY,receive,users);
[SecondPrecode.Z]=PrecodingBD(ChannelGE.ZZ,receive,users);

% SecondPrecodeMatrix=[SecondPrecode.X;SecondPrecode.Y;SecondPrecode.Z];
% ChannelSecondPrecode=blkdiag(ChannelGE.XX,ChannelGE.YY,ChannelGE.ZZ)...
%     *SecondPrecodeMatrix;

% The channel components with the second precoding 
ChannelBD.XX=ChannelGE.XX*SecondPrecode.X; 
ChannelBD.YY=ChannelGE.YY*SecondPrecode.Y;
ChannelBD.ZZ=ChannelGE.ZZ*SecondPrecode.Z;

end