function [ChannelUEPre]=Precoding_UECluster(transmit,receive,users,ChanPol)
% function: The precoding scheme based on user cluster, i.e., each
% polarized channel is assigned to different users.

ChannelUEPre.XX=ChanPol.XX;
ChannelUEPre.YY=ChanPol.YY;
ChannelUEPre.ZZ=ChanPol.ZZ;

user_num_x=zeros(users.num,1);
user_num_x(1:3:users.num,:)=1;
user_num_x=kron(user_num_x,ones(receive.totalNum,transmit.totalNum));
ChannelUEPre.XX=ChannelUEPre.XX.*user_num_x;

user_num_y=zeros(users.num,1);
user_num_y(2:3:users.num)=1;
user_num_y=kron(user_num_y,ones(receive.totalNum,transmit.totalNum));
ChannelUEPre.YY=ChannelUEPre.YY.*user_num_y;

user_num_z=zeros(users.num,1);
user_num_z(3:3:users.num)=1;
user_num_z=kron(user_num_z,ones(receive.totalNum,transmit.totalNum));
ChannelUEPre.ZZ=ChannelUEPre.ZZ.*user_num_z;

end