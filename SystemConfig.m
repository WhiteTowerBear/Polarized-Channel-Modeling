function [transmit,receive,users,nearFieldRegion]=SystemConfig(par)

% HMIMOS setting 
transmit.num.h=15;
transmit.num.v=15;
transmit.spacing.h=1/2*par.waveLambda; % horizontal spacing
transmit.spacing.v=1/2*par.waveLambda; % vertical spacing
transmit.len.h=transmit.num.h*transmit.spacing.h;
transmit.len.v=transmit.num.v*transmit.spacing.v;
transmit.aperture=sqrt(transmit.len.h^2+transmit.len.v^2);
transmit.patchArea=transmit.spacing.h*transmit.spacing.v;
transmit.totalNum=transmit.num.h*transmit.num.v;


receive.num.h=15;
receive.num.v=15;
receive.spacing.h=1/2*par.waveLambda;
receive.spacing.v=1/2*par.waveLambda;
receive.len.h=receive.num.h*receive.spacing.h;
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
users.start.h=0:receive.aperture:users.num*receive.aperture; 
users.start.v=0:receive.aperture:users.num*receive.aperture; 
% users.start.z=receive.aperture:1*receive.aperture:1*users.num*receive.aperture;
% users.start.z=0.2*par.waveLambda:0.2*par.waveLambda:0.2*par.waveLambda*users.num;
% users.start.z=0.5*nearFieldRegion;
users.start.z=[10]*par.waveLambda;

end