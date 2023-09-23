load LINE_LIST_1.DAT
in=LINE_LIST_1(:,2:3);
load POINT_LIST_1.DAT;
x=POINT_LIST_1(:,2);
y=POINT_LIST_1(:,3);
Jx=POINT_LIST_1(:,4);
Jy=POINT_LIST_1(:,5);

nn=length(x);

elist=zeros(nn,2);
for i=1:nn
    idx1=find(i==in(:,2));
    idx3=find(i==in(:,1));
    elist(i,1)=in(idx1,1);
    elist(i,2)=in(idx3,2);
end

DX=x(elist(:,2))-x(elist(:,1));
DY=y(elist(:,2))-y(elist(:,1));
S=sqrt(DX.*DX+DY.*DY);
NX=[DY./S];
NY=[-DX./S];

Jn=NX.*Jx+NY.*Jy;

S_line=sqrt((x(in(:,1))-x(in(:,2))).^2+(y(in(:,1))-y(in(:,2))).^2);
Jn1_plus_Jn2=Jn(in(:,1))+Jn(in(:,2));
Line_Integral=sum(0.5*S_line.*(Jn1_plus_Jn2));
fprintf('Your integral value for Part 1 is %f\n',Line_Integral);

load LINE_LIST_2.DAT
in=LINE_LIST_2(:,2:3);
load POINT_LIST_2.DAT;
x=POINT_LIST_2(:,2);
y=POINT_LIST_2(:,3);
Jx=POINT_LIST_2(:,4);
Jy=POINT_LIST_2(:,5);

nn=length(x);

elist=zeros(nn,2);
for i=1:nn
    idx1=find(i==in(:,2));
    idx3=find(i==in(:,1));
    elist(i,1)=in(idx1,1);
    elist(i,2)=in(idx3,2);
end

DX=x(elist(:,2))-x(elist(:,1));
DY=y(elist(:,2))-y(elist(:,1));
S=sqrt(DX.*DX+DY.*DY);
NX=[DY./S];
NY=[-DX./S];

Jn=NX.*Jx+NY.*Jy;

S_line=sqrt((x(in(:,1))-x(in(:,2))).^2+(y(in(:,1))-y(in(:,2))).^2);
Jn1_plus_Jn2=Jn(in(:,1))+Jn(in(:,2));
Line_Integral=sum(0.5*S_line.*(Jn1_plus_Jn2));
fprintf('Your integral value for Part II is %f\n',Line_Integral);
