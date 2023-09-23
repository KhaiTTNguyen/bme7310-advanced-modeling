% 1-D FDM exaample with Jacobi
figure(1);subplot(1,3,1);subplot(1,3,2);subplot(1,3,3);
Tol=1e-6;


% Solution 1 Efficient: Type 1 BC of V=1 @x=0, and V=0 @x=0.2

x=[0:.002:0.2];
xc=[0.001:0.002:0.2];
h=0.002;
V=zeros(1,101);
Vnew=V;
Sigma=ones(1,100);
Sigma(1,31:40)=4*Sigma(1,31:40);
Sigma(1,71:80)=0.1*Sigma(1,71:80);

V(1)=1;
V(101)=0;
Vnew(1)=1;
Vnew(101)=0;
Err=max(abs(Vnew));
itr=0;
Pitr=0;

while Err>Tol
  itr=itr+1;
  Pitr=Pitr+1;
  Vnew(2:100)= 1./(Sigma(1:99)+Sigma(2:100)).*(Sigma(1:99).*V(1:99)+Sigma(2:100).*V(3:101));
  Err=max(abs(V-Vnew));
  V(2:100)=Vnew(2:100);  
  if Pitr==100
      Pitr=0;
      figure(1);subplot(1,3,1);plot(x,V);xlabel('Position (m)');ylabel('Potential (V)');drawnow;
  end
end
title(['# of Iterations ' num2str(itr)]);

% Solution 2 Less Efficient: Type 1 BC of V=1 @x=0, and V=0 @x=0.2

x=[0:.002:0.2];
xc=[0.001:0.002:0.2];
h=0.002;
V=zeros(1,101);
Vnew=V;
Sigma=ones(1,100);
Sigma(1,31:40)=4*Sigma(1,31:40);
Sigma(1,71:80)=0.1*Sigma(1,71:80);

V(1)=1;
V(101)=0;
Vnew(1)=1;
Vnew(101)=0;
Err=max(abs(Vnew));
itr=0;
Pitr=0;

while Err>Tol
  itr=itr+1;
  Pitr=Pitr+1;
  for i=2:100
      Vnew(i)=1/(Sigma(i-1)+Sigma(i))*(Sigma(i-1)*V(i-1)+Sigma(i)*V(i+1));
  end
  Err=max(abs(V-Vnew));
  V(2:100)=Vnew(2:100);  
  if Pitr==100
      Pitr=0;
      figure(1);subplot(1,3,2);plot(x,V);xlabel('Position (m)');ylabel('Potential (V)');drawnow;
  end
end
title(['# of Iterations ' num2str(itr)]);

% Solution 3 Type 2 BC of J=# @x=0, and V=0 @x=0.2

% So the unused FDM equation would be
% Vshadow*(-sigma/h/h) + V1*(2*sigma/h/h) + V2*(-sigma/h/h)
% Taking the values from the previous solve, we can solve for Vshadow
V1=V(1);
V2=V(2);
Sig=1;
h=0.002;
Vshadow=(-h*h/Sig)*(V2*Sig/h/h -V1*2*Sig/h/h);
J=-Sig*(V2-Vshadow)/2/h;

x=[0:.002:0.2];
xc=[0.001:0.002:0.2];
h=0.002;
V=zeros(1,101);
Vnew=V;
Sigma=ones(1,100);
Sigma(1,31:40)=4*Sigma(1,31:40);
Sigma(1,71:80)=0.1*Sigma(1,71:80);

% Solution 3 Efficient: Type 2 BC of dV/dx=5 @x=0, and V=0 @x=0.2

V(101)=0;
Vnew(101)=0;
Err=max(abs(Vnew))+1;
itr=0;
Pitr=0;

while Err>Tol
  itr=itr+1;
  Pitr=Pitr+1;
  Vnew(1)=1./((Sigma(1)+Sigma(1))/h/h)*( ((2*Sigma(1).*V(2))/h/h) + J*2/h);
  Vnew(2:100)= 1./((Sigma(1:99)+Sigma(2:100))/h/h).*(Sigma(1:99).*V(1:99)/h/h+Sigma(2:100).*V(3:101)/h/h);
  Err=max(abs(V-Vnew));
  V(1:100)=Vnew(1:100);  
  if Pitr==100
      Pitr=0;
      figure(1);subplot(1,3,3);plot(x,V);xlabel('Position (m)');ylabel('Potential (V)');drawnow;
  end
end
title(['# of Iterations ' num2str(itr)]);

