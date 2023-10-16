clear all
h=0.05;
a = 50;
l = 12; % 12 for 0.6 cm, 8 for 0.4
tumor_idx = l;
k = 0.15;

y=[0:h:1.2]; % y - row - i
x=[0:h:3]; % x - col - j

% ---- BC -------
time_points = [10, 50, 100];
tolCheck = 10*eps;
dt=0.0005; % days
Dxx = 4.95e-4;
Dyy = 5.8e-5;
Dgm = 6.90e-5;

rxx = Dxx*dt/(h^2);
ryy = Dyy*dt/(h^2);
rgm = Dgm*dt/(h^2);

A=zeros(length(y), length(x));
A(:,1) = 0;
A(:,end) = 0;
A(tumor_idx,30) = 0.1;

error=1;
pitr=0;
time_step=0;
while (time_step*dt < 100)
    time_step=time_step+1;

     Aold=A;
     % first column condition
     A(:,1)=0;
     % columns in between
     for j=2:length(x)-1 % x - col - j  C(i,j-1)=C(i,j+1)
        
         for i=1:length(y) % y - row - i 
            if i==1 % bottom row
                A(i,j) = 2*ryy*Aold(i+1,j) + rxx*Aold(i,j+1)+ rxx*Aold(i,j-1)...
                    + (-2*rxx - 2*ryy)*Aold(i,j) + k*dt*Aold(i,j) - k*dt*Aold(i,j)^2 + Aold(i,j) ;
            elseif 1 < i & i < 6 % first white matter
                 A(i,j) = ryy*Aold(i+1,j) + ryy*Aold(i-1,j) + rxx*Aold(i,j+1)+ rxx*Aold(i,j-1)...
                    + (-2*rxx - 2*ryy)*Aold(i,j) + k*dt*Aold(i,j) - k*dt*Aold(i,j)^2 + Aold(i,j) ;
            
            elseif 6 <= i & i <= 18 % grey matter
                % if i==tumor_idx & j==30
                %     A(i,j) = 0.1;
                % else
                A(i,j) = rgm*(Aold(i+1,j) + Aold(i-1,j) + Aold(i,j+1) + Aold(i,j-1) - 4*Aold(i,j)) ...
                 + Aold(i,j) + dt*k*Aold(i,j) - dt*k*A(i,j);
                % end
            elseif 18 < i & i <length(y) % second white matter
                 A(i,j) = ryy*Aold(i+1,j) + ryy*Aold(i-1,j) + rxx*Aold(i,j+1)+ rxx*Aold(i,j-1)...
                    + (-2*rxx - 2*ryy)*Aold(i,j) + k*dt*Aold(i,j) - k*dt*Aold(i,j)^2 + Aold(i,j) ;
            elseif i==length(y) % last row
                 A(i,j) = 2*ryy*Aold(i-1,j) + rxx*Aold(i,j-1) + rxx*Aold(i,j-1)...
                    + (-2*rxx - 2*ryy)*Aold(i,j) + k*dt*Aold(i,j) - k*dt*Aold(i,j)^2 + Aold(i,j) ;
            end 
         end
    
     end  
     % last column condition
     A(:,end)=0;
     if any(abs(time_points - time_step*dt) <= tolCheck)
       figure(time_step)    
       imagesc(A), colorbar
       axis equal
    end   
end
