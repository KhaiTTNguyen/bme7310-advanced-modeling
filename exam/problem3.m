clear all
% ------------- l=0.6cm ---------------
h=0.05;
l = 12; % 12 for 0.6 cm, 8 for 0.4
tumor_idx = l;
k = 0.15;

y=[0:h:1.2]; % y - row - i
x=[0:h:3]; % x - col - j
time_points = [10, 50, 100];
tolCheck = 10*eps;
dt=0.0005; % days

% ---- BC -------
Dxx = 4.95e-4;
Dyy = 5.8e-5;
Dgm = 6.90e-5;

rxx = Dxx*dt/(h^2);
ryy = Dyy*dt/(h^2);
rgm = Dgm*dt/(h^2);

A=zeros(length(y), length(x));
A(:,1) = 0;
A(:,end) = 0;
A(tumor_idx,30) = 0.01;

time_step=0;
time_v = [];
seed_6 = [];
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
                A(i,j) = rgm*(Aold(i+1,j) + Aold(i-1,j) + Aold(i,j+1) + Aold(i,j-1) - 4*Aold(i,j)) ...
                 + Aold(i,j) + dt*k*Aold(i,j) - dt*k*Aold(i,j)^2;

            elseif 18 < i & i <length(y) % second white matter
                 A(i,j) = ryy*Aold(i+1,j) + ryy*Aold(i-1,j) + rxx*Aold(i,j+1)+ rxx*Aold(i,j-1)...
                    + (-2*rxx - 2*ryy)*Aold(i,j) + k*dt*Aold(i,j) - k*dt*Aold(i,j)^2 + Aold(i,j) ;
            elseif i==length(y) % last row
                 A(i,j) = 2*ryy*Aold(i-1,j) + rxx*Aold(i,j+1) + rxx*Aold(i,j-1)...
                    + (-2*rxx - 2*ryy)*Aold(i,j) + k*dt*Aold(i,j) - k*dt*Aold(i,j)^2 + Aold(i,j) ;
            end 
         end
    
     end  
     % last column condition
     A(:,end)=0;
     time_v(end+1) = time_step*dt;
     seed_6(end+1) = A(tumor_idx,30);
     
     if any(abs(time_points - time_step*dt) <= tolCheck)
       figure(time_step+tumor_idx)    
       imagesc(x,y,A), colorbar
       axis equal
       title(['GM/WM domain with tumor at l = 0.6cm, '+ string(time_step*dt) + ' days']);
    end   
end

[a,b]=find(A>0.01);
% since length(a) = length(b)
tumor_area = 0.05^2 * length(a);
fprintf('Tumor area for l=%f: %f\n', l*0.05, tumor_area)
fprintf('Total area percentage: %f\n', 100* tumor_area/ (x(end)*y(end)))


%% ------------ l = 0.4 cm ---------------
h=0.05;
l = 8; % 12 for 0.6 cm, 8 for 0.4
tumor_idx = l;
k = 0.15;

y=[0:h:1.2]; % y - row - i
x=[0:h:3]; % x - col - j
time_points = [10, 50, 100];
tolCheck = 10*eps;
dt=0.0005; % days

% ---- BC -------
Dxx = 4.95e-4;
Dyy = 5.8e-5;
Dgm = 6.90e-5;

rxx = Dxx*dt/(h^2);
ryy = Dyy*dt/(h^2);
rgm = Dgm*dt/(h^2);

A=zeros(length(y), length(x));
A(:,1) = 0;
A(:,end) = 0;
A(tumor_idx,30) = 0.01;

time_step=0;
time_v = [];
seed_4 = [];
     
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
                A(i,j) = rgm*(Aold(i+1,j) + Aold(i-1,j) + Aold(i,j+1) + Aold(i,j-1) - 4*Aold(i,j)) ...
                 + Aold(i,j) + dt*k*Aold(i,j) - dt*k*Aold(i,j)^2;

            elseif 18 < i & i <length(y) % second white matter
                 A(i,j) = ryy*Aold(i+1,j) + ryy*Aold(i-1,j) + rxx*Aold(i,j+1)+ rxx*Aold(i,j-1)...
                    + (-2*rxx - 2*ryy)*Aold(i,j) + k*dt*Aold(i,j) - k*dt*Aold(i,j)^2 + Aold(i,j) ;
            elseif i==length(y) % last row
                 A(i,j) = 2*ryy*Aold(i-1,j) + rxx*Aold(i,j+1) + rxx*Aold(i,j-1)...
                    + (-2*rxx - 2*ryy)*Aold(i,j) + k*dt*Aold(i,j) - k*dt*Aold(i,j)^2 + Aold(i,j) ;
            end 
         end
    
     end  
     % last column condition
     A(:,end)=0;
     time_v(end+1) = time_step*dt;
     seed_4(end+1) = A(tumor_idx,30);
     if any(abs(time_points - time_step*dt) <= tolCheck)
       figure(time_step+tumor_idx)    
       imagesc(x,y,A), colorbar
       axis equal
       title(['GM/WM domain with tumor at l = 0.4cm, '+ string(time_step*dt) + ' days']);
    end   
end

[a,b]=find(A>0.01);
% since length(a) = length(b)
tumor_area = 0.05^2 * length(a);
fprintf('Tumor area for l=%f: %f\n', l*0.05, tumor_area)
fprintf('Total area percentage: %f\n', 100* tumor_area/ (x(end)*y(end)))

%% ----------------- l=0.2cm -----------------
h=0.05;
l = 4; % 12 for 0.6 cm, 8 for 0.4, 4 for 0.2
tumor_idx = l;
k = 0.15;

y=[0:h:1.2]; % y - row - i
x=[0:h:3]; % x - col - j
time_points = [10, 50, 100];
tolCheck = 10*eps;
dt=0.0005; % days

% ---- BC -------
Dxx = 4.95e-4;
Dyy = 5.8e-5;
Dgm = 6.90e-5;

rxx = Dxx*dt/(h^2);
ryy = Dyy*dt/(h^2);
rgm = Dgm*dt/(h^2);

A=zeros(length(y), length(x));
A(:,1) = 0;
A(:,end) = 0;
A(tumor_idx,30) = 0.01;

time_step=0;
time_v = [];
seed_2 = [];
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
                A(i,j) = rgm*(Aold(i+1,j) + Aold(i-1,j) + Aold(i,j+1) + Aold(i,j-1) - 4*Aold(i,j)) ...
                 + Aold(i,j) + dt*k*Aold(i,j) - dt*k*Aold(i,j)^2;

            elseif 18 < i & i <length(y) % second white matter
                 A(i,j) = ryy*Aold(i+1,j) + ryy*Aold(i-1,j) + rxx*Aold(i,j+1)+ rxx*Aold(i,j-1)...
                    + (-2*rxx - 2*ryy)*Aold(i,j) + k*dt*Aold(i,j) - k*dt*Aold(i,j)^2 + Aold(i,j) ;
            elseif i==length(y) % last row
                 A(i,j) = 2*ryy*Aold(i-1,j) + rxx*Aold(i,j+1) + rxx*Aold(i,j-1)...
                    + (-2*rxx - 2*ryy)*Aold(i,j) + k*dt*Aold(i,j) - k*dt*Aold(i,j)^2 + Aold(i,j) ;
            end 
         end
    
     end  
     % last column condition
     A(:,end)=0;
     time_v(end+1) = time_step*dt;
     seed_2(end+1) = A(tumor_idx,30);
     if any(abs(time_points - time_step*dt) <= tolCheck)
       figure(time_step+tumor_idx)    
       imagesc(x,y,A), colorbar
       axis equal
       title(['GM/WM domain with tumor at l = 0.2cm, '+ string(time_step*dt) + ' days']);
    end   
end

[a,b]=find(A>0.01);
% since length(a) = length(b)
tumor_area = 0.05^2 * length(a);
fprintf('Tumor area for l=%f: %f\n', l*0.05, tumor_area)
fprintf('Total area percentage: %f\n', 100* tumor_area/ (x(end)*y(end)))

%%
figure(10)
plot(time_v,seed_6,'DisplayName',['l = 0.6cm'])
hold on
plot(time_v,seed_4,'DisplayName',['l = 0.4cm'])
plot(time_v,seed_2,'DisplayName',['l = 0.2cm'])
hold off

legend(gca,'show')
xlabel('Time (days)')
ylabel('Cell Concentration at Seed Point')
title(['Cell Concentration at Seed Point vs Time']);

