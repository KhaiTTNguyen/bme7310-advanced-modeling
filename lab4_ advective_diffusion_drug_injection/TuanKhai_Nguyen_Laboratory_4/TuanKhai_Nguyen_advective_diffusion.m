clear all

figure(1);

% SOR
clear;
h=0.05;
a = 50;
w=1.785;

x=[0:h:1];
y=[0:h:1];
A=zeros(length(y), length(x));

% ---- BC -------
A(:,1) = 1;
error=1;
itr=0;
pitr=0;
while (error > 1e-5 & itr < 10000)
     itr=itr+1;
     Aold=A;
     for i=1:length(y)
        for j=2:length(x)    
            if i==1 % bottom condition
                if j<length(x) % bottom condition
                    A(i,j)=w/4*( 2*Aold(i+1,j) + Aold(i,j+1) +A(i,j-1) )+(1-w)*Aold(i,j); 
                else % bottom right corner
                    A(i,j)=w/4*( 2*Aold(i+1,j) + 2*A(i,j-1) )+(1-w)*Aold(i,j); 
                end
            
            elseif 1<i & i<length(y) % middle
                if j < length(x) % in middle
                    A(i,j)=w/4*( Aold(i+1,j)+Aold(i,j+1)+A(i-1,j)+A(i,j-1) )+(1-w)*Aold(i,j);
                else % right condition
                    A(i,j)=w/4*(Aold(i+1,j)+A(i-1,j)+2*A(i,j-1))+(1-w)*Aold(i,j);
                end
            
            elseif i==length(y) % top condion segment by j
                if j<=11
                    A(i,j)=w/4*( 2*A(i-1,j) + Aold(i,j+1) + A(i,j-1) )+(1-w)*Aold(i,j); 
                elseif 11<j & j < 21
                    A(i,j)=w/4*( 2*A(i-1,j) + Aold(i,j+1) + A(i,j-1) + 2*a*h )+(1-w)*Aold(i,j);
                else % top right corner
                    A(i,j)=w/4*( 2*A(i-1,j) + 2*A(i,j-1) + 2*a*h)+(1-w)*Aold(i,j); 
                end   
            end 
        end
      end  
     error=max(max(abs(A-Aold)))/max(max(abs(A)));
    pitr=pitr+1;
    if pitr==5
       figure(1)
       contour(y,x,A,20); % y row, x col
       pitr=0;
    end   
end
 
fprintf('SOR Iterations %d\n',itr);
figure(1), clf
contour(y,x,A,20);
xlabel('x')
ylabel('y')
title(['SOR: # of Iterations ' num2str(itr)]);

fprintf('Value at point  x=0.7, y=0.7 %.12f \n', A(14,14))
%
Vx=zeros(length(y), length(x));
Vy=zeros(length(y), length(x));

% left
Vx(:,1) = 0;
Vy(:,1) = 0;

%right
Vx(:,end) = 0;
Vy(:,end) = 0; % calc below

% bottom 
Vx(1,:) = 0;  % calc below
Vy(1,:) = 0; 

% top
Vx(end,:) = 0; % calc below
Vy(end,1:11) = 0;
Vy(end,12:21) = -a;

% mid
Vx(:, 2:end-1)= - ( A(:, 3:end) - A(:, 1:end-2) )/ (2*h);
Vy(2:end-1, 2:end) = - ( A(3:end, 2:end) - A( 1:end-2, 2:end) )/ (2*h);

figure(2);
quiver(x, y, Vx, Vy)
xlim([0,1.05])
ylim([-0.05,1.05])
title('Velocity vector')
% ---------------- Problem 2 - Center difference -----------------
clearvars -except Vx Vy x y

%
h=0.05;
w=1.785;

x=[0:h:1];
y=[0:h:1];
A=zeros(length(y), length(x));
A(16,12)=1;
% ---- BC -------
error=1;
itr=0;
pitr=0; 
while (error > 1e-5 & itr < 10000)
     itr=itr+1;
     Aold=A;
     A(12,16)=1;
     for i=1:length(y) % row
        for j=1:length(x) % col 
            if i==1 % bottom condition
                if j==1 % bottom left
                    A(i,j)=w/8*( 4*Aold(i+1,j) + 4*Aold(i,j+1) ) + (1-w)*Aold(i,j); 
                elseif 1<j & j<length(x) % bottom condition
                    A(i,j)=w/8*( 4*Aold(i+1,j) + Aold(i,j+1)*(2-h*Vy(i,j)) + A(i,j-1)*(2+h*Vy(i,j)) ) + (1-w)*Aold(i,j); 
                elseif j==length(x) % bottom right corner
                    A(i,j)=w/8*( 4*Aold(i+1,j) + 4*A(i,j-1) ) + (1-w)*Aold(i,j); 
                end
            
            elseif 1<i & i<length(y) % middle
                if j==1 % left 
                    A(i,j)=w/8*( Aold(i+1,j)*(2-h*Vx(i,j)) + 4*Aold(i,j+1) + A(i-1,j)*(2+h*Vx(i,j)) ) + (1-w)*Aold(i,j); 
                elseif 1<j & j < length(x) % in middle
                    if i==16 & j==12 % set before loop but just set again here
                        A(i,j)=1;
                    else
                        A(i,j)=w/8*(    Aold(i+1,j)*(2-h*Vx(i,j)) + ...
                                        Aold(i,j+1)*(2-h*Vy(i,j)) + ...
                                        A(i-1,j) * (2+h*Vx(i,j)) + ...
                                        A(i,j-1)* (2+h*Vy(i,j))  ) + (1-w)*Aold(i,j);
                    end 
                elseif j==length(x)  % right condition
                    A(i,j)=w/8*( Aold(i+1,j) * (2-h*Vx(i,j)) + A(i-1,j) * (2+h*Vx(i,j)) + 4*A(i,j-1) ) + (1-w)*Aold(i,j);
                end
            
            elseif i==length(y) % top condion segment by j
                if j==1 % top left corner
                    A(i,j)=w/8*( 4*A(i-1,j) + 4*Aold(i,j+1) ) + (1-w)*Aold(i,j); 
                elseif 1< j & j<=11
                    A(i,j)=w/8*( 4*A(i-1,j) + Aold(i,j+1)*(2-h*Vy(i,j)) + A(i,j-1)*(2+h*Vy(i,j)) )+(1-w)*Aold(i,j); 
                elseif 11<j & j <= 21
                    A(i,j)=0;
                end   
            end 
        end
      end  
     error=max(max(abs(A-Aold))) / max(max(abs(A)));
     pitr=pitr+1;
    if pitr==5
       figure(3)
       contour(y,x,A,20);
       pitr=0;
    end   
end
 
fprintf('SOR Iterations %d\n',itr);
figure(3)
contour(y,x,A,20);
xlabel('x')
ylabel('y')
title(['SOR: # of Iterations ' num2str(itr)]);

%----------------------- Upstream------------------------
%
h=0.05;
w=1.785;

x=[0:h:1];
y=[0:h:1];
A=zeros(length(y), length(x));
A(16,12)=1;
% ---- BC -------
error=1;
itr=0;
pitr=0; 
while (error > 1e-5 & itr < 10000)
     itr=itr+1;
     Aold=A;
     A(12,16)=1;
     for i=1:length(y) % row
        for j=1:length(x) % col 
            if i==1 % bottom condition
                if j==1 % bottom left
                    A(i,j)=w/(4 -Vx(i,j)*h - Vy(i,j)*h ) *( (2-Vx(i,j)*h)*Aold(i+1,j) + (2-Vy(i,j)*h)*Aold(i,j+1) ) + (1-w)*Aold(i,j); 
                elseif 1<j & j<length(x) % bottom condition
                    A(i,j)=w/(4 -Vx(i,j)*h - Vy(i,j)*h ) *( (2-Vx(i,j)*h)*Aold(i+1,j) + Aold(i,j+1)*(1-h*Vy(i,j)) + A(i,j-1) ) + (1-w)*Aold(i,j); 
                elseif j==length(x) % bottom right corner
                    A(i,j)=w/(4 -Vx(i,j)*h - Vy(i,j)*h ) *( (2-Vx(i,j)*h)*Aold(i+1,j) + (2-Vy(i,j)*h)*A(i,j-1) ) + (1-w)*Aold(i,j); 
                end
            
            elseif 1<i & i<length(y) % middle
                if j==1 % left 
                    A(i,j)=w/(4 -Vx(i,j)*h - Vy(i,j)*h ) *( Aold(i+1,j)*(1-h*Vx(i,j)) + (2-Vy(i,j)*h)*Aold(i,j+1) + A(i-1,j) ) + (1-w)*Aold(i,j); 
                elseif 1<j & j < length(x) % in middle
                    if i==16 & j==12 % set before loop but just set again here
                        A(i,j)=1;
                    else
                        A(i,j)=w/(4 -Vx(i,j)*h - Vy(i,j)*h ) ...
                                *(  Aold(i+1,j)*(1-h*Vx(i,j)) + ...
                                    Aold(i,j+1)*(1-h*Vy(i,j)) + ...
                                    A(i-1,j) + ...
                                    A(i,j-1)  ) + (1-w)*Aold(i,j);
                    end 
                elseif j==length(x)  % right condition
                    A(i,j)=w/(4 -Vx(i,j)*h - Vy(i,j)*h ) *( Aold(i+1,j) * (1-h*Vx(i,j)) + A(i-1,j) + A(i,j-1)* (2-h*Vy(i,j)) ) + (1-w)*Aold(i,j);
                end
            
            elseif i==length(y) % top condion segment by j
                if j==1 % top left corner
                    A(i,j)=w/(4 -Vx(i,j)*h - Vy(i,j)*h ) *( A(i-1,j)*(2-Vx(i,j)*h) + Aold(i,j+1)*(2-Vy(i,j)*h) ) + (1-w)*Aold(i,j); 
                elseif 1< j & j<=11
                    A(i,j)=w/(4 -Vx(i,j)*h - Vy(i,j)*h ) *( A(i-1,j)*(2-Vx(i,j)*h) + Aold(i,j+1)*(1-h*Vy(i,j)) + A(i,j-1) )+(1-w)*Aold(i,j); 
                elseif 11<j & j <= 21
                    A(i,j)=0;
                end   
            end 
        end
      end  
     error=max(max(abs(A-Aold))) / max(max(abs(A)));
     pitr=pitr+1;
    if pitr==5
       figure(4)
       contour(y,x,A,20);
       pitr=0;
    end   
end
 
fprintf('SOR Iterations %d\n',itr);
figure(4)
contour(y,x,A,20);
xlabel('x')
ylabel('y')
title(['SOR: # of Iterations ' num2str(itr)]);

%
%----------------------- Downstream------------------------
h=0.05;
w=1.785;

x=[0:h:1];
y=[0:h:1];
A=zeros(length(y), length(x));
A(16,12)=1;
% ---- BC -------
error=1;
itr=0;
pitr=0; 
while (error > 1e-5 & itr < 10000)
     itr=itr+1;
     Aold=A;
     A(12,16)=1;
     for i=1:length(y) % row
        for j=1:length(x) % col 
            if i==1 % bottom condition
                if j==1 % bottom left
                    A(i,j)=w/(4 + Vx(i,j)*h + Vy(i,j)*h) *( (2+Vx(i,j)*h) *Aold(i+1,j) + (2+Vy(i,j)*h) *Aold(i,j+1) ) + (1-w)*Aold(i,j); 
                elseif 1<j & j<length(x) % bottom condition
                    A(i,j)=w/(4 + Vx(i,j)*h + Vy(i,j)*h) *( (2+Vx(i,j)*h) *Aold(i+1,j) + Aold(i,j+1) + A(i,j-1)*(1+h*Vy(i,j)) ) + (1-w)*Aold(i,j); 
                elseif j==length(x) % bottom right corner
                    A(i,j)=w/(4 + Vx(i,j)*h + Vy(i,j)*h) *( (2+Vx(i,j)*h) *Aold(i+1,j) + (2+Vy(i,j)*h) *A(i,j-1) ) + (1-w)*Aold(i,j); 
                end
            
            elseif 1<i & i<length(y) % middle
                if j==1 % left 
                    A(i,j)=w/(4 + Vx(i,j)*h + Vy(i,j)*h) *( Aold(i+1,j) + (2+Vy(i,j)*h) *Aold(i,j+1) + A(i-1,j)*(1+h*Vx(i,j)) ) + (1-w)*Aold(i,j); 
                elseif 1<j & j < length(x) % in middle
                    if i==16 & j==12 % set before loop but just set again here
                        A(i,j)=1;
                    else
                        A(i,j)=w/(4 + Vx(i,j)*h + Vy(i,j)*h) ...
                                *(    Aold(i+1,j) + ...
                                        Aold(i,j+1) + ...
                                        A(i-1,j) * (1+h*Vx(i,j)) + ...
                                        A(i,j-1)* (1+h*Vy(i,j))  ) + (1-w)*Aold(i,j);
                    end 
                elseif j==length(x)  % right condition
                    A(i,j)=w/(4 + Vx(i,j)*h + Vy(i,j)*h) *( Aold(i+1,j) + A(i-1,j) * (1+h*Vx(i,j)) + (2+h*Vy(i,j)) *A(i,j-1) ) + (1-w)*Aold(i,j);
                end
            
            elseif i==length(y) % top condion segment by j
                if j==1 % top left corner
                    A(i,j)=w/(4 + Vx(i,j)*h + Vy(i,j)*h) *(  (2+h*Vx(i,j))*A(i-1,j) + (2+h*Vy(i,j))*Aold(i,j+1) ) + (1-w)*Aold(i,j); 
                elseif 1< j & j<=11
                    A(i,j)=w/(4 + Vx(i,j)*h + Vy(i,j)*h) *( (2+h*Vx(i,j))*A(i-1,j) + Aold(i,j+1) + A(i,j-1)*(1+h*Vy(i,j)) )+(1-w)*Aold(i,j); 
                elseif 11<j & j <= 21
                    A(i,j)=0;
                end   
            end 
        end
      end  
     error=max(max(abs(A-Aold))) / max(max(abs(A)));
     pitr=pitr+1;
    if pitr==5
       figure(5)
       contour(y,x,A,20);
       pitr=0;
    end   
end
 
fprintf('SOR Iterations %d\n',itr);
figure(5)
contour(y,x,A,20);
xlabel('x')
ylabel('y')
title(['SOR: # of Iterations ' num2str(itr)]); 

%% ---------------- Diffusion - Problem 3-----------------
clear all
h=0.006;
w=1.785;
x=[0:h:0.3];
A= 2.5 ./ (sqrt(2*pi)) .* ( exp(-(30 .*x - 4.5).^2./2) );
figure(20), clf
hold all

% ---- BC -------
error=1;
itr=0;
time_points = [0, 0.2, 0.4, 0.6, 1.0, 1.8, 3, 4];
tolCheck = 10*eps;
dt=0.001;
D = 0.001;
r = D*dt/(h^2);
while (error > 1e-4 & itr < 10000)
     Aold=A;
     for i=1:length(x)
        if i==1
            A(i) = Aold(i)*(1-2*r) + 2*Aold(i+1)*r;
        elseif 1<i & i<length(x)
            A(i) = Aold(i)*(1-2*r) + Aold(i+1)*r + Aold(i-1)*r;
        elseif i==length(x)
            A(i) = Aold(i)*(1-2*r) + 2*Aold(i-1)*r;
        end
     end  
     error=max(max(abs(A-Aold)))/max(max(abs(A)));
    if any(abs(time_points - itr*dt) <= tolCheck)
       figure(20)    
       plot(x,A, 'DisplayName',['Time: '+  string(itr*dt)+' sec'])
    end   
    itr=itr+1;
end
 
fprintf('Explicit FDM Iterations %d\n',itr);
figure(20)
plot(x,A, 'DisplayName',['Time: '+  string(itr*dt)+' sec'])
legend(gca,'show')
hold off
xlabel('x')
ylabel('C')
title(['Explicit FDM: # of Iterations '+ string(itr)+ '. Total time: ' + string(itr*dt)+ ' sec. ' + 'Timestep: '+ string(dt)]);
