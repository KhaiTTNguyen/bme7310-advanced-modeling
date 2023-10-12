%% ---------------- Diffusion - Problem 3-----------------
clear all

r0 = 10; % cm
x0 = 50; % cm;
dx = 2; % cm
h=dx;
Lx = 100; %cm
c = 9.905; % cm/s
tau = 0.05; % sec
x=[0:dx:Lx];
dt = 0.005; %sec
theta = 0.5;

r = sqrt((x-x0).^2); 

A = 30 * cos(pi * r / (2*r0) ) ; 
A(r>r0) = 0;
A = [0,A,0]; % extend PBC
Aprev = A;
Anext_new = A;
Anext_old = A;
time_step = 0;
time_points = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
tolCheck = 10*eps;
k = c^2*dt^2 / (h^2);

figure(1), clf
plot(x,A(2:end-1),'DisplayName',['Time: '+  string(time_step*dt)+' sec'])
legend(gca,'show')
xlabel('x')
ylabel('C')
title(['Implicit FDM: # of Iterations '+ string(time_step)+ '. Total time: ' + string(time_step*dt)+ ' sec. ' + 'Timestep: '+ string(dt)]);

x0 = [A(2)];
xL = [A(end-1)];
t = [0];

while (time_step*dt < 10)
    time_step=time_step+1;
    if time_step > 1    
        Aprev = A;
        A = Anext_new;
    elseif time_step == 1
        Aprev = A;
        Anext_new = A;
    end
    Anext_old = Anext_new;
    itr=0;
    pitr=0;
    error=1;
    while (error > 1e-5 & itr < 10000)  
        i=2;
        Anext_new(i) = -1 / (1 + tau*dt/2 + k*theta + k*h/(c*dt) + k*theta*tau*h/(2*c)) * ...
            ( ...
            -k*theta* Anext_old(i+1) + ...
            (-2 + ( 2*k + h*k*tau/c)*(1-theta) )* A(i) + (-2*k*(1-theta))*A(i+1) + ...
            (1 - tau*dt/2 + k*theta - k*h/(c*dt) + k*theta*h*tau/(2*c))*Aprev(i) + (-k*theta)*Aprev(i+1) ...
            );
        
        for i=3:(length(A)-2)
            Anext_new(i) =  -1 / (1 + tau*dt/2 + k*theta) * ...
                (...
                (-k*theta/2)*Anext_old(i-1) + (-k*theta/2)*Anext_old(i+1) +...
                (-k*(1-theta))*A(i-1) + (-2 + 2*k*(1-theta))*A(i) + (-k*(1-theta))*A(i+1) +...
                (-k*theta/2)*Aprev(i-1) + (1 - tau*dt/2 + k*theta)*Aprev(i) + (-k*theta/2)*Aprev(i+1)...
                );
        end

        i=length(A)-1;
        Anext_new(i) = -1 / (1 + tau*dt/2 + k*theta + k*h/(c*dt) + k*theta*tau*h/(2*c)) * ...
            ( ...
            (-k*theta)* Anext_old(i-1) + ...
            (-2 + ( 2*k + h*k*tau/c)*(1-theta) )* A(i) + (-2*k*(1-theta))*A(i-1) + ...
            (1 - tau*dt/2 + k*theta - k*h/(c*dt) + k*theta*h*tau/(2*c))*Aprev(i) + (-k*theta)*Aprev(i-1) ...
            );    
          
         error=max(max(abs(Anext_new-Anext_old)))/max(max(Anext_old));
         itr=itr+1;
         Anext_old = Anext_new;
         pitr=pitr+1;
    end
    % compute shadow nodes
    i=2;
    Anext_new(i-1) = (-2*h/(c*theta*dt)) * ( ...
         (Anext_new(i) - Aprev(i)) ...
         - c*dt * ( theta/(2*h) * (Anext_new(i+1) + Aprev(i+1) - Aprev(i-1)) + (1-theta)/h*(A(i+1) - A(i-1)) ) ...
         + tau*dt *(theta*( Anext_new(i) + Aprev(i) )/2 + (1-theta)*A(i)) ...
         );
    
     i=length(A)-1;
     Anext_new(i+1) = (-2*h/(c*dt*theta)) * ( ...
         (Anext_new(i) - Aprev(i)) ...
         + c*dt * (theta/(2*h) * ( -Anext_new(i-1) + Aprev(i+1) - Aprev(i-1)) + ((1-theta)/h)*(A(i+1) - A(i-1)) ) ...
         + tau*dt *(theta*( Anext_new(i) + Aprev(i) )/2 + (1-theta)*A(i)) ...
         );
    
     t(end+1) = time_step*dt;
     x0(end+1) = Anext_new(2);
     xL(end+1) = Anext_new(end-1);
     
    % plot 
    if any(abs(time_points - time_step*dt) <= tolCheck)
           figure(time_step)    
           plot(x,Anext_new(2:end-1), 'DisplayName',['Time: '+  string(time_step*dt)+' sec'])
           legend(gca,'show')
           xlabel('x')
            ylabel('C')
            title(['Implicit FDM: # of Iterations '+ string(time_step)+ '. Total time: ' + string(time_step*dt)+ ' sec. ' + 'Timestep: '+ string(dt)]);

    end   
end


[PkAmp, PkTime] = findpeaks(x0); 
[~,idx] = sort(PkAmp,'descend');
max_x0 = PkAmp(idx(1)); %Amplitude of the peak
max_time_x0 = t(PkTime(idx(1))); %Time of the peak

[PkAmp, PkTime] = findpeaks(xL); 
[~,idx] = sort(PkAmp,'descend');
max_xL = PkAmp(idx(1)); %Amplitude of the peak
max_time_xL = t(PkTime(idx(1))); %Time of the peak


figure(1)
plot(t,x0)
xlabel('t')
ylabel('x0')
title(['x0 over time. Maximum peak at '+ string(max_x0) + ' at '+ string(max_time_x0) + ' sec']);


figure(2)
plot(t, xL)
xlabel('t')
ylabel('xL')
title(['xL over time. Maximum peak at '+ string(max_xL)+ ' at '+ string(max_time_xL) + ' sec']);
