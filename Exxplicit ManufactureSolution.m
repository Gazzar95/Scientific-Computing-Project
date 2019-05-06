%% Scientific computing Project 
%Manyfactured Solution Explicit Method
%Omar El Gazzar
clc
clear

%Iterations where delta_x = delta_y
Nx = 15;
Ny = 15;
nt = 30;

%Paramiter 
D = 1;

%Domain
ax = 0;
ay = ax;
bx = pi;
by = bx;

%Steps
dx = (bx-ax)/(Nx-1);
dy = (by-ay)/(Ny-1);
dt = (min([dx,dy])^(2)/4);


x= 0:dx:(bx-ax);      
y= 0:dy:(by-ay);
t = 0:dt:dt*nt;
x = repmat(x,1,1,nt+1);%convert to 3D
y = repmat(y,1,1,nt+1);%Convert to 3D
t = repmat(t,1,1,nt+1);%Convert to 3D

% pre-location 
v = 0;
U = zeros(Nx,Ny,nt+1);
U_true = zeros(Nx,Ny,nt+1);


%Phi and psi Drechle Boundry condition 
U(1,:,:) = 0;
U(Nx,:,:) = 0;
U(:,1,:)=0;
U(:,1,:)=0;
% U(:,:,1)=sin(2.*x).*sin(2.*x);
% U(1,k,1:nt+1) = repmat(50,1,Ny,nt+1);
% U(Nx,k,1:nt+1) = repmat(50,1,Ny,nt+1);






for n=1:nt
    
    %Visualisation<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
    h = surf(x(:,:,1),y(:,:,1),U(:,:,n)','EdgeColor','none');
    colorbar
    shading interp
    %axis ([0 2*pi 0 2*pi -180 70])
    title({['time (\itt) = ',num2str(n*dt)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('{\leftarrow} Spatial co-ordinate (y)')
    zlabel('Solution to Diffution Equation \rightarrow')
    drawnow; 
    colormap('gray');
    grid
    refreshdata(h)
    %Visualisation<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
    
v = 1;
for k=2:Ny-1
for j=2:Nx-1
        %Discritised explicit method equation 
        U(j,k,n+1) = U(j,k,n) + (D*dt)*...
            ((1/dx^(2))*(U(j-1,k,n)-2*U(j,k,n)+U(j+1,k,n))+...
            (1/dy^(2))*(U(j,k-1,n)-2*U(j,k,n)+U(j,k+1,n)))+7*exp(-t(n))*sin(2*x(j))*sin(2*y(k));
            U_true(j,k,n+1) = exp(-t(n))*sin(2*x(j))*sin(2*y(k));
            E(v) = U(j,k,n+1) - U_true(j,k,n+1);
            v = v + 1;
end
end
end
% 

  e1 = (1/(Nx*Ny*nt))* sum(U(:,:,:) - U_true(:,:,:)); %Absolute Error
  e2 = sqrt((1/(Nx*Ny*nt))* sum(U(:,:,:) - U_true(:,:,:)).^2); %Mean Square Error
  einf = max(abs(sum(U(:,:,:) - U_true(:,:,:))));
  L5 = (1/(Nx*Ny*nt)).*abs(E);
  G = 1:length(L5);
  semilogx(G,log(L5))
  xlabel('log(Iteration Number)');
  ylabel('log(Error)');

