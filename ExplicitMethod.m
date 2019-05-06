%% Scientific computing Project 

clc
clear
%Iterations where delta_x = delta_y
Nx = 10;
Ny = 10;
nt = 5;

%Paramiter 
D = 1;

%Domain
ax = 0;
ay = ax;
bx = 2*pi;
by = bx;

%Steps
dx = (bx-ax)/(Nx-1);
dy = (by-ay)/(Ny-1);
dt = (min([dx,dy])^(2)/4);


x= 0:dx:(bx-ax);      
y= 0:dy:(by-ay);
x = repmat(x,1,1,nt+1);%convert to 3D
y = repmat(y,1,1,nt+1);%Convert to 3D

%neumann boundry condition and pre-location 
v = 0;
U = zeros(Nx,Ny,nt+1);


%Phi and psi Drechle Boundry condition 
U(1,:,:) = (y-ay).^(2).*sin(pi.*(y-ay)./(2*(by-ay)));
U(Nx,:,:) = ((cos(pi.*(y-ay))-1).*cosh(by-y));
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
    refreshdata(h)
    %Visualisation<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
   
    
for k=2:Ny-1
for j=2:Nx-1
        %Discritised explicit method equation 
        U(j,k,n+1) = U(j,k,n) + (D*dt)*...
            ((1/dx^(2))*(U(j-1,k,n)-2*U(j,k,n)+U(j+1,k,n))+...
            (1/dy^(2))*(U(j,k-1,n)-2*U(j,k,n)+U(j,k+1,n)));
end
end

end

