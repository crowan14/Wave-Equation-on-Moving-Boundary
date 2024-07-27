clear; close all;

%FREE VIBRATION OF STRING ON MOVING DOMAIN

%maximum length of domain
L = 1;

%length of time simulation
T = 1; 

%number of global spatial basis functions
N = 5; 

%determines which number basis function is the initial displacement
which = 1; 

%number of spatial integration points
x_pts = 150; 

%number of time steps
t_pts = 5E4; 

%spatial integration and plotting grid
x_grid = linspace(0,L,x_pts); 
dx = x_grid(2)-x_grid(1);
x_grid = x_grid + dx/2;
x_grid = x_grid(1:end-1);

%time grid
t_grid = linspace(0,T,t_pts); 
dt = t_grid(2)-t_grid(1);

%compute spatial derivatives of shape functions required for K
x = sym('x');
t = sym('t');
i = sym('i');

sym_basis = w(x,i,L) * H(x,t,@a,L,T);
W = matlabFunction( sym_basis , 'vars', {x,t,i} );

diff_sym_basis = diff( sym_basis , x , 1 );
W_x = matlabFunction( diff_sym_basis , 'vars', {x,t,i} );

diff_sym_basis = diff( sym_basis , t , 1 );
W_t = matlabFunction( diff_sym_basis , 'vars', {x,t,i} );

%check static solution with particular basis (determined by the time)
tstar = T/5;
xstar = a(tstar,L,T);

%compute integration constants in exact solution
coeff = [ xstar , 1 ; L , 1  ] \ [ xstar^3/6 ; L^3/6 ];
exact = ( -x_grid.^3/6 + coeff(1)*x_grid + coeff(2) );
vec = x_grid > xstar;
exact = vec .* exact;

%stiffness matrix and force vector
K = zeros(N,N); F = zeros(N,1);
for i=1:N
    F(i) = dx * sum( b(x_grid,L) .* W(x_grid,tstar,i) );
    for j=1:N
        K(i,j) = dx * sum( W_x(x_grid,tstar,i) .* W_x(x_grid,tstar,j) );
    end
end

%solve system for degrees of freedom
u = K \ F;

%build solution
sol = zeros(x_pts-1,1);
for i=1:N
    sol = sol + u(i) * W(x_grid,tstar,i)';
end

% figure(1)
% plot( x_grid , sol )
% hold on
% plot( x_grid ,exact )
% legend('approximate','exact')
% xlabel('x')
% ylabel('u(x)')
% title('Comparison of Solutions')
% hold off
% drawnow

% %plot motion of heaviside function used to expand domain
% for i=1:t_pts
%     if mod(i,100) == 0
%         figure(2)
%         plot( x_grid , H(x_grid,t_grid(i),@a,L,T) )
%         xlabel('x')
%         ylabel('H(x,t)')
%         title('Expansion of Domain')
%         drawnow
%     end
% end

% %plot basis at particular time
% tstar = 0;
% for j=1:N
%     figure(3)
%     plot( x_grid , W(x_grid,tstar,j) )
%     xlabel('x')
%     ylabel('w_i(x)')
%     title('Basis Set at Specified Time')
%     hold on
%     drawnow
% end

%store time-dependent mass and stiffness matrices
M = zeros(t_pts,N,N);
K = zeros(t_pts,N,N);

%build mass and stiffness matrices at each point in time
for i=1:t_pts
    for j=1:N
        for k=1:N
            %use summetry
            if j >= k
                M(i,j,k) = dx * sum( W(x_grid,t_grid(i),j) .* W(x_grid,t_grid(i),k) );
                K(i,j,k) = dx * sum( W_x(x_grid,t_grid(i),j) .* W_x(x_grid,t_grid(i),k) );

                M(i,k,j) = M(i,j,k);
                K(i,k,j) = K(i,j,k);
            end
        end
    end
end

%initial condition
ui = zeros(N,1);
ui(which) = 1;

%initial strain energy
Uinit = 0.5 * ui' * squeeze( K(1,:,:) ) * ui;

% %plot initial condition
% figure(4)
% plot(x_grid,W(x_grid,0,which))
% xlabel('x')
% ylabel('u')
% title('Initial Condition')
% drawnow

%store degrees of freedom and implement displacement IC with no initial velocity
u = zeros(N,t_pts);
u(:,1) = ui;
u(:,2) = ui;

%store energy at each time
U = zeros(t_pts,1);
U(1) = Uinit;

%time integration
for i=2:t_pts

    %mass and stiffness matrices at current time
    Mi = squeeze( M(i,:,:) );
    Ki = squeeze( K(i,:,:) );
    
    %forward euler update of solution
    u(:,i+1) = -dt^2 * ( Mi \ Ki ) * u(:,i) + 2 * u(:,i) - u(:,i-1);

    %elastic strain energy
    strain_energy = 0.5 * u(:,i)' * Ki * u(:,i);

    %need to compute new matrices for kinetic energy
    T1 = zeros(N,N);
    T2 = zeros(N,N);
    for j=1:N
        for k=1:N
            T1(j,k) = dx * sum( W_t(x_grid,t_grid(i),j) .* W_t(x_grid,t_grid(i),k) );
            T2(j,k) = dx * sum( W_t(x_grid,t_grid(i),j) .* W(x_grid,t_grid(i),k) );
        end
    end

    du = ( u(:,i) - u(:,i-1) ) / dt;

    kinetic_energy = 0.5 * du' * Mi * du + 0.5 * u(:,i)' * T1 * u(:,i) +  u(:,i)' * T2 * du;
    
    %total energy
    U(i) = strain_energy + kinetic_energy;

end

%reconstruct the solution from time varying basis
sol = zeros(x_pts-1,t_pts);
for i=1:t_pts
    for j=1:N
        sol(:,i) = sol(:,i) + u(j,i) * W(x_grid,t_grid(i),j)';
    end
end

%plot the solution every S steps
S = 100;
for i=1:t_pts
    if mod(i,S)==0
        figure(5)
        plot(x_grid,sol(:,i))
        xlabel('x')
        ylabel('u(x)')
        title('Solution')
        ylim([ -max(max(abs(sol))) max(max(abs(sol))) ])
        grid on
        drawnow
    end
end

%plot total energy vs. time
figure(6)
plot(t_grid,U)
xlabel('time')
ylabel('total energy')
title('Checking energy conservation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time shift of position of heaviside function
function val = a(t,L,T)
    
    %fraction of domain exposed initially
    init = 0.5;
    
    %fraction of initial domain remaining by the end
    final = 0.5;

    %shift of heaviside position
    val = init * L * ( 1 - (1-final) * t/(T) );

end

%smooth heaviside function shifted as a function of time
function val = H(x,t,a,L,T)
    p = 20;
    pos = a(t,L,T);
    val = 0.5 * ( 1 + tanh( p * ( x - pos ) ) );
end

%fourier basis
function val = w(x,i,L)
    val = sin( i * pi * x / L );
end

%body force for static problem
function val = b(x,L)
    val = x/L;
end