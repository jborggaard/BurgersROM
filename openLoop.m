%  Compute Burgers equation matrices and perform open loop simulation

  n = 41;
  m = 1;   % a control input that we won't use
  
  epsilon = 0.005;
  tInfinity = 15;
  
  [E,A,B,N,zInit] = BurgersFEMControl(n,m);
  xNodes = linspace(0,1,n);

%   A = epsilon*(E\A);
%   B = E\B;
%   N = E\N;
  A = epsilon*A;
  options = odeset('Mass',E);

  zdot = @(t,z) A*z + N*kron(z,z);
  [T,Z] = ode23(zdot,[0 tInfinity],zInit,options);
  figure(10)
  mesh(xNodes,T,Z)
  xlabel('x'); ylabel('time')
  title('Open Loop Simulation')

%  Alternate by change of variable (to remove M matrix)
  sqM = sqrtm(full(M));
  sqMinv = inv(sqM);
  
  Ac = sqMinv*A*sqMinv;                  %#ok
  Bc = sqMinv*B;                         %#ok
  Nc = kroneckerRight(sqMinv*N,sqMinv);  %#ok
  Qc = eye(n);  R = r0*eye(m);
  
  rhs_open = @(t,x) [ Aopen*x(1:end-1) + Nopen*kron(x(1:end-1),x(1:end-1));...
                        x(1:end-1).'*Qc*x(1:end-1) ];
                        
  [t,x] = ode15s( rhs_open, [0 T], [x0;0], options );  
  figure(1); hold on
  Z = x(:,1:end-1)*sqMinv; Z(:,end+1) = Z(:,1);  %#ok
  mesh(xNodes,t,Z)
  xlabel('x'); ylabel('time')
  title('Open Loop Simulation')
  view([1 -1 1])
  axis([0 1 -0.1 5 -.4 .6]) 
  
  fprintf('Open Loop Cost (0,T) is %g\n\n',x(end,end));
