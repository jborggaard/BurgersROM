function [w_save,x,t,e_conn] = burgers_1d_periodic(epsilon,q1,q2,Nx,Nt)
%-------------------------------------------------------------------------------
%  burgers_1d_periodic - Solves 1D viscous Burgers equation on a periodic domain
%
%    Solution is performed using a standard Galerkin finite element formulation.
%
%           w_t = - w w_x + epsilon w_{xx}, \Omega = (0,1)
%
%                       { q1 \sin(2 \pi x) + q2 \sin(4 \pi x), 0\leq x \leq .5 
%              w(x,0) = {
%                       { 0                                  , otherwise
%
%              w(0,t) = w(1,t)  (periodic boundary conditions)
%
%    the solution is integrated from t=0 to t=T_{max}=10.
%
%
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.0
%
%  Usage:    [w,x,t,e_conn] = burgers_1d_periodic(epsilon,q1,q2,Nx,Nt)
%
%  Variables:
%              epsilon
%                       Value of the "viscosity" parameter
%                       ( default = 1/100 )
%              q1,q2
%                       Parameters describing initial conditions
%                       ( default: q1 = 0.5 and q2 = 0.2 )
%              Nx
%                       Number of *elements* in (0,1), N = 2*Nx (quadratic elem)
%                       ( default = 50 )
%              Nt
%                       Number of uniform timesteps from t=0 to T_{max}
%                       ( default = 101 )
%
%  Outputs:
%              w
%                       Nodal solution matrix:  (N+1) rows x Nt columns
%              x
%                       Nodal coordinates
%              t
%                       Vector of time integration points
%
%  Example usage:
%              [w,x,t] = burgers_1d_periodic(.01,0.5,0.2,50,101);
%              surf(x,t,w')
%
%  Copyright 2019, Jeff Borggaard, Virginia Tech
%
%  Permission is hereby granted, free of charge, to any person obtaining 
%  a copy of this software and associated documentation files (the "Software"),
%  to deal in the Software without restriction, including without limitation 
%  the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%  and/or sell copies of the Software, and to permit persons to whom the 
%  Software is furnished to do so, subject to the following conditions:
%  
%  The above copyright notice and this permission notice shall be included 
%  in all copies or substantial portions of the Software.
%
%  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
%  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
%  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
%  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
%  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
%  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
%  DEALINGS IN THE SOFTWARE.
%
%-------------------------------------------------------------------------------

  % addpath('/Volumes/borggaard/Software/FEM/fem_functions')

%-------------------------------------------------------------------------------
%  Define "Input" Parameters
%-------------------------------------------------------------------------------
  if ( nargin==0 )
    epsilon = 0.01;
    q1      = 0.5;
    q2      = 0.2;
    Nx      = 50;
    Nt      = 101;
  end

  % parameters for nonlinear equation solve
  max_iterations   = 20;
  tol_step         = 1.d-8;
  tol_resid        = 1.d-8;

  % parameters for time integration
  theta            = 0.5;   % (explicit)  0 ... 1  (implicit)
  t                = linspace(0,10,Nt);

%-------------------------------------------------------------------------------
%  Geometry Module
%-------------------------------------------------------------------------------
  xb(:,1)      = [0; .5; 1];
  e_connb(1,:) = [1 2 3];
  rho(1)       = Nx;
  [x,e_conn,index_u,index_c] = oned_mesh(xb,e_connb,rho);


%-------------------------------------------------------------------------------
%  Solver Module
%-------------------------------------------------------------------------------
  [n_nodes   , n_dimensions] = size(x     );
  [n_elements, nel_dof     ] = size(e_conn);


  %-----------------------------------------------------------------------------
  % Determine equation numbers and set up boundary condition information
  %-----------------------------------------------------------------------------
  ide(1:n_nodes-1) = 1:n_nodes-1;   % setting equation numbers
  ide(n_nodes)     = 1;             % enforcing periodic bcs.

  %-----------------------------------------------------------------------------
  %  Set Initial Conditions
  %-----------------------------------------------------------------------------
  w0 = zeros(n_nodes,1);
  for n_nd=1:n_nodes
    if (x(n_nd)<=0.5)
      w0(n_nd) = q1*sin(2*pi*x(n_nd)) + q2*sin(4*pi*x(n_nd));
    else
      w0(n_nd) = 0;
    end
  end

  %-----------------------------------------------------------------------------
  %  Set storage for time history
  %-----------------------------------------------------------------------------
  w_save  = zeros( length(w0), Nt );

  %-----------------------------------------------------------------------------
  %  Perform Time Integration
  %-----------------------------------------------------------------------------
  n_store           = 1;
  w_save(:,n_store) = w0;
  n_store           = n_store + 1;

  wc = w0;
  t_step = t(2)-t(1);

  for k=2:Nt
    w  = wc;  % a good initial guess for the solution at the next time
              % step is the solution at the current time step.

    iteration = 0;
    converged = 0;
    diverged  = 0;

    while (~converged && ~diverged)
      iteration = iteration + 1;

      %-------------------------------------------------------------------------
      %  Compute the Newton step
      %-------------------------------------------------------------------------
      [b,A] = weak_resid(x,e_conn,w,wc,theta,t_step,epsilon,ide);

 %     b = b - (.75/Nx^2)*cos(10*t(k));
      
      newton_step =-A\b;

%     if ( iteration<3 )  % unsophisticated relaxation
%       lambda = 0.5;
%     else
        lambda = 1;
%     end

      %  Implicit update
      w(1:n_nodes-1) = w(1:n_nodes-1) + lambda*newton_step;
      w(n_nodes)     = w(1);

      norm_step  = norm(newton_step,2);
      norm_w     = norm(w,2);
      norm_resid = norm(b,2);

      if ( iteration==1 )
        norm_resid0 = max(norm_resid,1e-20);
      end

      converged = ( (norm_resid/norm_resid0 < tol_resid | ...
                     norm_resid < tol_resid           ) & ...
                     norm_step /norm_w      < tol_step );

      diverged = (iteration>max_iterations);
    end % newton iteration

    wc = w;
%   if ( mod(time-t_initial,t_save)==0 )
      w_save(:,n_store) = w;
      n_store = n_store + 1;
%   end

  end % time integration

end
%-------------------------------------------------------------------------------
%  Supporting Functions
%-------------------------------------------------------------------------------

%  oned_bilinear
%  oned_f_int
%  oned_gauss
%  oned_mesh
%  oned_shape
%
%  weak_resid

function [b,A] = weak_resid(x,e_conn,w,wc,theta,dt,epsilon,ide)
%-------------------------------------------------------------------------------
%  Function to compute the weak residual of the Burgers equation
%-------------------------------------------------------------------------------

  [n_elements,nel_dof] = size(e_conn);

  n_gauss = 3;
  [r,wt] = oned_quadrature(n_gauss);
  n_equations = size(x,1)-1;

  b  = zeros(n_equations, 1);

  II = zeros(3*n_elements,1);
  JJ = zeros(3*n_elements,1);
  AA = zeros(3*n_elements,1);
  n_triplets = 0;

  delta_w = 1e-5;

  for n_el=1:n_elements
    % compute value of each test function and spatial derivatives
    % at the integration points (x_g - Gauss points, wt_g - Gauss weights)
    nodes_local           = e_conn(n_el,:);
    x_local               = x(nodes_local,:);
    [x_g, wt_g, phi, p_x] = oned_shape(x_local,r,wt);

    % compute the values of w and w_current at the Gauss points
    wc_local = wc(nodes_local);
    w_local  = w (nodes_local);

    wc_g     = phi*wc_local;
    w_g      = phi*w_local;

    wc_xg    = p_x*wc_local;
    w_xg     = p_x*w_local;

    %---------------------------------------------------------------------------
    %  Integrate the weak form of the equations
    %---------------------------------------------------------------------------
    b_loc =               oned_f_int(        w_g, phi, wt_g)  ...
          + dt*theta    *(oned_f_int(  w_g.*w_xg, phi, wt_g)  ...
          +       epsilon*oned_f_int(       w_xg, p_x, wt_g)) ...
          -               oned_f_int(       wc_g, phi, wt_g)  ...
          + dt*(1-theta)*(oned_f_int(wc_g.*wc_xg, phi, wt_g)  ...
          +       epsilon*oned_f_int(      wc_xg, p_x, wt_g));

    %---------------------------------------------------------------------------
    % Assemble contributions into the global system matrix
    %---------------------------------------------------------------------------
    for n_t=1:nel_dof
      n_test = ide(nodes_local(n_t));
      if (n_test > 0)  % this is an unknown, fill the row
        b(n_test) = b(n_test) + b_loc(n_t);

        % calculate the column of Jacobian corresponding to n_t
        w_local(n_t) = w_local(n_t) + delta_w;

        w_g      = phi*w_local;
        w_xg     = p_x*w_local;

        b_plus   =        oned_f_int(        w_g, phi, wt_g) ...
          +     dt*theta*(oned_f_int(  w_g.*w_xg, phi, wt_g) ...
          +       epsilon*oned_f_int(       w_xg, p_x, wt_g))...
          -               oned_f_int(       wc_g, phi, wt_g) ...
          + dt*(1-theta)*(oned_f_int(wc_g.*wc_xg, phi, wt_g) ...
          +       epsilon*oned_f_int(      wc_xg, p_x, wt_g));

        w_local(n_t) = w_local(n_t) - delta_w;

        A_loc = (b_plus-b_loc)/delta_w;

        for n_u=1:nel_dof
          n_unk = ide(nodes_local(n_u));
          if (n_unk > 0)
            % A(n_unk,n_test) = A(n_unk,n_test) + A_loc(n_u);
            n_triplets = n_triplets + 1;
            II(n_triplets) = n_unk;
            JJ(n_triplets) = n_test;
            AA(n_triplets) = A_loc(n_u);
          end
        end
      end
    end
  end

  A = sparse(II(1:n_triplets),JJ(1:n_triplets),AA(1:n_triplets),...
             n_equations, n_equations);

end
