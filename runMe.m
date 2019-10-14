%  A script to compute a POD model for Burgers equation.  This script drives
%  the following three steps:
%
%   (1)  compute the simulation data:  using burgers_1d_periodic.
%
%   (2)  use the simulation to compute a POD basis: using generate_pod.m
%
%   (3)  run test_model.m that builds the low-dimensional dynamical
%        system, simulates it, reconstructs an approximation to the solution
%        to Burgers equation, then computes the error in the approximation
%        (difference between FEM solution and POD solution).
%
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
%%

  % Create the simulation data at parameter "p1"
  %-----------------------------------------------------------------------------
  R  = 1000;
  q1 = 0.5;
  q2 = 0.0;
  r_dim  = 12;   % the number of POD basis functions to use...
    
  epsilon = 1/R;
  p1 = [ epsilon; q1; q2 ];
  
  % the 80 and 201 below are spatial and temporal FEM discretization parameters
  [z,x,t,e_conn] = burgers_1d_periodic(epsilon,q1,q2,80,201);
  
  %  Comment out the eye candy below for parametric experiments...
  %-----------------------------------------------------------------------------
  figure(1)
  mesh(x,t,z')
  view([.7 -.8 .6])
  pause(0.001)
  xlabel('x'); ylabel('t'); zlabel('z')
  title('Finite Element Solution')
    
  % Use the simulation data to compute a POD basis of dimension r_dim
  [M] = compute_mass_matrix(x,e_conn);
  
  [POD1,lam1] = generate_pod(z,r_dim,M);
        

  % Assemble the nonlinear ROM, perform simulation, and compute the error
  %-----------------------------------------------------------------------------
  [ROM_error] = test_model(POD1,r_dim,epsilon,q1,q2,M,z);

  fprintf('At P_1, p_1 = %g, p_2 = %g, and p_3 = %g,\n',epsilon,q1,q2);
  fprintf('  the ROM_error is %g\n',ROM_error);

  
