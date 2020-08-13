  function [POD,POD_eig,POD_energy] = generate_pod(snapshots,r,M,times)
%GENERATE_POD     Computes POD basis vectors from simulation snapshots.
%
%  Given snapshot data,
%  an (optional) maximum number of basis vectors, 
%  an (optional) weighting matrix M (usually the mass matrix in FEM), and
%  an (optional) set of times where snapshots are provided 
%  (the default assumes they are uniformly spaced), 
%
%  this function returns the solution to the discretized POD problem:
%
%        POD(:,i) = argmin_{POD(:,i).'*M*POD(:,j)=0, j<i}  E( POD(:,i) )
%
%  where E( POD(:,i) ) == \sum_{k=1}^n \| snapshots(k)-POD(:,i)*a(k) \|_M^2
%        a(k) = < snapshots(k), POD(:,i) >_M,
%        < POD(:,i), POD(:,i) >_M = 1, and < a, b >_M == a.'*M*b.
%
%  The only internal parameter the user should change is: "small_problem"
%  small_problem is a parameter that toggles the use of SVD or "the method of
%  snapshots" for computing approximations to the POD (Fredholm) eigenvalue
%  problem.
%
%
%  Copyright (c) 2019, Jeff Borggaard, Virginia Tech
%  Version: 3.0
%
%  Usage:   [POD,POD_eig,POD_energy] = GENERATE_POD(snapshots,r,M,times)
%
%  Variables:  POD          (output)
%                           the POD basis vectors
%
%              POD_eig      (output)
%                           eigenvalues corresponding to the Fredholm problem
%
%              POD_energy   (output)
%                           cumulative sum of the POD_energy, often used
%                           as a heuristic in choosing basis size
%
%              snapshots    (input)                                  size (NN,K)
%                           solution snapshots
%
%              r            (input, optional)
%                           number of POD basis vectors computed
%
%                           if unspecified, the default is K, the # of snapshots
%
%              M            (input, optional)              size (N,N) or (NN,NN)
%                           a positive definite, symmetric matrix used to
%                           define a weighted inner product,
%                           e.g. the finite element mass matrix
%
%                           if NN=N*c, and M is N*N, then the same M is used
%                           to weigh each of the c components.
%
%                           if unspecified, the default is speye(N)
%
%              times        (input, optional)
%                           vector of times at which snapshots are collected
%
%                           if unspecified, the default is to assume the
%                           snapshots are collected at uniformly spaced times
%
%  Example usage:
%              snapshots = [ u; v; w ];  % u is N x K  (K timesteps)
%
%  Note that either
%              MM                            = spalloc(3*N,3*N,3*nnz(M));
%              MM(    1:  N,    1:  N)       = M;
%              MM(  N+1:2*N,  N+1:2*N)       = M;
%              MM(2*N+1:3*N,2*N+1:3*N)       = M;
%              [POD,pod_energy,total_energy] = GENERATE_POD(snapshots,r,MM)
%  or
%              [POD,pod_energy,total_energy] = GENERATE_POD(snapshots,r,M)
%
%  produce the same output.
%
%  There is an internal parameter "small_problem" that toggles between
%  using the SVD and approximating solutions to the POD eigenvalue (Fredholm)
%  problem.
%%
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

  fprintf('Generating POD basis functions\n')
  small_problem = 5000;

  [n,q] = size(snapshots);

  if ( nargin<2 )
    warn = sprintf('r is not entered, attempting to compute %d POD vectors',q);
    warning(warn)
    r = q;
  end

  validateattributes(r,{'numeric'},{'size',[1,1]})
  validateattributes(r,{'numeric'},{'>',0})

  if ( nargin>2 )
    [nm,   ~] = size(M);

    if ( mod(n,nm)~=0 )
      error('snapshots and M have incompatible dimensions')
    end
  else
    M = speye(n);
  end

  if ( nargin<4 )     % equally weigh the snapshots
    Dt = speye(q);
  else                % weigh the snapshots consistent with the trapezoidal rule
    times = reshape(times,1,q);
    dd    = diff(times);
    Dt    = sparse(sqrt(diag([dd 0]/2 + [0 dd]/2)));
  end


  if ( n <= small_problem )   %  Solve "small problems" with the generalized SVD
    if ( nargin>=3 )
      %-------------------------------------------------------------------------
      %  Weight the snapshot data
      %-------------------------------------------------------------------------
      % get the dimension if same M is used for each direction...
      dim = n/nm;

      MM = sparse(n,n);
      for i=1:dim
        idx1 = 1+(i-1)*nm;
        idx2 = i*nm;
        MM(idx1:idx2,idx1:idx2) = M;
      end
      [Q,R] = mgs_weighted(snapshots,MM);

      clear M MM
      clear snapshots

      [U,Sigma,V] = svd(R,0);
      [n1,n2]  = size(U);
      POD = Q*U(:,1:r);
      clear U Q
      
    else

      [U,Sigma,V] = svd(snapshots,0);
      clear snapshots

      [n1,n2]  = size(U);

      POD = U(:,1:r);
      clear U

    end

    POD_eig   = diag(Sigma(1:min(n1,n2),1:min(n1,n2))).^2;
    tmp       = cumsum(POD_eig);
    POD_energy = tmp/tmp(end);

  else                                          % We have a "large size problem"
    dim  = n/nm;
    rows = 1:nm:n+1;

    % Compute the temporal autocorrelation matrix
    Rt = zeros(q,q);
    for i=1:dim
      Rt = Rt + snapshots(rows(i):rows(i+1)-1,:).'*M*...
                snapshots(rows(i):rows(i+1)-1,:);
    end

    Rt = 0.5*(Rt + Rt.');

    [Psi,D2] = eigs(Rt,r);
    [~,index] = sort(diag(D2),'descend');

    %sqrt(diag(D2(index,index)))

    POD = snapshots*Psi(:,index);

    for i=1:r
      normP2 = 0;
      for j=1:dim
        normP2 = normP2 + POD(rows(j):rows(j+1)-1,i).'*M*...
                          POD(rows(j):rows(j+1)-1,i);
      end
      POD(:,i) = POD(:,i) /sqrt( normP2 );
    end
    POD_eig   = diag( D2(index,index) );
    tmp       = cumsum(POD_eig);
    POD_energy = tmp/trace(Rt);
  end

end % function generate_pod
