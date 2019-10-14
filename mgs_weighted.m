function [A,R] = mgs_weighted(A,M)
%MGS_WEIGHTED     M-orthogonal version of the QR factorization.
%
%  [Q,R] = MGS_WEIGHTED(A,M), where A is m-by-n and M is an m-by-m SPD matrix.
%  on output, Q is an m-by-m matrix satisfying Q'*M*Q = eye(m) and R is
%  an m-by-n upper-triangular matrix..
%
%  The factorization is performed by the modified Gram-Schmidt algorithm 
%  where the SPD matrix, M, is used to define inner products.  No test is 
%  performed to ensure that M satisfies these properties.  We also assume 
%  that m>=n.
%
%  For memory efficiency, the input matrix A is reused to store the columns 
%  of Q.
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

  [m,n] = size(A);
%   if (n>m)
%     error('A must be a square or tall matrix');
%   end

  if ( nargin<2 )
    M = speye(m);
  end

  % Preallocate storage
  R = zeros(n,n);
  c = zeros(1,n);  % index array for implicit column pivots
  y = zeros(1,n);  % array for storing column norms

  for k=1:n
    c(k) = k;
    y(k) = A(:,k)'*M*A(:,k);
  end

  for k=1:n
    %  get the index of the largest remaining column
    [tmp,idx] = max(y(k:n));
    p = idx+k-1;

    if abs( tmp(1) )<1e-30     % ~eps^2, could add relative error for small M...
      break   %  A does not have full column rank, stop here.
    
    else
      %  perform (implicit) column pivoting
      [c(p),c(k)] = deal(c(k),c(p));
      [y(p),y(k)] = deal(y(k),y(p));

      %  perform MGS update for the pivoted matrix
      R(c(k),c(k)) = sqrt(A(:,c(k))'*M*A(:,c(k)));
      A(:,c(k)) = A(:,c(k))/R(c(k),c(k));
      for l=k+1:n
        R(c(k),c(l)) = A(:,c(k))'*M*A(:,c(l));
        A(:,c(l)) = A(:,c(l)) - R(c(k),c(l))*A(:,c(k));
      end

      %  update the column norms
      for l=k+1:n
        y(l) = A(:,c(l))'*M*A(:,c(l));    % replace with update formula
      end
    end

  end
