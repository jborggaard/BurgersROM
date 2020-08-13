function [M] = compute_mass_matrix(x,e_conn,flag)
%%----------------------------------------------------------------------
%  compute_mass_matrix - computes a mass matrix \int(h_i(x) h_j(x)) dx
%                        from given nodal coordinates and element 
%                        connectivity
%                        
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    [M] = compute_mass_matrix( x , e_conn );
%
%  Variables:  x      
%                      nodal coordinates
%              e_conn
%                      element connectivity
%              flag
%                      'hermite' for 1d Hermite cubic elements
%                      'hermite_periodic' similar, but for periodic domains
%
%              M
%                      (sparse) mass matrix
%% ---------------------------------------------------------------------


  [N,          dim    ] = size(x     );
  [n_elements, nel_dof] = size(e_conn);

  % Due to memory limitations, we will fill the mass matrix by integrating
  % over partitions of elements.  The constant (max_elem_per_partition)
  % should be adjusted based on available memory.
  max_elem_per_partition = 10000;
  
  n_part = floor( n_elements / (max_elem_per_partition+1) ) + 1;
  elem_segment = floor( linspace(0,n_elements,n_part+1) );
  max_part_size = max( diff( elem_segment ) );

  if (dim==1)
    if (nargin==2)
      if (nel_dof==2)
        [r,w]   = oned_quadrature(2);
      elseif ( nel_dof==3 )
        [r,w]   = oned_quadrature(3);
      else
        [r,w]   = oned_quadrature(7);
      end
      
      for n_pt = 1:n_part
        II = zeros( max_part_size*nel_dof^2,1 );  
        JJ = zeros( max_part_size*nel_dof^2,1 );
        XX = zeros( max_part_size*nel_dof^2,1 );
        entry_counter = 0;

        for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
          nodes_local      = e_conn(n_el,:);
          x_local          = x(nodes_local,:);
          [~, w_g, phi, ~] = oned_shape(x_local,r,w);

          one = ones(size(w_g));
          M_loc = oned_bilinear( one, phi, phi, w_g );

          %  Assemble into the global system matrix
          for i=1:nel_dof
            for j=1:nel_dof
              entry_counter = entry_counter + 1;
              II(entry_counter) = nodes_local(i);
              JJ(entry_counter) = nodes_local(j);
              XX(entry_counter) = M_loc(i,j);
            end
          end
        end
    
        if ( n_pt==1 )
          M = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      N, N );
        else
          M = M + ...
              sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      N, N );
        end
      
      end % n_pt loop
    
      clear II JJ XX
    elseif (strcmp(flag(1:7),'hermite'))
      if ( nel_dof~=2 )
        error('COMPUTE_MASS_MATRIX: Only Hermite Cubic Elements Are Supported')
      end
      
      [r,w]   = oned_quadrature(4);
    
      ide = zeros(N,2);
      
      if ( strcmp(flag,'hermite_periodic') )
        for n_nd=1:N-1
           ide(n_nd,1) = 2*n_nd-1;
           ide(n_nd,2) = 2*n_nd;
        end  
        ide(N,1) = ide(1,1);
        ide(N,2) = ide(1,2);
  
        n_equations = 2*(N-1);  

      elseif ( strcmp(flag,'hermite_periodic_separate') )
        % the same matrix as above, with a different unknown numbering
        for n_nd=1:N-1
          ide(n_nd,1) = n_nd;
          ide(n_nd,2) = n_nd+N-1;
        end
        
        ide(N,1) = ide(1,1);
        ide(N,2) = ide(1,2);
  
        n_equations = 2*(N-1);  

      else
        for n_nd=1:N
          ide(n_nd,1) = 2*n_nd-1;
          ide(n_nd,2) = 2*n_nd;
        end
      end
     
      for n_pt = 1:n_part
        II = zeros( max_part_size*nel_dof^2,1 );  
        JJ = zeros( max_part_size*nel_dof^2,1 );
        XX = zeros( max_part_size*nel_dof^2,1 );
        entry_counter = 0;

        for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
          nodes_local          = e_conn(n_el,:);
          x_local              = x(nodes_local,:);
          [~, w_g, phi0, phi1] = oned_shapeherm(x_local,r,w);

          one = ones(size(w_g));
          M0_loc = oned_bilinear( one, phi0, phi0, w_g );
          M1_loc = oned_bilinear( one, phi1, phi1, w_g );
          S_loc  = oned_bilinear( one, phi0, phi1, w_g );

          %  Assemble into the global system matrix
          unk_0    = [ide(n_el,1); ide(n_el+1,1)];
          unk_1    = [ide(n_el,2); ide(n_el+1,2)];

          for i=1:nel_dof
            for j=1:nel_dof
              entry_counter = entry_counter + 1;
              II(entry_counter) = unk_0(i);
              JJ(entry_counter) = unk_0(j);
              XX(entry_counter) = M0_loc(i,j);
              
              entry_counter = entry_counter + 1;
              II(entry_counter) = unk_0(i);
              JJ(entry_counter) = unk_1(j);
              XX(entry_counter) = S_loc(j,i);
              
              entry_counter = entry_counter + 1;
              II(entry_counter) = unk_1(i);
              JJ(entry_counter) = unk_0(j);
              XX(entry_counter) = S_loc(i,j);
              
              entry_counter = entry_counter + 1;
              II(entry_counter) = unk_1(i);
              JJ(entry_counter) = unk_1(j);
              XX(entry_counter) = M1_loc(i,j);
            end
          end
        end
    
        if ( n_pt==1 )
          M = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      n_equations, n_equations );
        else
          M = M + ...
              sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      n_equations, n_equations );
        end
      
      end % n_pt loop
      
      
    else
      error('COMPUTE_MASS_MATRIX: One Dimensional Option Not Recognized')
    end
  end % dim==1
  
  if (dim==2)
    if ( nel_dof==3 || nel_dof==6 || nel_dof==7 )
      [r,s,w] = twod_quadrature(7);
      one = ones(size(w));
    
      for n_pt = 1:n_part
        II = zeros( max_part_size*nel_dof^2,1 );
        JJ = zeros( max_part_size*nel_dof^2,1 );
        XX = zeros( max_part_size*nel_dof^2,1 );
        entry_counter = 0;

        for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
          nodes_local               = e_conn(n_el,:);
          x_local                   = x(nodes_local,:);
          [x_g, w_g, phi, p_x, p_y] = twod_shape(x_local,r,s,w);

          M_loc = twod_bilinear( one, phi, phi, w_g );

          %  Assemble into the global system matrix
          for i=1:nel_dof
            for j=1:nel_dof
              entry_counter = entry_counter + 1;
              II(entry_counter) = nodes_local(i);
              JJ(entry_counter) = nodes_local(j);
              XX(entry_counter) = M_loc(i,j);
            end
          end
        end

        if ( n_pt==1 )
          M = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      N, N );
        else
          M = M + ...
              sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      N, N );
        end

      end % n_pt loop
    
    else
      [r,s,w] = twod_quadrature_Q(3);
      one = ones(size(w));
    
      for n_pt = 1:n_part
        II = zeros( max_part_size*nel_dof^2,1 );
        JJ = zeros( max_part_size*nel_dof^2,1 );
        XX = zeros( max_part_size*nel_dof^2,1 );
        entry_counter = 0;

        for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
          nodes_local               = e_conn(n_el,:);
          x_local                   = x(nodes_local,:);
          [ ~ , w_g, phi,  ~ ,  ~ ] = twod_shape_Q(x_local,r,s,w);

          M_loc = twod_bilinear( one, phi, phi, w_g );

          %  Assemble into the global system matrix
          for i=1:nel_dof
            for j=1:nel_dof
              entry_counter = entry_counter + 1;
              II(entry_counter) = nodes_local(i);
              JJ(entry_counter) = nodes_local(j);
              XX(entry_counter) = M_loc(i,j);
            end
          end
        end

        if ( n_pt==1 )
          M = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      N, N );
        else
          M = M + ...
              sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      N, N );
        end

      end % n_pt loop

    end
    
    clear II JJ XX
    
  end % dim==2

  if (dim==3)
    if ( nel_dof==4 || nel_dof==10 )
      % For tetrahedral meshes
      [r,s,t,w] = threed_quadrature(11);

      for n_pt = 1:n_part
        II = zeros( max_part_size*nel_dof^2,1 );  
        JJ = zeros( max_part_size*nel_dof^2,1 );
        XX = zeros( max_part_size*nel_dof^2,1 );
        entry_counter = 0;

        for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
          nodes_local                    = e_conn(n_el,:);
          x_local                        = x(nodes_local,:);
          [x_g, w_g, phi, p_x, p_y, p_z] = threed_shape(x_local,r,s,t,w);

          one = ones(size(w_g));
          M_loc = threed_bilinear( one, phi, phi, w_g );

          %  Assemble into the global system matrix
          for i=1:nel_dof
            for j=1:nel_dof
              entry_counter = entry_counter + 1;
              II(entry_counter) = nodes_local(i);
              JJ(entry_counter) = nodes_local(j);
              XX(entry_counter) = M_loc(i,j);
            end
          end
        end
    
        if ( n_pt==1 )
          M = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      N, N );
        else
          M = M + ...
              sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      N, N );
        end
      
      end % n_pt loop

    elseif ( nel_dof==8 )
      % For hexahedral meshes
      [r,s,t,w] = threed_hex_quadrature(3);
      
      for n_pt = 1:n_part
        II = zeros( max_part_size*nel_dof^2,1 );  
        JJ = zeros( max_part_size*nel_dof^2,1 );
        XX = zeros( max_part_size*nel_dof^2,1 );
        entry_counter = 0;

        for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
          nodes_local                    = e_conn(n_el,:);
          x_local                        = x(nodes_local,:);
          [x_g, w_g, phi, p_x, p_y, p_z] = threed_hex_shape(x_local,r,s,t,w);

          one = ones(size(w_g));
          M_loc = threed_bilinear( one, phi, phi, w_g );

          %  Assemble into the global system matrix
          for i=1:nel_dof
            for j=1:nel_dof
              entry_counter = entry_counter + 1;
              II(entry_counter) = nodes_local(i);
              JJ(entry_counter) = nodes_local(j);
              XX(entry_counter) = M_loc(i,j);
            end
          end
        end
    
        if ( n_pt==1 )
          M = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      N, N );
        else
          M = M + ...
              sparse( II(1:entry_counter), JJ(1:entry_counter),...
                      XX(1:entry_counter),...
                      N, N );
        end
      
      end % n_pt loop

    else
      error('element type unknown')
      
    end
    
    clear II JJ XX
    
  end % dim==3
  
  
end % function compute_mass_matrix
