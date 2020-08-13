 function [POD,w] = create_POD_basis(epsilon,q1,q2)
 %% ---------------------------------------------------------------------------
 %  CREATE_POD_BASIS
 %       see the README file
 %%----------------------------------------------------------------------------

  if ( nargin<3 )
  	epsilon = 0.01;
    q1      = 0.5;
   	q2      = 0.2;
  end

  [w,x,t,e_conn] = burgers_1d_periodic(epsilon,q1,q2,80,200);

  %  Comment the eye candy out for the experiments...
  figure(1)
  mesh(x,t,w')
  view([.7 -.8 .6])

  [M] = compute_mass_matrix(x,e_conn);

  [POD] = generate_pod(w,20,M);

  %  Same here
  figure(2)
  for i=1:6
    subplot(3,2,i)
    plot(x,POD(:,i))
    title(['POD Mode ',int2str(i)])
  end

end
