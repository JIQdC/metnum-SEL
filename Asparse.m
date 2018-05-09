function A = Asparse(n)
  N=8*n^3;
  SD=sparse(1:N,1:N,6*ones(N,1),N,N);
  if1=setdiff(1:N-1, n:n:n*(8*n^2-1));
  SU1=sparse(if1,if1+1,-1*ones(1,(n-1)*8*n^2),N,N);
  [i, j]=find([ones(n*(2*n-1),4*n); zeros(n,4*n)]);
  if2=i+(j-1)*2*n^2;
  SU2=sparse(if2,if2+n,-1*ones(1,n*(2*n-1)*4*n),N,N);
  if3=1:2*n^2*(4*n-1);
  SU3=sparse(if3,if3+2*n^2,-1*ones(1,2*n^2*(4*n-1)),N,N);
  A=SD+SU1+SU1'+SU2+SU2'+SU3+SU3'; % Si se quiere la matriz en formato RALO
endfunction