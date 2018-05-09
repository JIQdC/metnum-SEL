err=zeros(20,1);
g=zeros(20,5);
c=zeros(8,3);

for n=1:10
  
  n
  N=8*n^3;
  A=Asparse(n);
  b = ones(N,1)/(n + 1)^2;
  tol=10^-6;
  x0=zeros(N,1);
  MAXIT=N;
  
  %Barra invertida%
  t3=0;
  for i=1:3
    t=time();
    xfijo=A\b;
    t3=t3+time()-t;
  endfor
  t3=t3/3;
 
  %PCG condicionado%
  t1=0;
  for i=1:3;
    t=time();
    [L U P]=ilu(A);
    x=pcg(A,b,tol,MAXIT,L,U,x0);
    t1=t1+time()-t;
  endfor
  t1=t1/3;
  
  
  %Descomposicion LU%
  t2=0;
  for i=1:3
    t=time();
    x=LUIntegrado(N,A,b);
    t2=t2+time()-t;
  endfor
  t2=t2/3;
  
  %PCG sin condicionar%
  t4=0;
  for i=1:3;
    t=time();
    x=pcg(A,b,tol,MAXIT,[],[],x0);
    t4=t4+time()-t;
  endfor
  t4=t4/3;
  
  %PCG sin condicionar%
  
  g(n,1)=N;
  g(n,2)=t1;
  g(n,3)=t2;
  g(n,4)=t3;
  g(n,5)=t4;
  if n>7
    continue
  endif
  
  %comparacion de numeros de condicion%
  sincond=eig(A);
  M=L*U;
  condic=eig(inverse(M)*A*inverse(M));
  kcond=max(condic)/min(condic);
  ksincond=max(sincond)/min(sincond);
  c(n,1)=N;
  c(n,2)=ksincond;
  c(n,3)=kcond;
  clear SD SU1 SU2 SU3 A b;
endfor
printf("\n\nN,Tpcgcond,Tlu,Tdirect,Tpcgsincond\n")
g
printf("\n\nN,kcond,ksincond\n")
c