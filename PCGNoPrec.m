err=zeros(20,1);
g=zeros(20,2);

for n=1:20
  N=8*n^3;
  %A=Afull(n);
  A=Asparse(n);
  b = ones(N,1)/(n + 1)^2;
  tol=10^-6;
  x0=zeros(N,1);
  MAXIT=N;
  
  
  t1=time();
  x=pcg(A,b,tol,MAXIT,[],[],x0);
  t1=time()-t1();
  
  xfijo=A\b;
  err(n,1)=norm(xfijo-x)
  
  g(n,1)=N;
  g(n,2)=t1
endfor

g
err