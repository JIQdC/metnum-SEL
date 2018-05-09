g=zeros(5,3);
err=zeros(5,1);
for p=4:15

  n=m=p^2;
  m
  L=zeros(m,n);

  L=L+diag(diag(ones(n)));
  rand("seed",2);

  A=rand(m,n);
  
  A(1,1)=10E-10;

  xfijo=rand(n,1);
  b=A*xfijo;

  t1=time();
  x=LUPartPivot(n,A,b);
  t1=time()-t1
  
  errpart(p-3,1)=norm(xfijo-x);
  
  t2=time();
  x=LUSimple(n,A,b);
  t2=time()-t2
  
  errsimple(p-3,1)=norm(xfijo-x);
  
  g(p-3,1)=m;
  g(p-3,2)=t1;
  g(p-3,3)=t2;
  prom=0;

endfor
g
errpart
errsimple