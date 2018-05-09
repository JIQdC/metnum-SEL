g=zeros(5,2);
err=zeros(5,1);
for p=4:15
  
  m=n=p^2;
  A=rand(m,n);

  xfijo=rand(n,1);
  b=A*xfijo;
  
  t1=time();
  x=LUSimple(n,A,b);
  t2=time()-t1;
  
  g(p-1,1)=m;
  g(p-1,2)=t2;

  err(p-1,1)=norm(xfijo-x);
  
  [m t2 err(p-1,1)]
endfor
g
err