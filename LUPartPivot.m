function x=LUPartPivot(n,A,b)
  m=n;
  
  L=zeros(m,n);

  L=L+diag(diag(ones(n)));
  rand("seed",2);
  
  y=zeros(n,1);
  x=zeros(n,1);

  for k=1:(n-1)
    max=k;
    for l=k:n
      if A(l,k)>A(max,k)
        max=l;
      endif
    endfor
    P=zeros(n,n);
    P=P+diag(diag(ones(n)));
    P(k,k)=P(max,max)=0;
    P(k,max)=P(max,k)=1;
    A=P*A;
    L=P*L;
    b=P*b;

    for i=(k+1):n
      L(i,k)=A(i,k)/A(k,k);
      for j=1:n
        A(i,j)=A(i,j)-L(i,k)*A(k,j);
      endfor
    endfor
  endfor

  for i=1:n
    y(i)=b(i);
    for j=1:(i-1)
      y(i)=y(i)-L(i,j)*y(j);
    endfor
  endfor

  for i=n:-1:1
    x(i)=y(i);
    for j=(i+1):n
      x(i)=x(i)-A(i,j)*x(j);
    endfor
    x(i)=x(i)/A(i,i);
  endfor
endfunction