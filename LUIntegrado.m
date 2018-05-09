function x=LUIntegrado(n,A,b)
  x=zeros(n,1);
  y=zeros(n,1);
  [L U]=lu(A);
  y=L\b;
  x=U\y;
endfunction
