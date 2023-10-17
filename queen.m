function f = queen(n,k) %generate spatial weight matrix according to Queen
m = n/k;                %choose n and k such that n is divisible by k.
permn = randperm(n);
R = reshape(permn, k, m); 
W = zeros(n,n);

%for inner units: each has eight neighbors
for i=2:(k-1); 
    for j=2:(m-1);
        id1 = R(i,j);
        id2 = R(i-1,j);
        id3 = R(i+1,j);
        id4 = R(i,j-1);
        id5 = R(i,j+1);
        id6 = R(i-1,j-1);
        id7 = R(i-1,j+1);
        id8 = R(i+1,j-1);
        id9 = R(i+1,j+1);        
      W(id1,id2) = 1;
      W(id1,id3) = 1;
      W(id1,id4) = 1;
      W(id1,id5) = 1;
      W(id1,id6) = 1;
      W(id1,id7) = 1;
      W(id1,id8) = 1;
      W(id1,id9) = 1;
    end;
end;

%for four coners: each has three neighbors
    id1=R(1,1);
    id2=R(1,2);
    id3=R(2,1);
    id4=R(2,2);
    W(id1,id2) = 1;
    W(id1,id3) = 1;
    W(id1,id4) = 1;

    id1=R(1,m);
    id2=R(1,m-1);
    id3=R(2,m);
    id4=R(2,m-1);
    W(id1,id2) = 1;
    W(id1,id3) = 1;
    W(id1,id4) = 1;

    id1=R(k,1);
    id2=R(k-1,1);
    id3=R(k,2);
    id4=R(k-1,2);
    W(id1,id2) = 1;
    W(id1,id3) = 1;
    W(id1,id4) = 1;

    id1=R(k,m);
    id2=R(k,m-1);
    id3=R(k-1,m);
    id4=R(k-1,m-1);
    W(id1,id2) = 1;
    W(id1,id3) = 1;
    W(id1,id4) = 1;

%for first and last rows: each has five neighbors
for j=2:(m-1);      
    id1 = R(1,j);   %first row neighbors: left, right and below
    id2 = R(1,j-1);
    id3 = R(1,j+1);
    id4 = R(2,j);
    id5 = R(2,j-1);
    id6 = R(2,j+1);

    W(id1,id2) = 1;
    W(id1,id3) = 1;
    W(id1,id4) = 1;
    W(id1,id5) = 1;
    W(id1,id6) = 1;
    
    id1 = R(k,j);   %las row neighbors: left, right and above
    id2 = R(k,j-1);
    id3 = R(k,j+1);
    id4 = R(k-1,j);
    id5 = R(k-1,j-1);
    id6 = R(k-1,j+1);

    W(id1,id2) = 1;
    W(id1,id3) = 1;
    W(id1,id4) = 1;
    W(id1,id5) = 1;
    W(id1,id6) = 1;
end;

%for first and last colmns: each has three neighbors
for i=2:(k-1);
    id1 = R(i,1);   %first column neighbors: above, below and right
    id2 = R(i-1,1);
    id3 = R(i+1,1);
    id4 = R(i,2);
    id5 = R(i-1,2);
    id6 = R(i+1,2);

    W(id1,id2) = 1;
    W(id1,id3) = 1;
    W(id1,id4) = 1;
    W(id1,id5) = 1;
    W(id1,id6) = 1;
    
    id1 = R(i,m);   %last column neighbors: above, below and right
    id2 = R(i-1,m);
    id3 = R(i+1,m);
    id4 = R(i,m-1);
    id5 = R(i-1,m-1);
    id6 = R(i+1,m-1);

    W(id1,id2) = 1;
    W(id1,id3) = 1;
    W(id1,id4) = 1;
    W(id1,id5) = 1;
    W(id1,id6) = 1;
end;
for i=1:n;
    W(i,:) = W(i,:)/sum(W(i,:));
end;

f = W;    % Spatial Weight Matrix