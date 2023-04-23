function [D] = diffMatrix(x,L,order)
%build L-diff matrix by order
%scheme and pattern for matrix build by order
maxprime = 0;
for i = 1:length(L)
    if L(i) ~=0
        maxprime = i;
    end
end

copyL = zeros(1,maxprime+order-1);
for i=1:maxprime
    copyL(i) = L(i);
end

pattern = 0:length(copyL)-1;
D = sparse(length(x),length(x));


if mod(length(pattern),2) == 0
    j_=1:length(copyL);
    for i = 1:(length(pattern)/2-1)
        c = computeCoeff(x, pattern, copyL, i);
        i_ =ones(1,length(copyL))*i;
        D = D + sparse(i_,j_,c.',length(x),length(x));
        pattern = pattern-1;
    end
    %pattern = -(length(pattern)/2-1):length(pattern)/2;
    %currow = length(pattern)/2;
    for i = length(pattern)/2:length(x)-length(pattern)/2
        c = computeCoeff(x, pattern, copyL, i);
        i_=ones(1,length(copyL))*i;
        
        D = D + sparse(i_,j_,c.',length(x),length(x));
        j_=j_+1;
    end
    j_=j_-1;
    for i = length(x)-length(pattern)/2 + 1:length(x)
       pattern = pattern - 1;
       i_=ones(1,length(copyL))*i;
       c = computeCoeff(x, pattern, copyL, i);
       D = D + sparse(i_,j_,c.',length(x),length(x));
    end
else
    for i = 1:((length(pattern)-1)/2)
        c = computeCoeff(x, pattern, copyL, i);
        i_ =ones(1,length(copyL))*i;
        j_=1:length(copyL);
        D = D + sparse(i_,j_,c.',length(x),length(x));
        pattern = pattern-1;
    end
    for i =((length(pattern)-1)/2 + 1):(length(x)-(length(pattern)-1)/2)
        c = computeCoeff(x, pattern, copyL, i);
        i_ =ones(1,length(copyL))*i;
        j_ = (0:(length(copyL)-1)) +(i-(length(pattern)-1)/2);
        D = D + sparse(i_,j_,c.',length(x),length(x));
    end
    for i =(length(x)-(length(pattern)-1)/2 + 1):length(x)
        pattern = pattern - 1;
        c = computeCoeff(x, pattern, copyL, i);
        i_= ones(1,length(copyL))*i;
        D = D + sparse(i_,j_,c.',length(x),length(x));
    end
end
