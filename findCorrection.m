function [v] = findCorrection(xv,yv,Mv,hv)
    D2 = diffMatrix(xv,[0,0,1],3);
    D1 = diffMatrix(xv,[0,1,0],3);
    D0 = diffMatrix(xv,[1,0,0],3);
    fi = zeros(Mv,1);
    for m=2:Mv-1
        D1(m,:) = -D1(m,:) * fyprime(xv(m),yv(m),(yv(m+1)-yv(m-1))/(2*hv));
        D0(m,:) = -D0(m,:) * fy(xv(m),yv(m),(yv(m+1)-yv(m-1))/(2*hv));
        fi(m) = f(xv(m),yv(m),(yv(m+1)-yv(m-1))/(2*hv)) - (yv(m+1)-2*yv(m)+yv(m-1))/hv^2;
    end % for m
    D = D2 + D1 + D0;

    % boundary conditions
    E = eye(Mv);
    D(1,:) = E(1,:);    fi(1) = 0;
    D(Mv,:) = E(Mv,:);  fi(Mv) = 0;
    % solve
    v = (D\fi)';
end

