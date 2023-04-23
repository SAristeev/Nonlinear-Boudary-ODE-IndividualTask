clear;
clc; clf;
% Nonlilear boundary problem
% second order
% y'' = f(x,y,y') x at [xL,xR]
% y(xL) = cL
% y(xR) = cR
% Look f.m to see the function
%
% Solving on coarest mesh M points
% mesh shifted to the left by h/2
% Finding Newton's corrections on finer mesh 
% and interpolate it to coarest


%% input parameters
xL = 0;
xR = 10;
cL = 0;
cR = tanh(10);
scope = 3;
eps = 1e-4;

%% Coarest mesh
M = 128;
h = (xR - xL)/(M - 1.5);

alpha = 0.5;
x = linspace(xL-h/2,xR,M);
yinit = (1-alpha)*tanh(x) + alpha * tanh(10) * x/10;
y = yinit;

%% display first plot
subplot(2,3,[2,3,5,6])
plot(x,yinit,'LineWidth',1, 'Color','R')
hold on
plot(x,tanh(x),'LineWidth',1,'Color','G')
hold off

xlabel("x");
ylabel("y");
ylim([0.0 1.1])
xlim([0.0 10.0])
title("Solution Initialization");
legend("y_{init}", "y_{analytic}", 'Location','southeast')

%% init deltas
accuracy = 1;
iter = 1;
while accuracy > eps
    %% init locals to internal loop - find v - Newtonian correction
    Mv = M; hv = h; xv = x; yv = y;
    mesh_num = 1;

    %% first solve
    v = findCorrection(xv,yv,Mv,hv);
    delta = norm(v,"inf");
    while delta > eps/2
        %% transform to fine mesh
        hv = hv/scope;
        Mv = (Mv - 1) * scope;

        xfine = linspace(xL-hv/2,xR,Mv);
        yfine = Lagrange(xv,yv,xfine,3);
        vprev = Lagrange(xv,v,xfine,3);
        
        %% find new correction
        v = findCorrection(xfine,yfine,Mv,hv);
        xv = xfine;
        yv = yfine;
        
        deltaprev = delta;
        delta = norm(vprev-v,"inf");
        reldelta = deltaprev / delta;

        %% display plot
        hold on
        subplot(2,3,1)
        plot(xv,v);
        hold off
        xlim([0.0 10.0])
        xlabel("x");
        ylabel("v");
        title(sprintf("Newton iteration %d mesh = %d M = %5d\nerr = %5.2e order = %f", iter, mesh_num, Mv, delta, reldelta));
        drawnow

        %% print info
        fprintf("iter = %d  mesh_num = %d   M = %5d  ", iter, mesh_num, Mv);
        fprintf("err = %5.2e  required accuracy = %5.2e   order = %f\n", delta, eps/2,reldelta);
        mesh_num = mesh_num + 1;
    end %while
    %% interpolate correction and change solution
    v_coarse = Lagrange(xfine,v,x,3);
    y = y + v_coarse;
    accuracy = norm(v_coarse,"inf");

    %% display plots
    % main plot
    hold on
    subplot(2,3,[2,3,5,6])
    plot(x,y,'DisplayName',['y_' num2str(iter)]);
    hold off

    title(sprintf("Newton iteration %d, ||v|| = %f",iter,accuracy));
    
    % Convergence plot
    hold on
    subplot(2,3,4)
    plot(iter,log10(accuracy),'Marker','o', 'Color','red');
    hold off

    xlabel("iteration");
    ylabel("log_{10}(||v||)");
    title(sprintf("Convergence of Newtonian iterations, iteration %d",iter));
    drawnow

    %% print info
    fprintf("For Newtonian iteration %2d norm of residual ||v|| = %6.4e  \n", iter, accuracy);
    iter = iter + 1;
end % for
fprintf("Newtonian iterations have SUCESSFULLY finished for alpha = %f \n", alpha);