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
yL = 0;
yR = tanh(10);


%% Base mesh and init 

M = 128;
h = (xR - xL)/(M - 1.5);

% eps = 1e-4;
% max 
% alpha = 1.78;
% min
% alpha = -8.53;

alpha = 1;
xbase = linspace(xL-h/2,xR,M);
yinit = (1-alpha)*tanh(xbase) + alpha * tanh(10) * xbase/10;
ybase = yinit;

eps = 1e-4;
scope = 2;
MaxNewtonsIter = 100;
MaxMeshPoints = M * scope^8;

%% create subplot's handles

mainPlot = subplot(2,3,[2,3,5,6]);
correctionPlot = subplot(2,3,1);
convergencePlot = subplot(2,3,4);

hold(mainPlot,'on');
hold(correctionPlot,'on');
hold(convergencePlot,'on');

%% display first plot

subplot(mainPlot);
plot(xbase,yinit,'LineWidth',2, 'Color','C')

plot(xbase,tanh(xbase),'LineWidth',2,'Color','B')
plot([xL,xR], [yL,yR], 'ob');

xlabel('x');
ylabel('y');
xlim([xL xR])
title('Solution Initialization');
legend('y_{init}', 'y_{analytic}', 'Boundary conditions', 'Location','southeast')

%% print headlines

fprintf('alpha = %f eps = %.1e\nxL = %f xR = %f\nyL = %f yR = %f\n', alpha, eps, xL, xR, yL, yR);
fprintf('Max Newtonian iterations = %d\n',MaxNewtonsIter);
fprintf('Max mesh points = %d\n',MaxMeshPoints);
fprintf('Start solve\n');
fprintf('|------|------|--------|----------|----------|-----------|\n');
fprintf('| iter | mesh | points |  delta   |  eps/2   |   ratio   |\n');
fprintf('|------|------|--------|----------|----------|-----------|\n');

%% prepare to loop

normv = 1;
iter = 1;
isSucessFindIter = true;
isSucessFindFunc = true;

tStart=tic;
while normv > eps/2
    %% init locals to internal loop - find v - Newtonian correction

    Mfine = M; hcoarse = h; xcoarse = xbase; ycoarse = ybase;
    mesh_num = 1;
    if iter > MaxNewtonsIter
        isSucessFindFunc = false;
        break
    end % if

    %% first solve

    vbase = findCorrection(xbase,ybase,M,h, yL);
    vcoarse = vbase;
    delta = norm(vbase,'inf');
    
    
    cla(correctionPlot);
    while delta > eps/2
        %% transform to fine mesh
        Mfine = Mfine * scope;
        hfine = (xR - xL)/(Mfine - 1.5);

        if Mfine > MaxMeshPoints
            isSucessFindIter = false;
            break
        end % if

        xfine = linspace(xL-hfine/2,xR,Mfine);
        yfine = Lagrange(xcoarse,ycoarse,xfine,4);
        vcoarse2fine = Lagrange(xcoarse,vcoarse,xfine,4);
        
        %% find new correction
        
        vfine = findCorrection(xfine,yfine,Mfine,hfine, yL);
        deltaprev = delta;
        
        delta = norm(vcoarse2fine-vfine,'inf');
        reldelta = deltaprev / delta;

        %% display plot and print info
        subplot(correctionPlot)
        
        plot(xfine,vfine,'DisplayName',[[['mesh ' num2str(mesh_num)] ', M = '] num2str(Mfine)]);
        legend()

        xlim([0.0 10.0])
        xlabel('x');
        ylabel('v');
        title(texlabel(sprintf('Newton iteration %d, mesh = %d, M = %5d\n delta_{%d} = %5.2e\n delta_{%d}/delta_{%d} = %f', iter, mesh_num, Mfine, mesh_num, delta,mesh_num,mesh_num-1, reldelta)));
        drawnow
        
        if mesh_num == 1
            fprintf('| %4d | %4d | %6d | %5.2e | %5.2e |   -----   |\n', iter, mesh_num, Mfine, delta, eps/2);
        else
            fprintf('| %4d | %4d | %6d | %5.2e | %5.2e | %9.4f |\n', iter, mesh_num, Mfine, delta, eps/2, reldelta);
        end
        xcoarse = xfine;
        ycoarse = yfine;
        vcoarse = vfine;
        Mcoarse = Mfine;
        hcoarse = hfine;
        mesh_num = mesh_num + 1;
    end %while

    %% interpolate correction and change solution

    vfine2base = Lagrange(xfine,vfine,xbase,4);
    ybase = ybase + vfine2base;
    normv = norm(vfine,'inf');

    %% display plots and print info
    
    % main plot
    subplot(mainPlot)
    plot(xbase,ybase,'DisplayName',[['y_{' num2str(iter)] '}']);
    title(sprintf('Current approximate\nNewton iteration %d, ||v|| = %4.2e',iter,normv));
    
    % Convergence plot
    subplot(convergencePlot)
    plot(iter,log10(normv),'Marker','o', 'Color','B', 'DisplayName','||v||');
    
    
    xlabel('Iteration number');
    ylabel('log_{10}(||v||)');
    title(sprintf('Convergence of Newtonian iterations, iteration %d',iter));

    drawnow

    if(isSucessFindIter)
        fprintf('| %4d |-------------------------------------------------|\n', iter);
        fprintf('| %4d | Find Newtonian iteration - success              |\n', iter);
        fprintf('| %4d | ||v|| = %4.2e   eps/2 = %5.2e             |\n', iter, normv, eps/2);
        fprintf('|------|------|--------|----------|----------|-----------|\n');
    else
        fprintf('| %4d |--------------------------------------------------\n', iter);
        fprintf('| %4d | Find Newtonian iteration - unsuccess            |\n', iter);
        fprintf('| %4d | Point limit exceeded                            |\n', iter);
        fprintf('| %4d | Maximum = %5d, current = %5d                |\n', iter, MaxMeshPoints, Mfine);
        fprintf('| %4d | ||v|| = %4.2e   eps/2 = %5.2e             |\n', iter, normv, eps/2);
        fprintf('|------|------|--------|----------|----------|-----------|\n');
        break;
    end %if
    
    iter = iter + 1;
end % for
tElapsed=toc(tStart);

if(isSucessFindFunc && isSucessFindIter)
    subplot(mainPlot)
    plot(xbase,ybase,'LineWidth',2,'Color','G','DisplayName','y_{result}');
    subplot(convergencePlot)
    plot(iter,log10(normv),'Marker','o', 'Color','G', 'DisplayName','||v||');
    
    fprintf('Solve finished SUCESSFULLY, time = %f\n', tElapsed);
else
    subplot(mainPlot)
    plot(xbase,ybase,'LineWidth',2,'Color','R','DisplayName','y_{result}');
    subplot(convergencePlot)
    plot(iter,log10(normv),'Marker','o', 'Color','R', 'DisplayName','||v||');
    
    
    fprintf('Solve finished UNSUCESSFULLY, time = %f\n', tElapsed);
    if(isSucessFindIter)
        fprintf('Exceeded limit of Newtonian iterations\n');
    end
    
end % if

hold(mainPlot,'off');
hold(correctionPlot,'off');
hold(convergencePlot,'off');
