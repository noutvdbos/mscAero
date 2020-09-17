classdef PanelMethod
    
% This code is based on the methods from Low Speed Aerodynamics, by Katz.
% Mostly from chapter 5 and 13. All notations such as xg, ac, nc etc. are
% explained in this book. Before going through this code, i recommend
% reading these chapters.

    
properties
    chord
    nodes
    xg
    xc
    ac
    nc
    tc
    Uinf
    rho
end

methods (Static)
      
function obj = PanelMethod(nodes,Uinf,rho)
    
    obj.nodes = nodes;
    obj.Uinf  = Uinf;
    obj.rho   = rho;
    obj       = obj.getGeom();
    obj.chord = norm(nodes(:,1)-nodes(:,end));
end

function U = vor2D(X,Xj,gammaj) 

   core = 0.0001;
   rj  = vecnorm(X-Xj);
   
   if (rj <core)
       rj = core;
   end
   mat = [0,1;-1,0];
   
   U = gammaj/(2*pi*rj^2).*mat*(X-Xj);

end
    
end  

methods

% This function calculates the control points, the vortex points, and
% properties of these points. These points represent the geometry.

function obj = getGeom(obj)
        
    obj.xg = 0.75*obj.nodes(:,1:end-1) + 0.25*obj.nodes(:,2:end);
    obj.xc = 0.25*obj.nodes(:,1:end-1) + 0.75*obj.nodes(:,2:end);
    
    obj.ac = -atan((obj.nodes(2,2:end)-obj.nodes(2,1:end-1))./...
                   (obj.nodes(1,2:end)-obj.nodes(1,1:end-1)));
          
    obj.nc = [sin(obj.ac);cos(obj.ac)];
    obj.tc = [cos(obj.ac);-sin(obj.ac)];
    
end   

%Steady state solver, for explanation of the algorithm, read section
%5.1-5.2.

function results = solveSteady(obj)

    mat    = obj.getCoefmat();
    rhs    = obj.getRHS(obj.Uinf);
    
    deltachord = vecnorm(obj.nodes(:,2:end)-obj.nodes(:,1:end-1))';
    
    results.gammas = linsolve(mat,rhs);
    results.Cls    = results.gammas./(0.5*norm(obj.Uinf)*obj.chord);
    results.Cps    = results.gammas./(0.5*norm(obj.Uinf).*deltachord);
    results.Cl     = sum(results.Cls);
    
end


%The unsteady solver is similar to the steady solver, with an additional
%matrix column. The small x's indicate a local reference frame, the large
%X/s indicate the global reference frame.

%thetafun and thetadotfun are the functions of how theta and thetadot
%change over time. They are user input but for the simple example case
%these functions are set to zero.

function results = solveUnsteady(obj,tend,dt,thetafun,thetadotfun)
    
    prevGamma = 0;
    t = dt:dt:tend;
    
    %initialise the coordianates of the wake, and the wake strength
    Xw = [];
    gw = [];
    
    %Initialise the result matrices
    Cls    = zeros(length(t),1);
    Cds    = zeros(length(t),1);
    gammas = zeros(length(t),length(obj.xg));
    
    %generate the coefficient matrix for unsteady case, it is modified from
    %the steady case
    xw   = obj.nodes(:,end)+0.2*dt*obj.Uinf;
    vecw = zeros(length(obj.xc),1);
    
    for i = 1:size(obj.xc,2)
        U = obj.vor2D(obj.xc(:,i),xw,1);
        vecw(i) = dot(U,obj.nc(:,i));
    end
    
    %similar to steady case, but an additional vector needs to be added to
    %the matrix.
    mat = obj.getCoefmat(); 
    mat = [mat vecw];
    mat = [mat;ones(1,length(mat))];
    
    %initialise new coordinate system
    theta0 = thetafun(0);
    transmat = [cos(theta0),sin(theta0);...
               -sin(theta0),cos(theta0)];
           
    Xc     = transmat*obj.xc;
    Xg     = transmat*obj.xg;
    Xnodes = transmat*obj.nodes;
    
    for idx = 1:length(t) 
        
        theta    = thetafun(t(idx));
        thetadot = thetadotfun(t(idx));
        
        %update the coordinate positions
        Xc     = Xc     - obj.Uinf*dt;
        Xg     = Xg     - obj.Uinf*dt;
        Xnodes = Xnodes - obj.Uinf*dt;

        [U,Uwake]   = obj.getVel(theta,thetadot,Xc,Xw,gw);
        
        %the RHS is also similar to that of the steady case, but another
        %entry needs to be added.
        rhs = obj.getRHS(U);
        rhs = [rhs;prevGamma];

        tempgammas = linsolve(mat,rhs);
        
        %save results
        gammas(idx,:)       = tempgammas(1:end-1);
        [Cls(idx),Cds(idx)] = obj.getClCd(U,Uwake,gammas,idx,dt);
        
        %update the wake
        for i = 1:size(Xw,2)
            indU  = obj.getIndVel(Xw(:,i),Xg,tempgammas(1:end-1));
            indUw = obj.getIndVel(Xw(:,i),Xw,gw);
            Xw(:,i) = Xw(:,i) + dt*(indU+indUw); 
        end
       
        Xwnew = Xnodes(:,end)+0.2*dt*obj.Uinf;
        Xw = [Xw Xwnew];            %add newest vortex position
        gw = [gw tempgammas(end)];  %add vortex strength
        
        results.XwHist{idx} = Xw;
        prevGamma = sum(tempgammas(1:end-1));
    end
    
    %save the results in a struct.
    results.gammas = gammas;
    results.Xg  = Xg;
    results.Xw  = Xw;
    results.gw  = gw;
    results.Cls = Cls;
    results.Cds = Cds;
    
end

%read section 13.8 for an explanation of this algorithm
function [Cl,Cd] = getClCd(obj,U,Uwake,gammas,idx,dt)
    
    dps    = zeros(length(obj.xg),1);
    Ds     = zeros(length(obj.xg),1);
    dchord = vecnorm(obj.nodes(:,2:end)-obj.nodes(:,1:end-1))';
    
    if (idx>1)
        gammaprev = gammas(idx-1,:); 
    else
        gammaprev = (gammas(idx,:));
    end
    
    if (size(U,2)==1)
        U = U.*ones(size(obj.nc));
    end
    
    for i = 1:length(dps)
        dgammadt = (sum(gammas(idx,1:i))-sum(gammaprev(1:i)))/dt;
        dps(i)   = obj.rho*(norm(U(:,i))*gammas(idx,i)/...
                   dchord(i) + dgammadt);
        Ds(i)    = obj.rho*(Uwake(2,i)*gammas(idx,i) + dgammadt*...
                   dchord(i)*sin(obj.ac(i)));
    end
    
    D = sum(Ds);
    L = sum(dps.*dchord.*cos(obj.ac)');
    
    Cl = L/(0.5*obj.rho*norm(obj.Uinf)^2*obj.chord);
    Cd = D/(0.5*obj.rho*norm(obj.Uinf)^2*obj.chord);
end


function indU = getIndVel(obj,x,xg,gammas)
    
    indU = [0;0];
    
    for i = 1:size(xg,2)
            indU = indU + obj.vor2D(x,xg(:,i),gammas(i));
    end
end


function [U,Uwake] = getVel(obj,theta,thetadot,Xc,Xw,gw)
    
    Uwake = zeros(2,length(Xc));
    for i = 1:length(Xc) 
        Uwake(:,i) = obj.getIndVel(Xc(:,i),Xw,gw);
    end
    
    mat = [cos(theta), -sin(theta); sin(theta), cos(theta)]; 
    U = Uwake + mat*obj.Uinf + [-thetadot*obj.xc(2,:);thetadot*obj.xc(1,:)];
    
end


function [mat] = getCoefmat(obj)
    
    mat = zeros(length(obj.xc),length(obj.xg));

    for i = 1:size(obj.xc,2)
        for j = 1:size(obj.xg,2)
            U = obj.vor2D(obj.xc(:,i),obj.xg(:,j),1);
            mat(i,j) = dot(U,obj.nc(:,i));
        end
    end
end


function [vec] = getRHS(obj,U)
    
    if (size(U,2)==1)
        U = U.*ones(size(obj.nc));
    end
    
    vec = zeros(length(obj.nc),1);
    
    for i = 1:size(obj.nc,2)
        vec(i) = -dot(U(:,i),obj.nc(:,i)); 
    end
end


%% Results calculations

function pres = getPressureField(obj,U,W,pref)
    
    % use bernouilli
    ptot = 0.5*obj.rho*norm(obj.Uinf)^2 + pref;
    pres = ptot - 0.5*obj.rho.*(U+W).^2;
    
end

function [U,W] =  getVelocityField(obj,x,y,xg,gammas)
   
    U = zeros(length(x),length(y));
    W = zeros(length(x),length(y));
    
    for i = 1:length(x)
        for j = 1:length(y)
            xvec    = [x(i);y(j)];
            indU = obj.getIndVel(xvec,xg,gammas);
            U(i,j) = indU(1) +obj.Uinf(1);
            W(i,j) = indU(2) +obj.Uinf(2);
        end
    end
    
end

%% plotting functions

function plotField(obj,X,Y,field,clab,offset,clims)
   
    figure();
    hold on;
    s = pcolor(X,Y,field');
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    colormap('jet')

    plot(obj.nodes(1,:)+offset(1),obj.nodes(2,:)+offset(2),...
         "k","LineWidth",2);
    
    c=colorbar();
    caxis(clims);
    c.Label.String = clab;
    
    xlabel("x [m]");
    ylabel("z [m]");
    
end

function saveVortex(obj,dt,Lpanel,results,plotName)
    
    hfig = figure('visible','off');

    loop = length(results.XwHist);
    for i=1:loop
        
        %plot the airfoil camber line
        plot(obj.nodes(1,:),obj.nodes(2,:),...
             'Color','red','LineWidth',2)
        hold on    
        
        %plot the waketrail at this point in time
        plot(results.XwHist{i}(1,:)+dt*i*obj.Uinf(1),...
            results.XwHist{i}(2,:),'Color','black','LineWidth',1.5)
        
        %add text about time and Cl
        text(0,Lpanel,['Time = ',num2str(round(i*dt,2)),' s'],'FontSize',14)
        text(5,Lpanel,['Cl = ',num2str(round(results.Cls(i),2)),' [-]'],...
            'FontSize',14)
        
        %make sure that the x-axis can fit the entire wake
        axis([-.5*Lpanel Lpanel + obj.Uinf(1)*loop*dt ...
              -2*Lpanel 2*Lpanel]);
          
        pbaspect([3 1 1])

        drawnow
        hold off
        F(i) = getframe(hfig);
    end

    %save the results as a video
    vid = VideoWriter(plotName);
    vid.FrameRate = 20;
    open(vid);
    writeVideo(vid,F);
    close(vid);

end
    
end
end