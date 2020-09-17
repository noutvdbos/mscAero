classdef LiftingLine
    
%This model is based on the lectures from Rotor/wake Aerodynamics at TU 
%Delft. The code is produced by Nout van den Bos, and Floyd van Steen. The
%Lifting Line model is inspired by Low Speed Aerodynamics, Katz.
    
properties
    Uinf
    N 
    c
    radius
    rTip
    rRoot
    omega
    distType
    propType
    excel
    Nblades
    Nrotors
    Drotors
    rotorPhase
    chordfun
    twistfun
end

methods
    
    %constructor
    function obj = LiftingLine(Uinf,N,radius,rTip,rRoot,TSR,distType,...
                               propType,polar,Nblades,Nrotors,Drotors,...
                               rotorPhase, chord,twist)
        obj.Uinf = Uinf;
        obj.N = N;        
        obj.radius = radius;
        obj.rTip = rTip;
        obj.rRoot = rRoot;
        obj.omega   = Uinf*TSR/radius;
        obj.distType = distType;
        obj.propType = propType;
        obj.excel = xlsread(polar);
        obj.Nblades = Nblades;
        obj.Nrotors = Nrotors;
        obj.Drotors = Drotors;
        obj.rotorPhase = deg2rad(rotorPhase); %convert to radians
        obj.chordfun = chord;
        obj.twistfun = twist;
    end
    
    % setup of spanwise distribution for uniform/cosine
    function [r] = spanwiseDistr(obj) 
        if obj.distType == "uniform"
            r = obj.rRoot:(obj.rTip-obj.rRoot)/(obj.N):obj.rTip;
        elseif obj.distType == "cosine"  
            angle = 0.0:pi/(obj.N):pi;
            r = flip(cos(angle));
            r = r*(obj.rTip-obj.rRoot)/2+(obj.rTip-obj.rRoot)/2+0.2;
        end
    end
    
    %this function calculates all the control points, as well as the
    %geometry of the wake. This is necessary to calculate the induced
    %velocities from the vortices in the wake.
    function [coord,control,r_cp] = modelSetUp(obj,LWake,aWake,deltaDeg,expansion)
        
        % Obtain spanwise distribution
        r = spanwiseDistr(obj);
        r_cp = zeros(length(r)-1,1);
        
        for j=1:length(r)-1
            r_cp(j) = (r(j) + r(j+1))/2;
        end       
        
        %Angle of each blade
        bladeAngle = 0:2*pi/(obj.Nblades):2*pi; 
        
        %Control points at quarterchord assuming no quarter chord sweep
        control = zeros(obj.N*obj.Nblades*obj.Nrotors,3);
        
        %setting up control points for all blades
        for i = 1:obj.Nblades*obj.Nrotors
            istart = (i-1)*obj.N +1;
            iend   = istart + obj.N -1;            
            if i <= obj.Nblades
                control(istart:iend,2) = obj.radius*r_cp*cos(bladeAngle(i));
                control(istart:iend,3) = obj.radius*r_cp*sin(bladeAngle(i));
            else 
                j = i - obj.Nblades;
                control(istart:iend,2) = obj.Drotors + obj.radius*...
                    r_cp*cos(obj.rotorPhase+bladeAngle(j));
                
                control(istart:iend,3) = obj.radius*r_cp*...
                    sin(obj.rotorPhase+bladeAngle(j));       
            end
        end
        
        %Initial positions of the blades
        coordInit = zeros(obj.N+1,3,obj.Nblades*obj.Nrotors);
        
        for i=1:obj.Nblades*obj.Nrotors
            if i <= obj.Nblades
                coordInit(:,2:3,i) = [obj.radius*r'*cos(bladeAngle(i)),...
                    obj.radius*r'*sin(bladeAngle(i))];
            else
                j = i-obj.Nblades;
                coordInit(:,2:3,i) = ...
                [obj.Drotors + obj.radius*r'*cos(obj.rotorPhase+...
                bladeAngle(j)),obj.radius*r'*sin(obj.rotorPhase+bladeAngle(j))];
            end
        end
        
        chord = obj.chordfun(r);
        twist = obj.twistfun(r);
        
        %Obtaining the angles in radians for the wake discretization
        Uwake =obj.Uinf*(1-aWake);                  % Wake convection speed
        t_OneRotation = 2*pi/obj.omega;             % time required for a single rotation
        L_OneRotation = t_OneRotation*Uwake;        % Length single rotation
        Nr_rotations = LWake./L_OneRotation;        %Number of rotations
        x = max(Nr_rotations)*L_OneRotation;        %Length of wake
        
        % wake discretisation in radians
        discr = 0:deg2rad(deltaDeg):2*pi*max(Nr_rotations)+deg2rad(deltaDeg);
        
        %First point should be at 1.25c from the leading edge (1c from control point)
        %Twist + pitch is taken to find the spanwise z coordinates of these points
        z0 = chord.*cos(deg2rad(twist));
        
        %The phase is when the point is reached by the blade after this the wake discretization starts
        phase = atan(z0./(r*obj.radius));
        
        %Setup 4D coordinates matrix
        %Columns are the positions of the blade in the wake for a spanwise point
        %Rows are spanwise positions
        %3rd dimension is x,y,z-coordinate
        %4th dimension is the bladenumber          
        coord = zeros(length(discr)+1,obj.N+1,3,obj.Nblades*obj.Nrotors);    
        
        xbar = (discr+phase(:))./(2*pi*Nr_rotations).*x./(r(:).*obj.radius);
        
        % Include wake expansion according to method explained in Wind 
        % Turbine Aerodynamics and Vorticity-Based Methods, section 15.1.
        
        Rwake     = zeros(length(r),length(discr));
        Vdecrease = zeros(1,length(discr));
        
        if (expansion == true)
            Rwake(1,:) = ones(1,length(discr))*r(1)*obj.radius;
            Rwake(2:end,:) = sqrt((1-aWake)*(1-aWake*(1+xbar(2:end,:)./...
                             sqrt(1+xbar(2:end,:).^2))).^-1).*r(2:end)'*obj.radius;
            Vdecrease = (obj.radius./Rwake(end,:)).^2;        
        else
           Rwake(:,:) = ones(1,length(discr)).*r(:).*obj.radius;
           Vdecrease = Vdecrease*0+1;
        end
        
        %Find coordinates of all discretised points in the wake
        for i=1:length(r)            
            for j=1:obj.Nblades*obj.Nrotors
                if j<=obj.Nblades
                    coord(2:end,i,1,j) = (discr+phase(i))/(2*pi*max...
                                         (Nr_rotations))*x.*Vdecrease;
                    coord(2:end,i,2,j) = Rwake(i,:).*cos(discr+phase(i)+...
                                         bladeAngle(j));
                    coord(2:end,i,3,j) = Rwake(i,:).*sin(discr+phase(i)+...
                                         bladeAngle(j));
                else
                    k = j-obj.Nblades;
                    coord(2:end,i,1,j) = (discr+phase(i))/(2*pi*max...
                                         (Nr_rotations))*x.*Vdecrease;
                    coord(2:end,i,2,j) = obj.Drotors + Rwake(i,:).*...
                                         cos(obj.rotorPhase+discr+...
                                         phase(i)+bladeAngle(k));
                    coord(2:end,i,3,j) = Rwake(i,:).*sin(obj.rotorPhase+...
                                         discr+phase(i)+bladeAngle(k));
                end         
            end
        end
        
        %Set first coordinates to the initial positions of the blade
        coord(1,:,:,:) = coordInit(:,:,:);
    
    end
    
    
    %calculates induced velocity due to single filament
    function [K,R12] = VelocitySingleFilament(obj,X1,X2,Xcp) 
       core = 0.0001;
       R1= norm(Xcp-X1);
       R2 = norm(Xcp-X2);
       R1vector = Xcp - X1;
       R2vector = Xcp - X2;
       
       R12 = cross(R1vector,R2vector);
       R12Magn = sum(R12.^2);
       R01 = dot((X2-X1),(Xcp-X1));
       R02 = dot((X2-X1),(Xcp-X2)); 
       
       if (R12Magn<core^2) 
           R12Magn = core^2;
       end
       if (R1<core) 
           R1 = core;
       end
       if (R2<core)
           R2 = core;
       end
       K = 1/(4*pi*R12Magn)*(R01/R1 - R02/R2);
    end
     
    
    function [u,v,w] = singleHshoe(obj,positions,cp)
       pos1 = positions(1:end-1,:);
       pos2 = positions(2:end,:);
       u = zeros(obj.N,size(pos1,1));
       v = zeros(obj.N,size(pos1,1));
       w = zeros(obj.N,size(pos1,1));
              
       for i=1:obj.N*obj.Nblades*obj.Nrotors %Loop through control points
           for j=1:size(pos1,1)              %Loop through each filament of the horsehoe
               
               [K,R12] = VelocitySingleFilament(obj, pos1(j,:),pos2(j,:),... 
                   cp(i,:));  

               u(i,j) = K*R12(1);       
               v(i,j) = K*R12(2);
               w(i,j) = K*R12(3);
           end           
       end
       
    end
    
    
    function [u,v,w] = getIndVel(obj,coordinates,cp)
        
        M = obj.N*obj.Nblades*obj.Nrotors;
        u = zeros(M,M);
        v = zeros(M,M);
        w = zeros(M,M);
        
        for j=1:obj.Nblades*obj.Nrotors          %Loop through the blades
            coord = coordinates(:,:,:,j);
            positions = zeros(2*size(coord,1),3);
            for i=1:obj.N
                
                %Constructing each horseshoe and order them in the right way
                for k=1:3
                    positions(:,k) = [flip(coord(:,i+1,k));coord(:,i,k)];
                end
                
                [uMatr,vMatr,wMatr] = singleHshoe(obj,positions,cp);
                
                %Sum all filaments and add this contribution to the induced V matrices
                u(:,i+(j-1)*obj.N) = sum(uMatr,2) + u(:,i+(j-1)*obj.N);
                v(:,i+(j-1)*obj.N) = sum(vMatr,2) + v(:,i+(j-1)*obj.N);
                w(:,i+(j-1)*obj.N) = sum(wMatr,2) + w(:,i+(j-1)*obj.N);
            end
        end
    end
    
   
    function [Fnorm,Ftan,gamma,alpha,inflowangle,cl] = ...
                getAirfoilForces(obj,Unorm,Utan,rCent)

        Vmag  = sqrt(Unorm^2+Utan^2);
        chord = obj.chordfun(rCent);
        twist = obj.twistfun(rCent);

        inflowangle = atand(Unorm/Utan);

        alpha = inflowangle - twist;

        if (obj.propType == "prop")
            if alpha > 10
                alpha = 10;
            elseif alpha <-25
                alpha = -25;
            end
            alphas = flipud(-obj.excel(:,1));    
            cls    = flipud(-obj.excel(:,2));    
            cds    = flipud(obj.excel(:,3));  
        
        elseif (obj.propType == "turb")
            alphas = obj.excel(:,1);    
            cls    = obj.excel(:,2);    
            cds    = obj.excel(:,3);    
        end 

        cl = interp1(alphas,cls,alpha);
        cd = interp1(alphas,cds,alpha);

        L = 0.5*cl*Vmag^2*chord;
        D = 0.5*cd*Vmag^2*chord;

        Fnorm = L*cosd(inflowangle)+D*sind(inflowangle);
        Ftan  = L*sind(inflowangle)-D*cosd(inflowangle);

        gamma = 0.5*cl*Vmag*chord;
    end
    
    
   function results = solveSystem(obj,U,V,W,cp,weight)
       
       %initialise the controlpoints on the blade
       r    = spanwiseDistr(obj);
       r_cp = zeros(length(r)-1,1);
       
       for j=1:length(r)-1
           r_cp(j) = (r(j) + r(j+1))/2;
       end
       
       if obj.Nrotors > 1
           cp = repmat(cp(1:obj.N*obj.Nblades,:),obj.Nrotors,1);
       end
       
       %establish the stopping criteria
       epsilon = 1E-5;
       maxIt   = 1E3;
       
       %initialise all the result arrays.
       resultsSize = obj.N*obj.Nblades*obj.Nrotors;
       gammaNew    = zeros(resultsSize,1);
       gamma       = zeros(resultsSize,1);
       Fnorm       = zeros(resultsSize,1);
       Ftan        = zeros(resultsSize,1);
       alpha       = zeros(resultsSize,1);
       cl          = zeros(resultsSize,1);
       inflowangle = zeros(resultsSize,1);
       a           = zeros(resultsSize,1);
       aprime      = zeros(resultsSize,1); 

       index = 0;
       while (1) && index<maxIt
           index = index + 1;

           for i = 1:size(gammaNew,1)
               gamma(i) = gammaNew(i);
           end

           for i = 1:resultsSize
               u = 0; v = 0; w = 0;
                
              for j = 1:resultsSize
                  u = u + U(i,j)*gamma(j);
                  v = v + V(i,j)*gamma(j);
                  w = w + W(i,j)*gamma(j);
              end
              
              %get index for r_cp 
              if i > obj.N
                if mod(i,obj.N) < 1E-4
                   idx = obj.N;
                else 
                   idx = mod(i,obj.N);
                end
              else
               idx = i;
              end
             
              %rotational velocity
              Urot  = cross([obj.omega,0,0],cp(i,:)); 

              %perceived velocity
              Uper = [obj.Uinf + u + Urot(1), v + Urot(2), w + Urot(3)]; 
              
              %Normal and tangential velocities
              Unorm = Uper(1);
              Utan  = sqrt(Uper(2)^2+Uper(3)^2); 
              
              [Fnorm(i),Ftan(i),gammaNew(i),alpha(i),inflowangle(i),cl(i)] = ...
                obj.getAirfoilForces(Unorm,Utan,r_cp(idx));
            
              a(i)      = -Unorm/obj.Uinf + 1;
              aprime(i) = Utan/(r_cp(idx)*obj.radius*+obj.omega)-1;
           end
           
           %check for convergence           
           error = abs(sum(gammaNew)-sum(gamma))/(max(max(gammaNew),epsilon));
           
           if (error<epsilon)
               disp("Number of iterations: " + num2str(index));
               break               
           end
           
           if (index == maxIt)
               disp("Solution did not converge.");
               break
           end
           
           %compose new gamma
           for i = 1:size(gammaNew,1)
               gammaNew(i) = gamma(i)*(1-weight) + gammaNew(i)*weight;               
           end
           
       end
       
       Cn = Fnorm/(0.5*obj.Uinf^2*obj.radius);
       Ct = Ftan/(0.5*obj.Uinf^2*obj.radius);     
       
       results = struct('a',a,'aprime',aprime,'Cn',Cn,'Ct',Ct,...
        'inflow',inflowangle,'AoA',alpha,'gamma',gamma,'Cl',cl);
   end
   
%% Coefficient calculating functions
   
   function CT = getCT(obj,r_cp,results)      
      %assume that it is uniformly distributed
      dr = r_cp(2)-r_cp(1);      
      CT =  dr*results.Cn.*obj.Nblades./pi; 
   end
   
   
   function CP = getCP(obj,r_cp,results)      
      %assume that it is uniformly distributed
      dr = r_cp(2)-r_cp(1);      
      rs = zeros(size(r_cp,1)*obj.Nblades*obj.Nrotors,1);
      for i = 1:obj.Nblades*obj.Nrotors
          sidx = (i-1)*obj.N +1;
          eidx = i*obj.N;
          rs(sidx:eidx) = r_cp;
      end      
      CP = dr.*results.Ct.*rs.*obj.Nblades*obj.radius*...
           obj.omega/(obj.Uinf*pi);
   end
   
   
   function Cgamma = getCGamma(obj,results)
       Cgamma = results.gamma*(obj.Nblades*obj.omega/(pi*obj.Uinf^2));
   end
   
%% Plotting functions

   function plotWake(obj,coordinates,colors)
        
        figure('DefaultAxesFontSize',12)
        
        for i=1:obj.Nblades*obj.Nrotors
            
            surf(squeeze(coordinates(:,:,1,i)),...
                 squeeze(coordinates(:,:,2,i)),...
                 squeeze(coordinates(:,:,3,i)),...
                 'FaceColor',string(colors{end+1-i,2}));
             
            hold on
        end

        axis 'equal'
        xlabel('x [m]')
        ylabel('y [m]')
        zlabel('z [m]')
   end
   
   
   function plotRotors(obj,coordinates,colors)
        figure('DefaultAxesFontSize',12)
        
        legstr = strings(1,obj.Nblades*obj.Nrotors);
        
        for i=1:obj.Nblades*obj.Nrotors
            
            legstr(i) = "Blade " + int2str(i);
            surf(squeeze(coordinates(1:2,:,1,i)),...
                 squeeze(coordinates(1:2,:,2,i)),...
                 squeeze(coordinates(1:2,:,3,i)),...
                 'FaceColor',string(colors{end+1-i,2}));
             
            hold on
            
        end

        ycoords = coordinates(:,:,2,:);
        zcoords = coordinates(:,:,3,:);
        ylim([min(ycoords(:)),max(ycoords(:))]);
        zlim([min(zcoords(:)),max(zcoords(:))]);
        
        legend(legstr);
        view(90,0)
        axis 'equal'
        xlabel('x [m]')
        ylabel('y [m]')
        zlabel('z [m]')
   end
   
    
   function plotResults(obj,bladeidx,r_cp,result,ylab,leg,colors)
        
        figure('DefaultAxesFontSize',12)
        
        for j = 1:size(result,2)
            for i = 1:size(bladeidx,2)
                
                idx      = bladeidx(i);
                startidx = 1 + (idx-1)*obj.N;
                endidx   = idx*obj.N;
            
                if size(leg)
                    bladestr = leg(j) + ": Blade " + int2str(idx);
                else
                    bladestr = "Blade " + int2str(idx);
                end
                
                plot(r_cp,result(startidx:endidx,j),'-x','DisplayName',...
                     bladestr,'Color',string(colors{end+1-i,2}));
                hold on; 
            end
        end

        xlabel('r/R [-]')
        ylabel(ylab)
        legend('show',"interpreter", "latex")
   end
   
   
end
end


