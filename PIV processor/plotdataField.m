%% Made by Tullio Natta

function [res] = plotdataField(folder,file_name,plot_points,mask,quivr)

    res = readtable(strcat(folder, file_name));
    res = res{:,1:end};

    %apply the mask for bad data
    indices = find(res(:,end)==0);
    res(indices,3:4) = NaN;

    %perform interpolation to match the plot points
    interp_u = scatteredInterpolant(res(:,1),res(:,2),res(:,3),'linear','none');
    interp_v = scatteredInterpolant(res(:,1),res(:,2),res(:,4),'linear','none');
    xs = linspace(min(res(:,1)),max(res(:,1)), round(plot_points*2.027));
    ys = linspace(min(res(:,2)),max(res(:,2)), round(plot_points));
    xs_full = linspace(min(res(:,1)),max(res(:,1)), round(plot_points*2.027));
    ys_full = linspace(min(res(:,2)),max(res(:,2)), round(plot_points));

    %create the meshgrid needed for plotting
    [X,Y] = meshgrid(xs,ys);
    [X_full,Y_full] = meshgrid(xs_full,ys_full);
    U = interp_u(X,Y);
    V = interp_v(X,Y);
    U_full = interp_u(X_full,Y_full);

    %insert the mask for the airfoil shape
    mask(:,1) = -mask(:,1);
    pgon = polyshape(mask);
    pgon = rotate(pgon,180);
    pgon = scale(pgon,0.08);
    pgon = translate(pgon,min(res(:,1)),max(res(:,2)));

    figure();
    
    %create the contours
    s =  pcolor(X_full,Y_full,U_full);
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    colormap('jet')
    c = colorbar;
    c.Label.String = 'U [m/s]';
    hold on

    %create the quivers
    quiver(X(1:quivr:end,1:quivr:end),Y(1:quivr:end,1:quivr:end),...
           U(1:quivr:end,1:quivr:end),...
           V(1:quivr:end,1:quivr:end),1, 'color',[0 0 0])
       
    set(gca, 'YDir','reverse')
    hold on
    plot(pgon,'FaceColor','black','FaceAlpha',1);
    xlabel('x [mm]'); ylabel('y [mm]')
    title('Velocity vectors and color contours of streamwise component')
end

