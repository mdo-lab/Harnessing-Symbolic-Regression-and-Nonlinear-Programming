% simulation.m
function [md, mass, h] = simulation(d1, d2, d3, mat, doPlot)
    if nargin < 5, doPlot = false; end
    warning('off','all');

    % Material properties
    switch mat
        case 1 % Aluminum
            d  = 2700;   ym = 70e9;   pr = 0.33;
        case 2 % Steel
            d  = 7800;   ym = 205e9;  pr = 0.29;
        case 3 % Nylon
            d  = 1150;   ym = 2e9;    pr = 0.40;
        otherwise
            error('mat must be 1 (Al), 2 (Steel), or 3 (Nylon).');
    end

    % Base curves
    t = linspace(0,2*pi,500)';
    xin = 15*cos(t);  yin = 9*sin(t);

    % Trim inner ellipse to |x|<11
    keep = xin < 11 & xin > -11;
    x1 = xin(keep);  y1 = yin(keep);

    % Other arcs
    x2 = 7 + d1*cos(t);   y2 = d2*sin(t);
    x3 = d1*cos(t) - 7;   y3 = d2*sin(t);
    x4 = d3*cos(t);       y4 = 6 + d3*sin(t);
    x5 = d3*cos(t);       y5 = d3*sin(t) - 6;
    x6 = d3*cos(t);       y6 = d3*sin(t);

    % Build polygon & mesh
    pgon = polyshape({x1,x2,x3,x4,x5,x6},{y1,y2,y3,y4,y5,y6});
    tr = triangulation(pgon);
    model = createpde('structural','static-planestress');
    geometryFromMesh(model,tr.Points',tr.ConnectivityList');
    generateMesh(model,'Hmax',0.25);

    % Physics
    structuralProperties(model,'YoungsModulus',ym,'PoissonsRatio',pr);
    structuralBC(model,'Edge',3,'Constraint','fixed');
    structuralBoundaryLoad(model,'Edge',2,'SurfaceTraction',[0;2000]);

    % Solve
    R = solve(model);

    % Outputs
    mass = area(pgon)*d;
    md   = max(R.Displacement.Magnitude);

        % ---- Plotting (avoid tiledlayout to keep pdeplot happy) ----
    h = struct();
    if doPlot
        h.fig = figure('Color','w','Name','Simulation');

        % Left axes: geometry & mesh
        h.ax1 = axes('Parent', h.fig, 'Position', [0.06 0.15 0.40 0.75]);
        axes(h.ax1); %#ok<LAXES>
        hold on;
        plot(pgon,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0.25 0.25 0.25]);
        pdemesh(model.Mesh);
        axis equal; box on;
        title('Geometry & Mesh');
        xlabel('x'); ylabel('y');
        hold off;

        % Right axes: displacement magnitude (let pdeplot draw into this axes)
        h.ax2 = axes('Parent', h.fig, 'Position', [0.56 0.15 0.38 0.75]);
        axes(h.ax2); %#ok<LAXES>
        pdeplot(model, ...
            'XYData', R.Displacement.Magnitude, ...
            'Deformation', R.Displacement, ...
            'DeformationScaleFactor', 1, ...
            'Mesh','off');
        axis equal; box on;
        colorbar(h.ax2);
        title(sprintf('Displacement Magnitude (max = %.3g m)', md));
        xlabel('x'); ylabel('y');

        % Footer text
        h.massText = annotation('textbox',[0.02 0.02 0.96 0.05], ...
            'String', sprintf('Mass = %.3f kg (density × area)', mass/1000), ...
            'EdgeColor','none', 'HorizontalAlignment','center');
    end
end