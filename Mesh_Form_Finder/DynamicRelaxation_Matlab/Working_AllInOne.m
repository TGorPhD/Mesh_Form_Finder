clear all
close all

% Main Script
Filename = 'Planes_5';

% Load Geometry and Connectivity
[verts, edges, fixed, Xdata_All, Ydata_All, Zdata_All] = load_geometry(Filename);

% Parameters and Setup
params = struct( ...
    'MaximumIterations', 150, ...
    'tolerance', 5e-1, ...
    'dt', 0.49, ... % Time interval
    'm', 1, ...       % Fictitious mass
    'E_elastic', 1000, ... % Elastic modulus (MPa)
    'Area', (0.1^2) * pi, ... % Cross-sectional area
    'constantForce', [0, 0, 0], ... % External forces
    'meshType', 'triangular_prismatic', ... % Mesh type
    'length_threshold', 0.5 ... % Threshold for long edges
);

nodeCount = 10;
nodelength = 100;
meshType= 'quad'; %tri, quad, tri_both
fixedNodeOption = 'all'; %'all', 'everyOther'%everother doesnt work as the node indexing is off over columns
[verts, edges, fixed, fixed_indices, midlines, diagonals] = generate_mesh_Square(nodeCount, nodelength, meshType, fixedNodeOption);

edgeVectors = verts(edges(:,2),:) - verts(edges(:,1),:);
restlengthOG = vecnorm(edgeVectors,2,2);
params.restlengthOG = restlengthOG;

% Apply sinusoidal shape to the diagonals
maxHeight = 50;
plotOption = true;
LinesChanged = diagonals; %diagonals.main, midlines.horizontal, struct('main', midlines.horizontal) # speciic ones: struct('main', [70;71;210;211])
shapeFunc = 'gaussian'; %'gaussian', 'pyramid' , 'sinusoidal', 'userdefined'
[verts,NewFixedIndices]  = apply_shape_function(verts, LinesChanged, shapeFunc, maxHeight, plotOption);
fixed = [fixed;verts(NewFixedIndices,:)];

%%
%{
%updating fixed for every other?
fixed = verts(fixed_indices(mod(fixed_indices, 2) == 1),:);

%deciding control indices
control_indices = [53:58];
verts(control_indices,3) = [10,18,25,25,18,10]; %applying random disp to one node to height of 40
fixed = [fixed;verts(control_indices,:)];
%}
[xyz, vertexCount, freeIndices] = initialize_geometry(verts, fixed);

AppliedLinkForce = zeros(vertexCount, 3); % External forces on nodes

%try random node force: 
%indtry = [1,15,33,34,35,36,200];
%AppliedLinkForce([indtry],:) = repmat([300,300,0],[size(indtry,2),1]);


% Initial Visualization
plot_geometry(xyz, edges, fixed);

% Dynamic Relaxation with Animation
[LinkTension, xyz, KEAllStore, edges, history] = dynamic_relaxation(params, ...
    xyz, edges, freeIndices, fixed, AppliedLinkForce);

% Post-processing and Visualization
plot_geometry(xyz, edges, fixed); % Final deformed geometry
plot_kinetic_energy(KEAllStore);
plot_final_tensions(edges, LinkTension);

% Create Animation and GIF
filename = 'dynamic_relaxation';
%create_animation(history, edges,filename, fixed);
%animate_plot(history, edges, fixed);


%% Function Definitions

function [updatedVerts,NewFixedIndices] = apply_shape_function(verts, lines, shapeFunction, maxHeight, plotOption)
    % Applies a shape function to the vertices along specified lines.
    %
    % Inputs:
    %   verts - Vertex coordinates (n x 3 matrix)
    %   lines - Struct with diagonal indices or other line definitions
    %   shapeFunction - String specifying the shape function:
    %       'pyramid', 'gaussian', 'sinusoidal', etc.
    %   maxHeight - Maximum height for the shape function
    %   plotOption - Boolean, if true, plot before and after
    %
    % Outputs:
    %   updatedVerts - Modified vertex coordinates (n x 3 matrix)

    % Validate shape function
    validShapes = {'pyramid', 'gaussian', 'sinusoidal','userdefined'};
    if ~ismember(shapeFunction, validShapes)
        error('Invalid shape function. Choose from "pyramid", "gaussian", "sinusoidal", "userdefined".');
    end

    % Combine all lines into a single list of unique indices
    %lineIndices = unique([lines.main; lines.anti]);
    %dynamicedit:
    % Assuming lines is a struct with various fields (e.g., 'main', 'anti', ...)
    fieldNames = fieldnames(lines);  % Get all field names in the 'lines' struct
    
    % Concatenate all values from the fields
    lineIndices = [];
    for i = 1:length(fieldNames)
        lineIndices = [lineIndices; lines.(fieldNames{i})];  % Concatenate each field's values
    end
    
    % Get unique indices
    lineIndices = unique(lineIndices);

    NewFixedIndices = lineIndices;

    % Extract vertices corresponding to the line indices
    lineVerts = verts(lineIndices, :);

    % Calculate distances from the line's center for each vertex
    lineCenter = mean(lineVerts(:, 1:2), 1); % Center in the x-y plane
    distances = vecnorm(lineVerts(:, 1:2) - lineCenter, 2, 2);

    % Normalize distances to [0, 1]
    maxDistance = max(distances);
    normalizedDistances = distances / maxDistance;

    % Compute new z values based on the selected shape function
    switch shapeFunction
        case 'pyramid'
            % Linear height distribution (inverted V shape)
            newZ = maxHeight * (1 - normalizedDistances);

        case 'gaussian'
            % Gaussian curve
            newZ = maxHeight * exp(-4 * (normalizedDistances).^2);

        case 'sinusoidal'
            % Half sinusoidal wave (arch-like shape)
            %newZ = maxHeight * (sin(pi * normalizedDistances / 2).^2);
            newZ = maxHeight * sin(pi * (1 - normalizedDistances) / 2).^2;

        case 'userdefined'
            %just apply heigth to nodes?
            newZ = maxHeight;

        otherwise
            error('Unsupported shape function.');
    end

    % Assign new z values to the vertices
    updatedVerts = verts;
    updatedVerts(lineIndices, 3) = newZ;

    % Plot before and after if plotOption is enabled
    if plotOption
        figure;

        % Original vertices
        subplot(1, 2, 1);
        scatter3(verts(:, 1), verts(:, 2), verts(:, 3), 20, 'filled');
        hold on;
        scatter3(lineVerts(:, 1), lineVerts(:, 2), lineVerts(:, 3), 50, 'r', 'filled');
        title('Original Mesh');
        xlabel('X'); ylabel('Y'); zlabel('Z');
        grid on; view(3);

        % Updated vertices
        subplot(1, 2, 2);
        scatter3(updatedVerts(:, 1), updatedVerts(:, 2), updatedVerts(:, 3), 20, 'filled');
        hold on;
        scatter3(updatedVerts(lineIndices, 1), updatedVerts(lineIndices, 2), ...
                 updatedVerts(lineIndices, 3), 50, 'r', 'filled');
        title(['Mesh After Applying ', shapeFunction, ' Shape']);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        grid on; view(3);
        axis equal;
    end
end



function [verts, edges, fixed, Xdata_All, Ydata_All, Zdata_All] = load_geometry(Filename)
    % Load geometry and connectivity data
    load(Filename);
    Xdata_All = DatStruct.Xdata;
    Ydata_All = DatStruct.Ydata;
    Zdata_All = DatStruct.Zdata;
    s1 = DatStruct.Connectivity(:, 1);
    s2 = DatStruct.Connectivity(:, 2);
    
    verts = [Xdata_All, Ydata_All, Zdata_All];
    fixed = verts((verts(:, 1) == Xdata_All(end)) | (verts(:, 1) == Xdata_All(1)), :);
    edges = [s1, s2];
end

function [xyz, vertexCount, freeIndices] = initialize_geometry(verts, fixed)
    % Initialize geometry
    xyz = verts;
    vertexCount = size(verts, 1);
    freeIndices = find(~ismember(verts, fixed, 'rows'));
end

function [verts, edges, fixed, fixed_indices, midlines, diagonals] = generate_mesh_Square(nodeCount, length, meshType, fixedNodeOption)
    % Generate a square mesh with specified properties
    %
    % Inputs:
    %   nodeCount - Number of nodes per side (scalar, >= 2)
    %   length - Total length of the square side (scalar, > 0)
    %   meshType - 'quad' for quadrilateral or 'tri' for triangular
    %   fixedNodeOption - Control for fixed nodes:
    %       'all'  - All edge nodes fixed
    %       'everyOther' - Every other edge node fixed
    %       'none' - No nodes fixed
    %
    % Outputs:
    %   verts - Nodal coordinates (n x 3 matrix)
    %   edges - Connectivity matrix (m x 2)
    %   fixed - Fixed nodes (p x 3 matrix)
    %   fixed_indices - Indices of fixed nodes
    %   midlines - Struct with indices for horizontal and vertical midlines
    %   diagonals - Struct with indices for main and anti-diagonals

    % Validate inputs
    if nodeCount < 2
        error('nodeCount must be at least 2.');
    end
    if length <= 0
        error('length must be greater than 0.');
    end
    if ~ismember(meshType, {'quad', 'tri', 'tri_both'})
        error('meshType must be "quad", "triboth" or "tri".');
    end
    if ~ismember(fixedNodeOption, {'all', 'everyOther', 'none'})
        error('fixedNodeOption must be "all", "everyOther", or "none".');
    end

    % Step 1: Create nodes
    spacing = length / (nodeCount - 1);
    [x, y] = meshgrid(0:spacing:length, 0:spacing:length);
    verts = [x(:), y(:), zeros(numel(x), 1)]; % All nodes lie in the z = 0 plane initially

    % Step 2: Define connectivity
    edges = [];
    for i = 1:nodeCount-1
        for j = 1:nodeCount-1
            % Get indices of the current square's corners
            n1 = sub2ind([nodeCount, nodeCount], i, j);
            n2 = sub2ind([nodeCount, nodeCount], i, j+1);
            n3 = sub2ind([nodeCount, nodeCount], i+1, j);
            n4 = sub2ind([nodeCount, nodeCount], i+1, j+1);

            if strcmp(meshType, 'quad')
                % Quadrilateral connectivity
                edges = [edges; n1, n2; n2, n4; n4, n3; n3, n1];
            elseif strcmp(meshType, 'tri')
                % Triangular connectivity
                edges = [edges; n1, n2; n2, n4; n4, n1; n1, n3; n3, n4];
            elseif strcmp(meshType, 'tri_both')
                % Triangular connectivity (both diagonals per quad)
                % First diagonal (n1-n4) 
                edges = [edges; n1, n2; n2, n4; n4, n1];  % Triangle 1
                % Second diagonal (n2-n3)
                edges = [edges; n2, n3; n3, n4; n4, n2];  % Triangle 2
            end
        end
    end

    % Remove duplicate edges
    edges = unique(sort(edges, 2), 'rows');

    % Step 3: Determine fixed nodes
    edgeIndices = unique([1:nodeCount, ... % Bottom edge
                          nodeCount:nodeCount:(nodeCount^2), ... % Right edge
                          (nodeCount^2):-1:(nodeCount^2-nodeCount+1), ... % Top edge
                          (nodeCount^2-nodeCount+1):-nodeCount:1]); % Left edge

    if strcmp(fixedNodeOption, 'all')
        fixed_indices = edgeIndices;
    elseif strcmp(fixedNodeOption, 'everyOther')
        fixed_indices = edgeIndices(1:2:end); % Select every other node
        %this is how get first row: test1 = verts((verts(:,1) == 0),:)
        %this is how get last row: test2 = verts((verts(:,1) == Length),:)
        %get every other element: test1(1:2:end), test2(1:2:end)
    else
        fixed_indices = [];
    end
    fixed = verts(fixed_indices, :);

    % Step 4: Identify midlines and diagonals
    midlines.horizontal = find(abs(verts(:, 2) - (length / 2)) < spacing / 2);
    midlines.vertical = find(abs(verts(:, 1) - (length / 2)) < spacing / 2);

    diagonals.main = find(abs(verts(:, 1) - verts(:, 2)) < spacing / 2);
    diagonals.anti = find(abs(verts(:, 1) + verts(:, 2) - length) < spacing / 2);
end

function plot_geometry(xyz, edges, fixed)
    % Plot geometry in 3D
    figure;
    hold on;
    for i = 1:size(edges, 1)
        plot3([xyz(edges(i, 1), 1), xyz(edges(i, 2), 1)], ...
              [xyz(edges(i, 1), 2), xyz(edges(i, 2), 2)], ...
              [xyz(edges(i, 1), 3), xyz(edges(i, 2), 3)], 'b');
    end
    scatter3(fixed(:, 1), fixed(:, 2), fixed(:, 3), 'ro', 'filled');
    scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'k');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Geometry');
    axis equal;
    grid on;
end

function h = plot_geometry2(xyz, edges, fixed)
    % Plot geometry in 3D and return the plot handle
    
    hold on;
    h = gca;  % Get the current axes handle
    for i = 1:size(edges, 1)
        plot3([xyz(edges(i, 1), 1), xyz(edges(i, 2), 1)], ...
              [xyz(edges(i, 1), 2), xyz(edges(i, 2), 2)], ...
              [xyz(edges(i, 1), 3), xyz(edges(i, 2), 3)], 'b');
    end
    scatter3(fixed(:, 1), fixed(:, 2), fixed(:, 3), 'ro', 'filled');
    scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'k');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Geometry');
    axis equal;
    grid on;
end


function [LongEdge] = identify_longitudinal_edges(edges, xyz, params)
    % Identify longitudinal edges dynamically
    switch params.meshType
        case 'triangular_prismatic'
            % For prismatic meshes, edges with dominant Z-direction components
            vectors = xyz(edges(:, 2), :) - xyz(edges(:, 1), :);
            zDominance = abs(vectors(:, 3)) ./ vecnorm(vectors, 2, 2);
            LongEdge = zDominance > 0.8; % Threshold for "vertical" edges
        case 'cable_net'
            % Handle different logic for cable nets if needed
            LongEdge = false(size(edges, 1), 1); % Placeholder
        otherwise
            error('Unsupported mesh type: %s', params.meshType);
    end
end

function [LinkTension, xyz, KEAllStore, edges, history] = dynamic_relaxation(params, ...
        xyz, edges, freeIndices, fixed, AppliedLinkForce)

    %% Extract parameters
    dt = params.dt;
    %K_S = params.E_elastic;
    Area = params.Area;
    restlengthOG = params.restlengthOG;
    tolerance = params.tolerance;
    MaximumIterations = params.MaximumIterations;

    %% Initialize variables
    vertexCount = size(xyz, 1);
    edgeCount = size(edges, 1);
    
    %edgeVectors = xyz(edges(:,2),:) - xyz(edges(:,1),:);
    %restlengthOG = vecnorm(edgeVectors,2,2);
    K_S = (params.E_elastic * params.Area) ./ restlengthOG;

    v = zeros(vertexCount, 3); % Initializing velocities
    S = zeros(vertexCount, 3); % Forces on nodes
    R = zeros(vertexCount, 3); % Residual forces
    KE = 0; % Initial kinetic energy
    ResetFlag = 0;

    KEAllStore = {}; % Store kinetic energy for all iterations
    history = struct; % History to track values per iteration
    history.xyz = {};
    history.v = {};
    history.KE = {};
    history.m = {};

    %% Main loop
    t = 1;
    difference = inf;
    KEStore = [];

    while (difference > tolerance) && (t <= MaximumIterations)
        %% Step 1: Initialization
        vp = v; % Store previous velocities
        xyz0 = xyz;  % Store previous positions
        KE_0 = KE;  % Store previous kinetic energy

        %% Step 2: Compute Link Lengths and Forces
        LinkVectors_t = xyz(edges(:,2),:) - xyz(edges(:,1),:);
        LinkLengths_t = vecnorm(LinkVectors_t, 2, 2);

        LinkTension = repmat(params.constantForce(2), [length(LinkLengths_t), 1]);
        %LongEdge = ismember(edges,[LongMem;fliplr(LongMem)],"rows");
        %LinkTension(:) = constantForce(3);
        % Compute tensions
        LinkTension = LinkTension + K_S .* (LinkLengths_t - restlengthOG);
        %LinkTension = K_S .* (LinkLengths_t - restlengthOG);
        LinkTension(LinkTension < 0) = 0;

        %% Step 3: Node-Specific Mass Computation
        m = zeros(vertexCount, 1);
        GeomStiffness = (LinkTension ./ LinkLengths_t);
        for nodes = 1:vertexCount
            AdjoiningEdges = sum(edges == nodes, 2);
            stiffness = sum((K_S ./ restlengthOG) .* AdjoiningEdges + GeomStiffness .* AdjoiningEdges);
            m(nodes) = stiffness;
        end

        %% Step 4/5: Force Resolution
        EdgeForceS = zeros(edgeCount, 3);
        for i = 1:edgeCount
            for j = 1:3
                EdgeForceS(i,j) = LinkVectors_t(i,j) / LinkLengths_t(i) * LinkTension(i);
            end
            S(edges(i,1),:) = S(edges(i,1),:) + EdgeForceS(i,:);
            S(edges(i,2),:) = S(edges(i,2),:) - EdgeForceS(i,:);
        end

        R(freeIndices,:) = AppliedLinkForce(freeIndices,:) + S(freeIndices,:);

        %% Step 6: Update Velocities and Positions
        A = 1; % Damping constant
        B = repmat((dt ./ m), [1, 3]);

        v(freeIndices,:) = A * vp(freeIndices,:) + B(freeIndices) .* R(freeIndices,:);
        xyz(freeIndices,:) = xyz0(freeIndices,:) + dt * v(freeIndices,:);

        KE = sum(0.5 * m .* vecnorm(v, 2, 2).^2);
        KEStore(t) = KE;

        %% Step 7/8: Kinetic Energy Correction
        if t <= 2 || ResetFlag == 1
            ResetFlag = 0;
        else
            if KE > KE_0 && ResetFlag == 0
                % Backtracking correction
                E = KEStore(t-1) - KEStore(t);
                D = KEStore(t-2) - KEStore(t-1);
                q = max(0, min(E / (E - D), 1)); % Clamp q between 0 and 1
                
                R_t = (m .* (v - vStore{t-1})) / dt;
                xyz = xyz - (dt * (1 + q) * vStore{t-1}) + (dt^2 / 2) * q * (R_t ./ m);

                % Reset variables
                KE = 0;
                v = zeros(vertexCount, 3);
                R = zeros(vertexCount, 3);
                S = zeros(vertexCount, 3);
                ResetFlag = 1;
                warning(['Timestep ', num2str(t), ': KE spike, resetting process']);
            end
        end

        %% Step 9: Convergence Check
        difference = norm(xyz - xyz0, 'fro');
        vStore{t} = v; % Store velocities for backtracking
        history.xyz{t} = xyz; % Store position history
        history.v{t} = v; % Store velocity history
        history.KE{t} = KE; % Store kinetic energy history
        history.m{t} = m;
        t = t + 1;
    end
    %store how many iterations:
    history.iters = t-1;
    % Store final kinetic energy history
    KEAllStore = KEStore;

end

function plot_kinetic_energy(KEAllStore)
    figure;
    plot(KEAllStore, 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('Kinetic Energy');
    title('Kinetic Energy Convergence');
end

function plot_final_tensions(edges, LinkTension)
    figure;
    plot(1:length(LinkTension), LinkTension, 'o');
    xlabel('Edge Index');
    ylabel('Tension (N)');
    title('Final Link Tensions');
end

function create_animation(history, edges, filename, fixed)
    filename = [filename, '.gif'];
    for t = 1:length(history)
        plot_geometry(history{t}, edges, fixed);
        drawnow;
        
        view([130 15])

        % Capture the current frame
        frame = getframe(gcf);
        img = frame2im(frame);
        [imind, cm] = rgb2ind(img, 256);
        
        % Write to GIF (create or append)
        if t == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
        
        % Close the figure after saving the frame to avoid opening many windows
        close(gcf);
    end
end

function animate_plot(history, edges, fixed)
    % Animate the plot in real-time without saving it as a GIF, on a single plot
    figure;
    % Initial plot setup
    %h = plot_geometry2(history{1}, edges, fixed); % Plot the initial geometry
    h = plot_geometry2(history.xyz{1}, edges, fixed); % Correct dot and cell indexing

    
    % Loop over the history and update the plot
    for t = 1:history.iters
        % Clear the current axes
        cla;
        
        % Plot the geometry at the current time step
        h = plot_geometry2(history.xyz{t}, edges, fixed);
        
        % Update plot title or other properties if needed
        title(['Iteration: ', num2str(t)]);
        
        % Redraw and pause to create animation effect
        drawnow;
        view([130 15])
        %pause(0.001); % Adjust this for desired animation speed
    end
end

