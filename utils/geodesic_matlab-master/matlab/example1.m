% this simple exampe shows the general principles of geodesic toolbox
% Danil Kirsanov, 09/2007 

global geodesic_library;                
geodesic_library = 'geodesic_debug';      %"release" is faster and "debug" does additional checks
%geodesic_library = 'geodesic_matlab_api';
rand('state', 0);                         %comment this statement if you want to produce random mesh every time



[vertices,faces] = read_ply('C:\Users\mraiola\Desktop\RSGs_4D_SPL_t_8.0000.ply');
plotmesh(vertices,faces)

mesh = geodesic_new_mesh(vertices,faces);         %initilize new mesh
algorithm = geodesic_new_algorithm(mesh, 'exact');      %initialize new geodesic algorithm

vertex_id = find(vertices(:,1:3)== cursors(1).Position,1);                                %create a single source at vertex #1
source_points = {geodesic_create_surface_point('vertex',vertex_id,vertices(vertex_id,:))};

geodesic_propagate(algorithm, source_points);   %propagation stage of the algorithm (the most time-consuming)

vertex_id = find(vertices(:,1:3)== cursors(2).Position,1);                                %create a single source at vertex #1
destination = geodesic_create_surface_point('vertex',vertex_id,vertices(vertex_id,:));
path = geodesic_trace_back(algorithm, destination);     %find a shortest path from source to destination

distances = zeros(N,1);              %find distances to all vertices of the mesh (actual pathes are not computed)

[source_id, distances] = geodesic_distance_and_source(algorithm);     %find distances to all vertices of the mesh; in this example we have a single source, so source_id is always equal to 1

geodesic_delete;                            %delete all meshes and algorithms

%-----------------plotting------------------------
hold off;
colormap('default');
trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3),distances, 'FaceColor', 'interp', 'EdgeColor', 'k');       %plot the mesh
daspect([1 1 1]);

hold on;
plot3(source_points{1}.x, source_points{1}.y, source_points{1}.z, 'or', 'MarkerSize',3);    %plot sources

plot3(destination.x, destination.y, destination.z, 'ok', 'MarkerSize',3);       %plot destination 
[x,y,z] = extract_coordinates_from_path(path);                                  %prepare path data for plotting
h = plot3(x*1.001,y*1.001,z*1.001,'k-','LineWidth',2);    %plot path
legend(h,'geodesic curve');