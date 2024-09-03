function nsrf = morphSurf3D(osrf, res, orig, spacing,im)
	F = mirt3D_F(res.okno);
	[Xx,Xy,Xz] = mirt3D_nodes2grid(res.X, F, res.okno);

	nx = size(Xx,2);
	ny = size(Xy,1);
	nz = size(Xz,3);
	[x, y, z]  = meshgrid(1:nx,1:ny,1:nz);


	dx = (Xx - x).*spacing.x;
	dy = (Xy - y).*spacing.y;
	dz = (Xz - z).*spacing.z;
	
    x = (x-1).*spacing.x + orig.x;
    y = (y-1).*spacing.y + orig.y;
    z = (z-1).*spacing.z + orig.z;
    
	disp(:,1) = interp3(x,y,z,dx,...
        osrf.vertices(:,1),osrf.vertices(:,2),osrf.vertices(:,3),'linear');
	disp(:,2) = interp3(x,y,z,dy,...
        osrf.vertices(:,1),osrf.vertices(:,2),osrf.vertices(:,3),'linear');
	disp(:,3) = interp3(x,y,z,dz,...
        osrf.vertices(:,1),osrf.vertices(:,2),osrf.vertices(:,3),'linear');
    

    if (nnz(isnan(disp))>0)
        fprintf('Warning: some part of surface lies outside of domain\n');
        disp(isnan(disp)) = 0;
    end
    
    nsrf.vertices = osrf.vertices -disp;
%     nsrf.faces = int32(osrf.faces);
% 
%         fid = fopen('temp_imdisp.vtk','w');
%         fprintf(fid, '# vtk DataFile Version 3.0\n');
%         fprintf(fid, '3D Displacement Data\n');
%         fprintf(fid, 'ASCII\n');
%         fprintf(fid, 'DATASET RECTILINEAR_GRID\n');
%         fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
% 
%         xc = reshape(x(1,:,1),[1, nx]);
%         yc = reshape(y(:,1,1),[1, ny]);
%         zc = reshape(z(1,1,:),[1, nz]);
%         fprintf(fid, 'X_COORDINATES %d double\n', nx);
%         for i=1:nx
%             fprintf(fid, '%f ', xc(i));
%         end
%         fprintf(fid, '\nY_COORDINATES %d double\n', ny);
%         for j=1:ny
%             fprintf(fid, '%f ', yc(j));
%         end
%         fprintf(fid, '\nZ_COORDINATES %d double\n', nz);
%         for k=1:nz
%             fprintf(fid, '%f ', zc(k));
%         end
%         
%         temp = zeros(ny,nx,nz);
%         temp(1:min(ny,size(im,1)), 1:min(nx,size(im,2)), 1:min(nz,size(im,3))) = ...
%             im(1:min(ny,size(im,1)), 1:min(nx,size(im,2)), 1:min(nz,size(im,3)));
%         fprintf(fid, '\nPOINT_DATA %d\n',nx*ny*nz);
%         fprintf(fid, 'SCALARS Intensity double\nLOOKUP_TABLE default\n');
%         temp(isnan(temp)) = 0;
%         for k=1:nz
%             for j=1:ny
%                 for i=1:nx
%                     fprintf(fid,'%f\n',temp(j,i,k));
%                 end
%             end
%         end
%         
%         fprintf(fid, 'VECTORS Displacement double\n');
%         for k=1:nz
%             for j=1:ny
%                 for i=1:nx
%                     fprintf(fid,'%f %f %f\n',...
%                         -dx(j,i,k),-dy(j,i,k),-dz(j,i,k));
%                 end
%             end
%         end
%         fclose(fid);
%         
%         nno = size(nsrf.vertices,1);
%         nel = size(nsrf.faces,1);
%         fid = fopen('temp_srfdisp.vtk','w');
%         fprintf(fid, '# vtk DataFile Version 3.0\n');
%         fprintf(fid, '3D Displacement Data\n');
%         fprintf(fid, 'ASCII\n');
%         fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
%         
%         fprintf(fid, 'POINTS %d double\n', nno);
%         for i=1:nno
%             for j=1:3
%                 fprintf(fid, '%f ',osrf.vertices(i,j));
%             end
%             fprintf(fid,'\n');
%         end
%         
%         fprintf(fid, 'CELLS %d %d\n', nel, nel*4);
%         for i=1:nel
%             fprintf(fid,'3 ');
%             for j=1:3
%                 fprintf(fid, '%d ', osrf.faces(i,j)-1);
%             end
%             fprintf(fid,'\n');
%         end
%         
%         fprintf(fid, 'CELL_TYPES %d\n', nel);
%         for i=1:nel
%             fprintf(fid, '5\n');
%         end
%         
%         fprintf(fid, '\nPOINT_DATA %d\n',nno);
%         fprintf(fid, 'VECTORS Displacement double\n');
%         for i=1:nno
%             for j=1:3
%                 fprintf(fid,'%f ',-disp(i,j));
%             end
%             fprintf(fid,'\n');
%         end
%         fclose(fid);
end
