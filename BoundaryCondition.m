function U = BoundaryCondition(nx,ny,U,Northtype,Easttype,Southtype,Westtype,k)
if isequal(Southtype,'wall')
    U(1:2,:,1) = U(4:-1:3,:,1);
    U(1:2,:,2) = U(4:-1:3,:,2);
    U(1:2,:,3) = -U(4:-1:3,:,3);
    U(1:2,:,4) = U(4:-1:3,:,4);
elseif isequal(Southtype, 'periodic')
    U(1:2,:,:) = U(ny:ny+1,:,:);
else
    if k == 1
        disp([Southtype,' boundary condition is not available yet!'])
        disp('wall boundary condition is imposed on south boundary')
    end
    U(1:2,:,1) = U(4:-1:3,:,1);
    U(1:2,:,2) = U(4:-1:3,:,2);
    U(1:2,:,3) = -U(4:-1:3,:,3);
    U(1:2,:,4) = U(4:-1:3,:,4);
end
if isequal(Westtype,'wall')
    U(:,1:2,1) = U(:,4:-1:3,1);
    U(:,1:2,2) = -U(:,4:-1:3,2);
    U(:,1:2,3) = U(:,4:-1:3,3);
    U(:,1:2,4) = U(:,4:-1:3,4);
elseif isequal(Westtype, 'periodic')
    U(:,1:2,:) = U(:,nx:nx+1,:);
else
    if k == 1
        disp([Westtype,' boundary condition is not available yet!'])
        disp('wall boundary condition is alternatively imposed on west boundary')
    end
    U(:,1:2,1) = U(:,4:-1:3,1);
    U(:,1:2,2) = -U(:,4:-1:3,2);
    U(:,1:2,3) = U(:,4:-1:3,3);
    U(:,1:2,4) = U(:,4:-1:3,4);
end
if isequal(Northtype,'wall')
    U(ny+3:ny+4,:,1) = U(ny+2:-1:ny+1,:,1);
    U(ny+3:ny+4,:,2) = U(ny+2:-1:ny+1,:,2);
    U(ny+3:ny+4,:,3) = -U(ny+2:-1:ny+1,:,3);
    U(ny+3:ny+4,:,4) = U(ny+2:-1:ny+1,:,4);
elseif isequal(Northtype, 'periodic')
    U(ny+3:ny+4,:,:) = U(4:5,:,:);
else
    if k == 1
        disp([Northtype,' boundary condition is not available yet!'])
        disp('wall boundary condition is alternatively imposed on north boundary')
    end
    U(ny+3:ny+4,:,1) = U(ny+2:-1:ny+1,:,1);
    U(ny+3:ny+4,:,2) = U(ny+2:-1:ny+1,:,2);
    U(ny+3:ny+4,:,3) = -U(ny+2:-1:ny+1,:,3);
    U(ny+3:ny+4,:,4) = U(ny+2:-1:ny+1,:,4);
end
if isequal(Easttype,'wall')
    U(:,nx+3:nx+4,1) = U(:,nx+2:-1:nx+1,1);
    U(:,nx+3:nx+4,2) = -U(:,nx+2:-1:nx+1,2);
    U(:,nx+3:nx+4,3) = U(:,nx+2:-1:nx+1,3);
    U(:,nx+3:nx+4,4) = U(:,nx+2:-1:nx+1,4);
elseif isequal(Easttype, 'periodic')
    U(:,nx+3:nx+4,:) = U(:,4:5,:);
else
    if k == 1
        disp([Easttype,' boundary condition is not available yet!'])
        disp('wall boundary condition is alternatively imposed on east boundary')
    end
    U(:,nx+3:nx+4,1) = U(:,nx+2:-1:nx+1,1);
    U(:,nx+3:nx+4,2) = -U(:,nx+2:-1:nx+1,2);
    U(:,nx+3:nx+4,3) = U(:,nx+2:-1:nx+1,3);
    U(:,nx+3:nx+4,4) = U(:,nx+2:-1:nx+1,4);
end