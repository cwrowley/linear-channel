classdef LinearChannel < handle
    % LINEARCHANNEL  Streamwise-constant linear channel flow
    %   Detailed explanation goes here
    
    properties
        nz = 0;  % number of Fourier modes in spanwise direction
        ny = 0;  % degree of Chebyshev polynomial in wall-normal direction
        npts = 0; % total number of grid points = (ny-1) * nz
        Re = 0;  % Reynolds number
        A = 0; % dynamics given by x' = A x
        M = 0;
        y = 0; % collocation points in wall-normal direction
        z = 0; % collocation points in spanwise direction
    end
    
    methods
        function obj = LinearChannel(ny, nz, Re)
            % LINEARCHANNEL  Create an object for streamwise-constant channel flow
            % linchan = LINEARCHANNEL(ny, nz, Re)
            %
            % ny: number of Chebyshev modes in y-direction
            % nz: number of Fourier modes in z-direction
            % Re: Reynolds number
            obj.ny = ny;
            obj.nz = nz;
            obj.npts = (ny-1) * nz;
            obj.Re = Re;
            obj.define_equations()
        end
        
        function define_equations(obj)
            % Chebyshev modes in y-dir
            Ny = obj.ny;
            Nz = obj.nz;
            Npts = obj.npts;
            [Dy11,obj.y] = cheb(Ny);
            Dy1 = Dy11(2:Ny,2:Ny);
            %Dyy1 = Dy11^2;
            %Dyy1 = Dyy1(2:Ny,2:Ny); % homogeneous boundary conditions
            %
            % fourth-order boundary conditions
            S = diag([0; 1./(1-obj.y(2:Ny).^2); 0]);
            D4y = (diag(1-obj.y.^2)*Dy11^4 - 8*diag(obj.y)*Dy11^3 - 12*Dy11^2)*S;
            D4y1=D4y(2:Ny,2:Ny);
            %clamped BC's on v
            D2y = (diag(1-obj.y.^2)*Dy11^2 - 4*diag(obj.y)*Dy11 - 2*eye(Ny+1))*S;
            Dyy1 = D2y(2:Ny,2:Ny);
            
            % stack into matrix
            Dy = zeros(Npts);
            Dyy = zeros(Npts);
            for i=1:Ny-1
                for j=1:Ny-1
                    Dy((i-1)*Nz+1:i*Nz, (j-1)*Nz+1:j*Nz) = diag(Dy1(i,j)*ones(Nz,1));
                    Dyy((i-1)*Nz+1:i*Nz, (j-1)*Nz+1:j*Nz) = diag(Dyy1(i,j)*ones(Nz,1));
                end
            end
            
            % stack into matrix
            D4y = zeros(Npts);
            %Dyy = zeros(Npts);
            for i=1:Ny-1
                for j=1:Ny-1
                    D4y((i-1)*Nz+1:i*Nz, (j-1)*Nz+1:j*Nz) = diag(D4y1(i,j)*ones(Nz,1));
                    %        Dyy((i-1)*Nz+1:i*Nz, (j-1)*Nz+1:j*Nz) = diag(Dyy1(i,j)*ones(Nz,1));
                end
            end
            
            % Fourier modes in z-dir
            [Dz1,obj.z] = fourier(Nz);
            Dzz1 = Dz1^2;
            
            % stack into matrix
            Dz = zeros(Npts);
            Dzz = zeros(Npts);
            for i=1:Ny-1
                Dz((i-1)*Nz+1:i*Nz, (i-1)*Nz+1:i*Nz) = Dz1;
                Dzz((i-1)*Nz+1:i*Nz, (i-1)*Nz+1:i*Nz) = Dzz1;
            end
            
            % Define grid
            zvar = zeros(Npts,1); yvar = zeros(Npts,1);
            for j=1:Ny-1
                yvar((j-1)*Nz+1:j*Nz) = obj.y(j+1);
                zvar((j-1)*Nz+1:j*Nz) = obj.z;
            end
            
            % Mean flow: U(y) = 1 - y^2
            Uprime = -2*yvar;            
            
            obj.M = [-(Dyy+Dzz), zeros(Npts); zeros(Npts), eye(Npts)];
            % the A matrix for streamwise-constant perturbations
            obj.A = obj.M \ [-1/obj.Re* (Dzz^2 + 2*Dyy*Dzz + D4y), zeros(Npts);
                -diag(Uprime) * Dz, 1/obj.Re * (Dyy + Dzz)];
            
            %obj.M = (obj.M + obj.M')/2;
        end
        
        function vec = optimal_perturbation(obj, t)
            % OPTIMAL_PERTURBATION  Initial perturbation that maximizes energy growth
            %
            %   vec = chan.OPTIMAL_PERTURBATION(t)
            
            expA = expm(obj.A * t);
            [~,~,V] = svd(expA);
            
            % return a linear combination of the first two singular vectors so that
            % vorticity = zero at z = 0 or 2*pi
            vec1 = V(:,1);
            vec2 = V(:,2);
            
            % index corresponding to vorticity at a point with z = 2*pi
            jmid = obj.ny/2;
            ind = obj.npts + obj.nz * jmid;
            
            a = vec1(ind);
            b = vec2(ind);
            r = sqrt(a^2 + b^2);
            a = a / r;
            b = b / r;
            
            vec = b * vec1 - a * vec2;            
        end
        
        function [vel, vort] = vec2field(obj, vec)
            % VEC2FIELD  Convert a vector to velocity & vorticity fields
            %
            %  [vel, vort] = VEC2FIELD(vec)
            
            % Note: duplicate endpoints of periodic BCs for symmetric plots
            vel = zeros(obj.ny+1, obj.nz+1);
            vort = zeros(obj.ny+1, obj.nz+1);
            
            for j = 2:obj.ny
                vel_offset = obj.nz * (j-2);
                vel(j,2:obj.nz+1) = vec(vel_offset+1:vel_offset+obj.nz);
                vel(j,1) = vel(j,obj.nz+1);
                vort_offset = obj.npts + obj.nz * (j-2);
                vort(j,2:obj.nz+1) = vec(vort_offset+1:vort_offset+obj.nz);
                vort(j,1) = vort(j,obj.nz+1);
            end
        end
        
        function plotfield(obj, field, clev)
            if nargin < 3
                clev = 14;
            end
            % add a zero to beginning of z coordinates, for plotting
            zz = zeros(obj.nz+1, 1);
            zz(1) = 0;
            zz(2:end) = obj.z;
            contourf(zz, obj.y, field, clev)
            axis equal
            xlabel('z')
            ylabel('y')
        end
    end
    
end
