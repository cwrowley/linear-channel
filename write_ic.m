function [init_vec]=write_ic(Ny,Nz,z,y,ic)

Npts=Nz*(Ny-1)

if (ic ==2 )
   
    os_data = load('/home/milak/balpod/matlab/OS_eigenfuns/os16a1b0R100.dat');
    
end


    %write delta function IC
    icwidth=1
    y_mid=0
    z_mid=z(Nz/2)
    for k=1:(Ny-1)
        for m=1:Nz
            if (ic ==1) 
               init_cond(k,m)=0.5*(cos(pi*y(k+1))+1)*exp(-2*( (y(k+1)-y_mid)^2 + (z(m)-z_mid)^2) ...
                /icwidth^2);
              %  init_cond(k,m)=(cos(pi*y(k+1))+1)*sin(z(m));
            else
             %  init_cond(k,m)=(os_data((k+1),4)+sqrt(-1)*os_data((k+1),5))*exp(sqrt(-1)*2*z(m)); 
            init_cond(k,m)=(os_data((k+1),4)+sqrt(-1)*os_data((k+1),5));%*exp(sqrt(-1)*(m)); 
            end
                
        end
    end

    %re-write as vector
    for k=1:(Ny-1)
        init_vec((k-1)*Nz+1:k*Nz) = init_cond(k,:);
    end
    %for m=1:Nz
    %    init_vec((m-1)*(Ny-1)+1:m*(Ny-1)) = init_cond(k,:);
    %end

%Bmat
    %plot initial condition to test
    % init_cond=[zeros(1, Nz) ; init_cond' ; zeros(1, Nz)];
    init_cond=[zeros(Nz,1)  init_cond'  zeros(Nz,1)];

    %figure(2)
    %surf(y,z,init_cond)

 %   if adj==0
        init_vec=[init_vec'; zeros(Npts,1)];
 %   else
        %the next line gives IC to vorticity only

   %     init_vec=[zeros(Npts,1); init_vec'];
%    end


%end

