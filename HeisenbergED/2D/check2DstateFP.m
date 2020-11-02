function Ns = check2DstateFP(s,mx,my,Nx,Ny,numtype)
   Fx = 1;
   Fy = 1;
   class = zeros(1,Nx*Ny-1);
   
   for x = 0:Nx-1
       for y = 0:Ny-1
           t = s;
           for a = 1:x
              t = Tx(t,Nx,Ny,numtype);
              if t < s
                 Ns = -1; return; 
              end
           end
           for b = 1:y
              t = Ty(t,Nx,Ny,numtype);
              if t < s
                 Ns = -1; return; 
              end
           end
           i = x + y*Nx;
           class(i+1) = t;
           if t == s
              if x ~= 0 
                  phaseX = -2*pi*mx*x/Nx;
                  Fx = Fx + exp(-1i*phaseX);
              end
              if y ~= 0 
                  phaseY = -2*pi*my*y/Ny;
                  Fy = Fy + exp(-1i*phaseY);
              end
           end
       end
   end
   F = norm(Fx*Fy);
   class = unique(class);
   D = length(class);
   Ns = D*F^2;
end