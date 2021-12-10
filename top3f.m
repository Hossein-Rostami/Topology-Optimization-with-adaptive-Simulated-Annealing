
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
%  Topology optimization with simulated annealing
%  It is the monotop2 version with all function embedded
%  and some parameter to execute and save the files 
%
%
function [Ji]=top3f(nelx, nely, penal, near, N, vold)
  initime=cputime;
  Tmax=100;% can be less or more. In the last run had not logical pattern till 50 degree
  Tmin=0.001;
  nov=nelx*nely; %number of elements, it will be the number of variable
  maxloop=N;
  C_i=ones(nely, nelx);
  alphaa=0.85; 
  te=Tmax;
  %generetae a current solution
  xi(1:nely, 1:nelx)=0.2;
  aux=FE(nelx, nely, xi, penal);
  Ji=computeCompliance(nelx, nely, xi, penal, aux);
  fun=[];
  temp=[];
  acepted=[];
  rejected=[];
  mean_J=[]; %store mean value of objective function for each temperature 
  max_J=[];
  min_J=[];
  mean_C=[]; %store mean value of cristalization factor
  mn=near; % qtd de elementos ao lado ser mudado
  cont=1;
  Fil=0;
  while te>Tmin
    a1=0; %acepted solution count for each temperature
    a2=0; %rejected solution for each temperature
    aux_J=[];
    aux_C=[];
    for ii=1:N
      [xj, alex, aley]=newpeturb(xi, C_i, nelx, nely,mn,vold);
      if a2 > 0.8*N;
          if Fil < 1;
          display('Filter started'); 
          Fil=Fil+1;
          end
      xj=filter3(xj, nelx, nely,2); %higher value for refined mesh
      end 
      aux2=FE(nelx, nely, xj, penal);
      Jn=computeCompliance(nelx, nely, xj, penal,aux2);
      dE=Jn-Ji;
      p=exp(-dE/te);
%       aux_C(i,1)=C_i(round(nely/2) , round(nelx/2)); %to plot the crystalization factor
%       aux_C(i, 2)=C_i( round(nely/2) , round(nelx*0.75) );
%       aux_C(i, 3) = C_i( round(nely/2), round(nelx*0.20));
% %       aux_C(i, 4) = C_i( round(nely*0.75), round(nelx/2));
      if dE<=0
        xi=xj; %the current solution receiver the newsolution
        Ji=Jn;
        a1=a1+1;
        aux_J(ii)=Ji;
        for i=aley-mn:aley+mn
          for j=alex-mn:alex+mn
            C_i(i, j)=C_i(i, j)-1;
            if C_i(i,j)<=0
              C_i(i, j)=1;
            end
          end
        end
      elseif rand()<=p
        xi=xj; %the current solution receiver the newsolution
        Ji=Jn;
        a1=a1+1;
        aux_J(i)=Ji;
        for i=aley-mn:aley+mn
          for j=alex-mn:alex+mn
            C_i(i, j)=C_i(i, j)-1;
            if C_i(i,j)<=0
              C_i(i, j)=1;
            end
          end
        end
      else
        for i=aley-mn:aley+mn
          for j=alex-mn:alex+mn
            if C_i(i, j)<20
              C_i(i, j)=C_i(i, j)+1;
            end
          end
        end
        aux_J(i)=Ji;
        a2=a2+1;
      end    
      if a1>maxloop
        display('break,')
        break
      end
    end
    acepted(cont)=a1;
    rejected(cont)=a2;
    max_J(cont)=max(aux_J);
    min_J(cont)=min(aux_J);
    mean_J(cont)=mean(aux_J);
%     mean_C1(cont)=mean(aux_C(:,1));
%     mean_C2(cont)=mean(aux_C(:,2));
%     mean_C3(cont)=mean(aux_C(:,3));
%     mean_C4(cont)=mean(aux_C(:,4));
    fun(cont)=Ji;
    temp(cont)=te;
    te=te*alphaa
    cont=cont+1;
    %figure(1);
    colormap(gray); imagesc(-xi); axis equal; axis tight; axis off;pause(1e-5);
    end
  %###### display the best value and plot Function vs Temperature
  
  volfracf=sum(sum(xi))/(nelx*nely);
  fprintf('the best value %f \n', Ji);
  fprintf('The volume frac of the best solution %f \n', volfracf)
  
  endtime=cputime;
  runtime=endtime-initime;
  fprintf('the Running time %f \n', runtime);
  %save top.mat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function to analysis thedeslocation of nods
%
function [U]=FE(nelx,nely,x,penal)
  [KE] = lk; 
  K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
  F = sparse(2*(nely+1)*(nelx+1),1);
  U = zeros(2*(nely+1)*(nelx+1),1);
  for elx = 1:nelx
    for ely = 1:nely
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
      K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
    end
  end
  % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
  pp=2*(nelx+1)*(nely+1);
%   F(2,1) = -1;
%   fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
  % Short Cantilever problem (commetn two previous lines)
  F(2*(nelx+1)*(nely+1),1)=-1;
  fixeddofs = [1:2*(nely+1)];
  alldofs     = [1:2*(nely+1)*(nelx+1)];
  freedofs    = setdiff(alldofs,fixeddofs);
  % SOLVING
  U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
  U(fixeddofs,:)= 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% function to create the stiffines matrix 
function [KE]=lk
  E = 1.; 
  nu = 0.3;
  k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
     -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
  KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
% Function to compute the compliance
%
function [J]=computeCompliance(nelx,nely,x,penal, U) %elemsIn,voidEps,
  %Compute the compliance of the entire mesh
  [KE] = lk;J = 0;
  for elx = 1:nelx
      for ely = 1:nely
          n1 = (nely+1)*(elx-1)+ely; 
          n2 = (nely+1)* elx   +ely;
          edof = [2*n1-1 2*n1 2*n2-1 2*n2 2*n2+1 2*n2+2 2*n1+1 2*n1+2]';
          %alpha = (1-elemsIn(ely,elx))*voidEPs + elemsIn(ely,elx);
          Ue = U(edof);
          J = J + x(ely, elx)^penal*Ue'*KE*Ue; % transposta deslocamento* matriz rigidez * deslocamento
      end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function to generate a new solution 
%
function [xj, alex, aley]=newpeturb(xc, C_i, nelx, nely,nn,vold);
  inte=0;
  crit=0;
  xj=xc;
  while crit ~=1 || inte<=20
    delr=1;
    e_i=1;
    delx=nelx-nn;
    dely=nely-nn;
    alex=randi([nn+1, delx]);
    aley=randi([nn+1, dely]);
    aux=0;
    inte=inte+1;
    %display(inte)
    for ij=1:C_i(aley, alex)
      aux=aux+(rand()-0.5);
    end
    
    for ik=(aley-nn):(aley+nn)
      for jk=(alex-nn):(alex+nn)
        delta=delr*e_i*aux/C_i(ik, jk);
        xj(ik,jk)=xc(ik,jk)+delta;
        if xj(ik, jk)>1
          xj(ik, jk)=1;
        elseif xj(ik, jk)<=0
          xj(ik, jk)=0.001;
        end
      end
    end
    crit=res(nelx, nely, xj,vold);
  end
end %end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function with the resctriction  
%
%
function [aux]=res(nelx, nely, x, vold)
  volu=sum(sum(x))/(nelx*nely);
  if volu<vold 
      aux=1;
  else
      aux=0;
  end
end



function [x_fil]=filter3(xj, nelx, nely, radius)
  x_fil=zeros(nely, nelx);
  %Loop for each i element
  for row = 1:nely 
    for columm = 1:nelx
        sums=0;
        w_i_j=0;
        w_i_js=0;
        % loop to determine j elements
        for y_el = max(row - round(radius), 1):min(row + round(radius), nely)
          for x_el = max(columm -round(radius), 1 ):min(columm + round(radius), nelx)
            if sqrt((x_el-columm)^2+(y_el-row)^2) < radius
              d_i_j = radius - sqrt((y_el - row)^2 + (x_el - columm)^2); %distance between i and j
              w_i_js= w_i_js + d_i_j;
              sums = sums + d_i_j*xj(y_el, x_el);
            end
          end
        end
        aux=sums/w_i_js;
%         if aux < 0.5
%             aux=0.001;
%         else aux=1;
%         end
x_fil(row, columm) = aux;
        
    end
  end
  %function end
end