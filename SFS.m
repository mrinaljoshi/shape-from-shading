%% Part 1:  Storing Values of p q and z at boundaries and  then calculating E 
close all
clear all
[x,y] = meshgrid(-10:0.5:10);
z = 30 -x.^2/4- y.^2/9;
figure, mesh(x,y,z)
title('Surface to be reconstructed');
% source is at (0,0)
ps = 0;
qs = 0;

% calculating E and Boundary conditions  
[p,q] = gradient(z);
E = (p*ps+q*qs+1)./sqrt((p.^2+q.^2+1)*(ps^2+qs^2+1));

f = 2*(sqrt(p.^2+q.^2+1)-1)./(p.^2+q.^2+eps).*p;
g = 2*(sqrt(p.^2+q.^2+1)-1)./(p.^2+q.^2+eps).*q;

p_bndy = zeros(size(p));
q_bndy = zeros(size(q));
f_bndy = zeros(size(f));
g_bndy = zeros(size(q));
z_bndy = zeros(size(z));

% Take a closed curve by changing boundary condition on x and y by b_x and b_y
b_x = 36; 
b_y = 25;
bc1 = x.^2==b_x & y.^2<=b_y;
bc2 = y.^2==b_y & x.^2<=b_x;

p_bndy(bc1) = p(bc1);
p_bndy(bc2) = p(bc2);
q_bndy(bc1) = q(bc1);
q_bndy(bc2) = q(bc2);
z_bndy(bc1) = z(bc1);
z_bndy(bc2) = z(bc2);
f_bndy(bc1) = f(bc1);
f_bndy(bc2) = f(bc2);
g_bndy(bc1) = g(bc1);
g_bndy(bc2) = g(bc2);

f_origin = f;
g_origin = g;
p_origin = p;
q_origin = q;
z_origin = z;
%% Part 2A: Estimating of p, q , z and comparing with original

weight = -0.1;
lambda = 10;
[p0,q0,z] = cal_pqz(E,p_bndy,q_bndy,z_bndy,ps,qs,lambda,weight,bc1,bc2);
figure()
mesh(x,y,z)
title('Estimated depth map');
error = sqrt(mean((z_origin(:)-z(:)).^2));
error_p = sqrt(mean((p_origin(:)-p0(:)).^2));
error_q = sqrt(mean((q_origin(:)-q0(:)).^2));
fprintf('The mean square error in estimation of Z is %f for lambda %d \n ',error,lambda );
fprintf('The mean square error in estimation of p is %f for lambda %d \n ',error_p,lambda );
fprintf('The mean square error in estimation of q is %f for lambda %d \n ',error_q,lambda );

%% Part 2B: Repeating for different lambda and weight  

weight = -0.1;
error = zeros(1,7); i=0;
for lambda = [0.01,0.1,1,10,100,1000,10000];
[p1,q1,z] = cal_pqz(E,p_bndy,q_bndy,z_bndy,ps,qs,lambda,weight,bc1,bc2);
i =i+1;
error(i) = sqrt(mean((z_origin(:)-z(:)).^2));
fprintf('The mean square error in estimation of Z is %f for lambda %d\n ',error,lambda);
end
figure()
plot(error)
title('error plot versus increasing lambda values')

% For lambda > 1, the error becomes almost a constant. It decreases
% monotonically as lambda is increased.

%% Part 2C: Adding AWGN to E for reconstruction

E_new = E + wgn(size(E,1),size(E,2),0.1);
figure()
mesh(x,y,E_new)
title('E with added wgn of variance 0.1')
weight = -0.1;
lambda = 10;
[p2,q2,z] = cal_pqz(E_new,p_bndy,q_bndy,z_bndy,ps,qs,lambda,weight,bc1,bc2);
error = sqrt(mean((z_origin(:)-z(:)).^2));
fprintf('The mean square error in estimation of Z is %f for variance 0.1 \n ',error);

%% Part 2D: Adding AWGN of different variances to E for reconstruction

i = 0;
error1 = zeros(1,9);
for var = [0.001,0.01,0.1,1,10,20,30,40,50]
E_new = E + wgn(size(E,1),size(E,2),var);
weight = -0.1;
lambda = 10;
i=i+1;
[p2,q2,z] = cal_pqz(E_new,p_bndy,q_bndy,z_bndy,ps,qs,lambda,weight,bc1,bc2);
error1(i) = sqrt(mean((z_origin(:)-z(:)).^2));
fprintf('The mean square error in estimation of Z is %f for variance %f \n ',error1(i),var);
end
figure()
plot(error1)
title('Error plot vs various noise variances')

% Error increases on increasing the variance of AWGN. This can be seen in
% the plot.
%% Part 2E: Estimate depth using f and g 

weight = -0.1;
lambda = 10;
[f0,g0,z] = cal_fgz(E,f_bndy,g_bndy,z_bndy,ps,qs,lambda,weight,bc1,bc2);
figure()
mesh(x,y,z)
title('Estimated z from f,g');
error = sqrt(mean((z_origin(:)-z(:)).^2));
fprintf('The mean square error in estimation of Z is %f for lambda %d \n ',error,lambda );

%% Part 2F: Perturb S and reconstruct

perturb = 0.51;
ps_prtb = ps + perturb;
qs_prtb = qs + perturb;

weight = -0.1;
lambda = 10;
[p3,q3,z] = cal_pqz(E_new,p_bndy,q_bndy,z_bndy,ps_prtb,qs_prtb,lambda,weight,bc1,bc2);
figure()
mesh(x,y,z)
title('Estimated z with a perturbed source vector');
error = sqrt(mean((z_origin(:)-z(:)).^2));
fprintf('The mean square error in estimation of Z is %f for ps %f and qs %f \n ',error,ps_prtb,qs_prtb );

%% Part 2G: Reconstructing z with a different viewing direction 

ps_new = 10;
qs_new = 9.75;

weight = -0.1;
lambda = 10;
[p4,q4,z] = cal_pqz(E_new,p_bndy,q_bndy,z_bndy,ps_new,qs_new,lambda,weight,bc1,bc2);
figure()
mesh(x,y,z)
title('Estimated z with a different source vector');
error = sqrt(mean((z_origin(:)-z(:)).^2));
fprintf('The mean square error in estimation of Z is %f for ps %f and qs %f \n ',error,ps_new,qs_new );

%% Part 2H: Recovering z from two sources with relative strength alpha and beta 

ps1 = 0.25; qs1 = 0.4;
ps2 = 0.5; qs2 = 0.8;
alpha = 0.6; beta = 0.4;

E_new = (alpha*(p*ps1+q*qs1+1)./sqrt((p.^2+q.^2+1)*(ps1^2+qs1^2+1)) )+ (beta*(p*ps2+q*qs2+1)./sqrt((p.^2+q.^2+1)*(ps2^2+qs2^2+1)));
kernel = [0 1 0; 1 0 1; 0 1 0]/4; % filter for averaging
p = p_bndy;
q = q_bndy;
R = zeros(size(p));

for i = 1:1000
    mod_pq = sqrt(p.^2+q.^2+1);
    mod_psqs1 = sqrt(ps1^2+qs1^2+1);
    mod_psqs2 = sqrt(ps2^2+qs2^2+1);
    Rp1 = ps1./(mod_pq*mod_psqs1) - (p*ps1+q*qs1+1).*p./(mod_pq.^3 *mod_psqs1.^2);
    Rq1 = qs1./(mod_pq*mod_psqs1) - (p*ps1+q*qs1+1).*q./(mod_pq.^3* mod_psqs1.^2);
    Rp2 = ps2./(mod_pq*mod_psqs2) - (p*ps2+q*qs2+1).*p./(mod_pq.^3 *mod_psqs2.^2);
    Rq2 = qs2./(mod_pq*mod_psqs2) - (p*ps2+q*qs2+1).*q./(mod_pq.^3* mod_psqs2.^2);
    Rp = alpha*Rp1 + beta*Rp2 ;
    Rq = alpha*Rq1 + beta*Rq2 ;
    
    p_avg = conv2(p,kernel,'same');
    q_avg = conv2(q,kernel,'same');
    p_up = p_avg - 1/lambda*(E_new-R).*Rp1;
    q_up = q_avg - 1/lambda*(E_new-R).*Rq1;
    
    % restoring values at bc 
    p_up(bc1) = p_bndy(bc1);
    p_up(bc2) = p_bndy(bc2);
    q_up(bc1) = q_bndy(bc1);
    q_up(bc2) = q_bndy(bc2);
    
    R1 = (p_up*ps1+q_up*qs1+1)./sqrt((p_up.^2+q_up.^2+1)*(ps1^2+qs1^2+1));
    R2 = (p_up*ps2+q_up*qs2+1)./sqrt((p_up.^2+q_up.^2+1)*(ps2^2+qs2^2+1));
    R = alpha*R1 + beta*R2 ;
    p = p_up;
    q = q_up;   
end

figure()
mesh(x,y,abs(R))

[px,py] = gradient(p_up);
[qx,qy] = gradient(q_up);
z   = z_bndy;
kernel = [0 1 0; 1 0 1; 0 1 0]/4;
for i = 1:1000;
    z_avg = conv2(z,kernel,'same');
    z_up = z_avg + weight*(px+qy);
    z = z_up;
    z(bc1) = z_bndy(bc1);
    z(bc2) = z_bndy(bc2);
end
figure()
mesh(x,y,z)
title('Estimated z with two source vectors');
error = sqrt(mean((z_origin(:)-z(:)).^2));
fprintf('The mean square error in estimation of Z for source1 at (%f,%f) and source2 at (%f,%f) is %f \n ',ps1,qs1,ps2,qs2,error );

%% Part 2I: Experiment with real images

load moz256.depth
z= reshape(moz256,256,256);
[x,y] = meshgrid(-127:1:128);
figure()
mesh(x,y,z);
title('Original depth map of mozart bust image');
z_origin = z;
ps = 0;
qs = 0;

[p,q] = gradient(z);
E = (p*ps+q*qs+1)./sqrt((p.^2+q.^2+1)*(ps^2+qs^2+1));
figure()
mesh(x,y,E)
title('E')
f = 2*(sqrt(p.^2+q.^2+1)-1)./(p.^2+q.^2+eps).*p;
g = 2*(sqrt(p.^2+q.^2+1)-1)./(p.^2+q.^2+eps).*q;

% for getting boundary condition detect edges and store those values in p_bndy and q_bndy.
B = bwboundaries(z);
edge1 = edge(z);

p_bndy = p.*edge1;
q_bndy = q.*edge1;
f_bndy = f.*edge1;
g_bndy = g.*edge1;
z_bndy = z.*edge1;

weight = -0.1;
lambda = 10;
kernel = [0 1 0; 1 0 1; 0 1 0]/4; % filter for averaging
p = p_bndy;
q = q_bndy;
R = zeros(size(p));

for i = 1:1000
    mod_pq = sqrt(p.^2+q.^2+1);
    mod_psqs = sqrt(ps^2+qs^2+1);
    Rp = ps./(mod_pq*mod_psqs) - (p*ps+q*qs+1).*p./(mod_pq.^3 *mod_psqs.^2);
    Rq = qs./(mod_pq*mod_psqs) - (p*ps+q*qs+1).*q./(mod_pq.^3* mod_psqs.^2);
    p_avg = conv2(p,kernel,'same');
    q_avg = conv2(q,kernel,'same');
    p_up = p_avg - 1/lambda*(E-R).*Rp;
    q_up = q_avg - 1/lambda*(E-R).*Rq;
    
    % restoring values at edge 
    p_up(find(edge1)) = p_bndy(find(edge1));
    q_up(find(edge1)) = q_bndy(find(edge1));
    
    R = (p_up*ps+q_up*qs+1)./sqrt((p_up.^2+q_up.^2+1)*(ps^2+qs^2+1));
    p = p_up;
    q = q_up;
    
end

[px,py] = gradient(p_up);
[qx,qy] = gradient(q_up);
z   = z_bndy;

kernel = [0 1 0; 1 0 1; 0 1 0]/4;
for i = 1:1000;
    z_avg = conv2(z,kernel,'same');
    z_up = z_avg + weight*(px+qy);
    z = z_up;
    z(find(edge1)) = z_bndy(find(edge1));
    
end

figure()
mesh(x,y,z)
title('Reconstructed depth map for mozart bust image');
error = sqrt(mean((z_origin(:)-z(:)).^2));
fprintf('The mean square error in estimation of Z is %f for lambda %d \n ',error,lambda );
