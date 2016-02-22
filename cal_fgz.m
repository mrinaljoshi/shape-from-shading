function [f,g,z] = cal_fgz(E,f_bndy,g_bndy,z_bndy,ps,qs,lambda,weight,bc1,bc2)
kernel = [0 1 0; 1 0 1; 0 1 0]/4; % filter for averaging

f = f_bndy;
g = g_bndy;
R = zeros(size(f));

for i = 1:1000

    Rf = 4*ps./(sqrt(1+ps^2+qs^2)*(4+f.^2+g.^2))-16*f./(sqrt(1+ps^2+qs^2)*(4+f.^2+g.^2).^2);
    Rg = 4*qs./(sqrt(1+ps^2+qs^2)*(4+f.^2+g.^2))-16*g./(sqrt(1+ps^2+qs^2)*(4+f.^2+g.^2).^2);
    f_avg = conv2(f,kernel,'same');
    g_avg = conv2(g,kernel,'same');
    f_up = f_avg - 1/lambda*(E-R).*Rf;
    g_up = g_avg - 1/lambda*(E-R).*Rg;
    
    % restoring values at bc 
    f_up(bc1) = f_bndy(bc1);
    f_up(bc2) = f_bndy(bc2);
    g_up(bc1) = g_bndy(bc1);
    g_up(bc2) = g_bndy(bc2);
    
    R = (4*ps*f_up+4*qs*g_up+4-f_up.^2-g_up.^2)./((4+f_up.^2+g_up.^2).*sqrt(1+ps^2+qs^2));
    f = f_up;
    g = g_up;
    
end

% figure()
% mesh(x,y,abs(R))

p_up = 4*f./(4-f.^2-g.^2);
q_up = 4*g./(4-f.^2-g.^2);

[px,py] = gradient(p_up);
[qx,qy] = gradient(q_up);


z   = z_bndy;

kernel = [0 1 0; 1 0 1; 0 1 0]/4;
for i = 1:1000;
    z_avg = conv2(z,kernel,'same');
    z_uf = z_avg + weight*(px+qy);
    z = z_uf;
    z(bc1) = z_bndy(bc1);
    z(bc2) = z_bndy(bc2);
end
