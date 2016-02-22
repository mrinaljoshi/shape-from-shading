function [p,q,z] = cal_pqz(E,p_bndy,q_bndy,z_bndy,ps,qs,lambda,weight,bc1,bc2)
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
    
    % restoring values at bc 
    p_up(bc1) = p_bndy(bc1);
    p_up(bc2) = p_bndy(bc2);
    q_up(bc1) = q_bndy(bc1);
    q_up(bc2) = q_bndy(bc2);
    
    R = (p_up*ps+q_up*qs+1)./sqrt((p_up.^2+q_up.^2+1)*(ps^2+qs^2+1));
    p = p_up;
    q = q_up;
    
end

% figure()
% mesh(x,y,abs(R))

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