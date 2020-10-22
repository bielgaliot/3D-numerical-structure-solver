function l = computeL(L1, L2, M, Me, g)

q1 = @(x) 0.85-0.15.*cos(pi.*x./L1);

q2 = @(y) -(1./L2.^2).*(L1-L2-y).*(L1+L2-y);

lambda1 = @(z) (3.*M./(2.*L1.^2)).*(L1-z)+M./(4.*(L1+L2));

lambda2 =  (L2+L1-L1).*(M./(4.*(L1+L2)));



Lift = integral(q1,0,L1) + integral(q2,L1,L2+L1);

Mass = integral(lambda1,0,L1) + lambda2 + Me;

l=g*Mass/Lift;

end