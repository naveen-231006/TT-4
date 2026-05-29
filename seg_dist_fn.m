function [dist, dir] = seg_dist_fn(p1,p2,p3,p4)
% Minimum distance between two line segments p1-p2 and p3-p4
u=p2-p1; v=p4-p3; w=p1-p3;
a=dot(u,u); b=dot(u,v); c=dot(v,v); d=dot(u,w); e=dot(v,w);
D=a*c-b*b;
if D<1e-8, sc=0; if b>c, tc=d/b; else, tc=e/c; end
else, sc=(b*e-c*d)/D; tc=(a*e-b*d)/D; end
sc=max(0,min(1,sc)); tc=max(0,min(1,tc));
pt1=p1+sc*u; pt2=p3+tc*v; dv=pt1-pt2;
dist=norm(dv);
if dist<1e-9, dir=[0 0 1]; else, dir=dv/dist; end
end
