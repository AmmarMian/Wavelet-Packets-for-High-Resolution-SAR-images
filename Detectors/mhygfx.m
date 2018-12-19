function hf = mhygfx(x,a,b,c)
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).

%     ============================================================
%     Purpose: This program computes the hypergeometric function
%     F(a,b,c,x) using subroutine HYGFX
%     Input :  a --- Parameter
%     b --- Parameter
%     c --- Parameter, c <> 0,-1,-2,...
%     x --- Argument ( x ף 1 )
%     Output:  HF --- F(a,b,c,x)
%     Example:
%     b = 3.30,  c = 6.70
%     a     F(a,b,c,.25)     F(a,b,c,.55)    F(a,b,c,.85)
%     ------------------------------------------------------
%     -2.5   .72356129D+00    .46961432D+00   .29106096D+00
%     -0.5   .93610145D+00    .85187390D+00   .75543187D+00
%     0.5   .10689695D+01    .11795358D+01   .13510497D+01
%     2.5   .14051563D+01    .23999063D+01   .57381566D+01

%     a = 3.30,  b = 6.70
%     c     F(a,b,c,.25)     F(a,b,c,.55)    F(a,b,c,.85)
%     ------------------------------------------------------
%     -5.5   .15090670D+05    .10170778D+11   .58682088D+19
%     -0.5  -.21631479D+04   -.30854772D+07  -.10217370D+13
%     0.5   .26451677D+03    .11967860D+06   .92370648D+10
%     4.5   .41946916D+01    .58092729D+02   .20396914D+05
%     ============================================================

hf = zeros(size(x));

for k=1:length(x),
    
  hf(k)=hygfx(a,b,c,x(k),[]);
  
end


function hf=hygfx(a,b,c,x,hf);

%     ====================================================
%     Purpose: Compute hypergeometric function F(a,b,c,x)
%     Input :  a --- Parameter
%     b --- Parameter
%     c --- Parameter, c <> 0,-1,-2,...
%     x --- Argument   ( x < 1 )
%     Output:  HF --- F(a,b,c,x)
%     Routines called:
%     (1) GAMMA for computing gamma function
%     (2) PSI for computing psi function
%     ====================================================



gc=[];gcab=[];gca=[];gcb=[];g1=[];g2=[];g3=[];ga=[];gb=[];gam=[];gbm=[];pa=[];pb=[];gabc=[];
hw=0.0;
pi=3.141592653589793d0;
el=.5772156649015329d0;
l0=c == fix(c)&c < 0.0;
l1=1.0d0-x < 1.0d-15&c-a-b <= 0.0;
l2=a == fix(a)&a < 0.0;
l3=b == fix(b)&b < 0.0;
l4=c-a == fix(c-a)&c-a <= 0.0;
l5=c-b == fix(c-b)&c-b <= 0.0;
if (l0|l1);
'the hypergeometric series is divergent',
return;
end;
eps=1.0d-15;
if (x > 0.95) eps=1.0d-8; end;
if (x == 0.0|a == 0.0|b == 0.0);
hf=1.0d0;
return;
elseif (1.0d0-x == eps&c-a-b > 0.0);
[gc]=gamma(c,gc);
[gcab]=gamma(c-a-b,gcab);
[gca]=gamma(c-a,gca);
[gcb]=gamma(c-b,gcb);
hf=gc.*gcab./(gca.*gcb);
return;
elseif (1.0d0+x <= eps&abs(c-a+b-1.0) <= eps);
g0=sqrt(pi).*2.0d0.^(-a);
[g1]=gamma(c,g1);
[g2]=gamma(1.0d0+a./2.0-b,g2);
[g3]=gamma(0.5d0+0.5.*a,g3);
hf=g0.*g1./(g2.*g3);
return;
elseif (l2|l3);
if (l2) nm=fix(abs(a)); end;
if (l3) nm=fix(abs(b)); end;
hf=1.0d0;
r=1.0d0;
for  k=1:nm;
r=r.*(a+k-1.0d0).*(b+k-1.0d0)./(k.*(c+k-1.0d0)).*x;
hf=hf+r;
end;
return;
elseif (l4|l5);
if (l4) nm=fix(abs(c-a)); end;
if (l5) nm=fix(abs(c-b)); end;
hf=1.0d0;
r=1.0d0;
for  k=1:nm;
r=r.*(c-a+k-1.0d0).*(c-b+k-1.0d0)./(k.*(c+k-1.0d0)).*x;
hf=hf+r;
end;
hf=(1.0d0-x).^(c-a-b).*hf;
return;
end;
aa=a;
bb=b;
x1=x;
if (x < 0.0d0);
x=x./(x-1.0d0);
if (c > a&b < a&b > 0.0);
a=bb;
b=aa;
end;
b=c-b;
end;
if (x >= 0.75d0);
gm=0.0d0;
if (abs(c-a-b-fix(c-a-b)) < 1.0d-15);
m=fix(c-a-b);
[ga]=gamma(a,ga);
[gb]=gamma(b,gb);
[gc]=gamma(c,gc);
[gam]=gamma(a+m,gam);
[gbm]=gamma(b+m,gbm);
[a,pa]=psi(a,pa);
[b,pb]=psi(b,pb);
if (m ~= 0) gm=1.0d0; end;
for  j=1:abs(m)-1;
gm=gm.*j;
end;
rm=1.0d0;
for  j=1:abs(m);
rm=rm.*j;
end;
f0=1.0d0;
r0=1.0d0;
r1=1.0d0;
sp0=0.d0;
sp=0.0d0;
if (m >= 0);
c0=gm.*gc./(gam.*gbm);
c1=-gc.*(x-1.0d0).^m./(ga.*gb.*rm);
for  k=1:m-1;
r0=r0.*(a+k-1.0d0).*(b+k-1.0)./(k.*(k-m)).*(1.0-x);
f0=f0+r0;
end;
for  k=1:m;
sp0=sp0+1.0d0./(a+k-1.0)+1.0./(b+k-1.0)-1.0./k;
end;
f1=pa+pb+sp0+2.0d0.*el+log(1.0d0-x);
for  k=1:250;
sp=sp+(1.0d0-a)./(k.*(a+k-1.0))+(1.0-b)./(k.*(b+k-1.0));
sm=0.0d0;
for  j=1:m;
sm=sm+(1.0d0-a)./((j+k).*(a+j+k-1.0))+1.0./ (b+j+k-1.0);
end;
rp=pa+pb+2.0d0.*el+sp+sm+log(1.0d0-x);
r1=r1.*(a+m+k-1.0d0).*(b+m+k-1.0)./(k.*(m+k)).*(1.0-x);
f1=f1+r1.*rp;
if (abs(f1-hw) < abs(f1).*eps) break; end;
hw=f1;
end;
      hf=f0.*c0+f1.*c1;
elseif (m < 0);
m=-m;
c0=gm.*gc./(ga.*gb.*(1.0d0-x).^m);
c1=-(-1).^m.*gc./(gam.*gbm.*rm);
for  k=1:m-1;
r0=r0.*(a-m+k-1.0d0).*(b-m+k-1.0)./(k.*(k-m)).*(1.0-x);
f0=f0+r0;
end;
for  k=1:m;
sp0=sp0+1.0d0./k;
end;
f1=pa+pb-sp0+2.0d0.*el+log(1.0d0-x);
for  k=1:250;
sp=sp+(1.0d0-a)./(k.*(a+k-1.0))+(1.0-b)./(k.*(b+k-1.0));
sm=0.0d0;
for  j=1:m;
sm=sm+1.0d0./(j+k);
end;
rp=pa+pb+2.0d0.*el+sp-sm+log(1.0d0-x);
r1=r1.*(a+k-1.0d0).*(b+k-1.0)./(k.*(m+k)).*(1.0-x);
f1=f1+r1.*rp;
if (abs(f1-hw) < abs(f1).*eps) break; end;
hw=f1;
end;
      hf=f0.*c0+f1.*c1;
end;
else;
[ga]=gamma(a,ga);
[gb]=gamma(b,gb);
[gc]=gamma(c,gc);
[gca]=gamma(c-a,gca);
[gcb]=gamma(c-b,gcb);
[gcab]=gamma(c-a-b,gcab);
[gabc]=gamma(a+b-c,gabc);
c0=gc.*gcab./(gca.*gcb);
c1=gc.*gabc./(ga.*gb).*(1.0d0-x).^(c-a-b);
hf=0.0d0;
r0=c0;
r1=c1;
for  k=1:250;
r0=r0.*(a+k-1.0d0).*(b+k-1.0)./(k.*(a+b-c+k)).*(1.0-x);
r1=r1.*(c-a+k-1.0d0).*(c-b+k-1.0)./(k.*(c-a-b+k)) .*(1.0-x);
hf=hf+r0+r1;
if (abs(hf-hw) < abs(hf).*eps) break; end;
hw=hf;
end;
     hf=hf+c0+c1;
end;
else;
a0=1.0d0;
if (c > a&c < 2.0d0.*a& c > b&c < 2.0d0.*b);
a0=(1.0d0-x).^(c-a-b);
a=c-a;
b=c-b;
end;
hf=1.0d0;
r=1.0d0;
for  k=1:250;
r=r.*(a+k-1.0d0).*(b+k-1.0d0)./(k.*(c+k-1.0d0)).*x;
hf=hf+r;
if (abs(hf-hw) <= abs(hf).*eps) break; end;
hw=hf;
end;
   hf=a0.*hf;
end;
if (x1 < 0.0d0);
x=x1;
c0=1.0d0./(1.0d0-x).^aa;
hf=c0.*hf;
end;
a=aa;
b=bb;
if (k > 120) ; end;
  %format(1x,'warning% you should check the accuracy');



function [ga]=gamma(x,ga);

%     ==================================================
%     Purpose: Compute gamma function ג(x)
%     Input :  x  --- Argument of ג(x)
%     ( x is not equal to 0,-1,-2,תתת)
%     Output:  GA --- ג(x)
%     ==================================================



pi=3.141592653589793d0;
if (x == fix(x));
if (x > 0.0d0);
ga=1.0d0;
m1=x-1;
for  k=2:m1;
ga=ga.*k;
end;
else;
ga=1.0d+300;
end;
else;
if (abs(x) > 1.0d0);
z=abs(x);
m=fix(z);
r=1.0d0;
for  k=1:m;
r=r.*(z-k);
end;
z=z-m;
else;
z=x;
end;
g=[1.0d0,0.5772156649015329d0,-0.6558780715202538d0,-0.420026350340952d-1,0.1665386113822915d0,-.421977345555443d-1,-.96219715278770d-2,.72189432466630d-2,-.11651675918591d-2,-.2152416741149d-3,.1280502823882d-3,-.201348547807d-4,-.12504934821d-5,.11330272320d-5,-.2056338417d-6,.61160950d-8,.50020075d-8,-.11812746d-8,.1043427d-9,.77823d-11,-.36968d-11,.51d-12,-.206d-13,-.54d-14,.14d-14,.1d-15];
gr=g(26);
for  k=25:-1:1;
gr=gr.*z+g(k);
end;
ga=1.0d0./(gr.*z);
if (abs(x) > 1.0d0);
ga=ga.*r;
if (x < 0.0d0) ga=-pi./(x.*ga.*sin(pi.*x)); end;
end;
end;



function [x,ps]=psi(x,ps);

%     ======================================
%     Purpose: Compute Psi function
%     Input :  x  --- Argument of psi(x)
%     Output:  PS --- psi(x)
%     ======================================



xa=abs(x);
pi=3.141592653589793d0;
el=.5772156649015329d0;
s=0.0d0;
if (x == fix(x)&x <= 0.0);
ps=1.0d+300;
return;
elseif (xa == fix(xa));
n=xa;
for  k=1 :n-1;
s=s+1.0d0./k;
end;
ps=-el+s;
elseif (xa+.5 == fix(xa+.5));
n=xa-.5;
for  k=1:n;
s=s+1.0./(2.0d0.*k-1.0d0);
end;
ps=-el+2.0d0.*s-1.386294361119891d0;
else;
if (xa < 10.0);
n=10-fix(xa);
for  k=0:n-1;
s=s+1.0d0./(xa+k);
end;
xa=xa+n;
end;
x2=1.0d0./(xa.*xa);
a1=-.8333333333333d-01;
a2=.83333333333333333d-02;
a3=-.39682539682539683d-02;
a4=.41666666666666667d-02;
a5=-.75757575757575758d-02;
a6=.21092796092796093d-01;
a7=-.83333333333333333d-01;
a8=.4432598039215686d0;
ps=log(xa)-.5d0./xa+x2.*(((((((a8.*x2+a7).*x2+a6).*x2+a5).*x2+a4).*x2+a3).*x2+a2).*x2+a1);
ps=ps-s;
end;
if (x < 0.0) ps=ps-pi.*cos(pi.*x)./sin(pi.*x)-1.0d0./x; end;

