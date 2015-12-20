using ApproxFun, SingularIntegralEquations, DualNumbers, Holomorphic


f=stieltjes(w)

f(.1+.2im)




Γ=ComplexPlane()\Interval()
S=Space(Γ)


x=Fun()
w=sqrt(1-x^2)


f(z)=stieltjes(w,z)

g=Fun(Fun(x->f(x+0im)-f(x-0im),S.space).coefficients/(-2π*im),S)

g(.1+.2im)
f(.1+.2im)

stieltjes(w).coefficients

g.coefficients

stieltjes(1/w,1.0001)




d=Interval(-1.-1im,1+1.im)

fromcanonical(d,0.1)-fromcanonical(d,0.1+0.im)
fromcanonical(d,0.1)

fromcanonical(d,0.-0.im)
fromcanonical(d,0.+0.im)

using DualNumbers

fromcanonical(d,dual(0.,1.im))

Pkg.checkout("DualNumbers")

