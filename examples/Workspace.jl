using ApproxFun, SingularIntegralEquations, DualNumbers, Holomorphic

Disk()

points(Circle()

x=Fun()
w=sqrt(1-x^2)
g=stieltjes(w)
g/(-2π*im)
g|(2Circle())
imag(g|+)
g|-
∂(g)

domain(g)

-∂
cauchy(+1

f(.1+.2im)




Γ=ComplexPlane()\Interval()
S=Space(Γ)




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
fromcanonical(d,dual(0.,-1.im))

Dual(a::Real,b::Complex)=Dual(a+0*b,b)

sqrt(Dual(-1.,1.im))
sqrt(Dual(

sqrt(-1.-0.im)
sqrt(dual(-1.,-1im))
sqrt(dual(-1.+0.im,+1im))
sqrt(dual(-1.-0.im,-1im))

realpart(sqrt(dual(-1.,-im)))
realpart(sqrt(dual(-1.,im)))

f(z)=sqrt((1.+1.im)*z)

realpart(f(dual(-1/(1.+1.im),-1.im)))

f(-1/(1.+1.im)-0.000001im)



sqrt(-1.)
Pkg.checkout("DualNumbers")

