__precompile__()
module Holomorphic
    using Base, ApproxFun, SingularIntegralEquations


import ApproxFun: UnivariateSpace, domain, evaluate, ComplexBasis, spacescompatible,
                    Space, defaultFun, union_rule, ConstantSpace, points, checkpoints,
                    transform, plan_transform
import SingularIntegralEquations: stieltjes


export ℂ

# represents the whole Planess
immutable ComplexPlane <: Domain{Complex128,2} end

const ℂ=ComplexPlane()

Base.intersect(a::ComplexPlane,b::ComplexPlane)=a


# represents S\A
immutable Complement{SS,AA,T,d} <: Domain{T,d}
    full::SS
    del::AA
end

Base.intersect(a::Complement,b::Complement)=Complement(a.full ∩ b.full,a.del ∪ b.del)

Complement(A::Domain,B::Domain)=Complement{typeof(A),typeof(B),eltype(A),ndims(A)}(A,B)
Base.setdiff(a::Domain,b::Domain)=Complement(a,b)


immutable StieltjesSpace{S,DD} <:  Space{ComplexBasis,Complement{ComplexPlane,DD},2}
    space::S
    function StieltjesSpace(sp::S)
        @assert isa(domain(sp),DD)
        new(sp)
    end
end

StieltjesSpace(sp::UnivariateSpace)=StieltjesSpace{typeof(sp),typeof(domain(sp))}(sp)
spacescompatible(a::StieltjesSpace,b::StieltjesSpace)=spacescompatible(a.space,b.space)

domain(sp::StieltjesSpace)=ℂ\domain(sp.space)
evaluate(v::AbstractVector,sp::StieltjesSpace,z)=v[1]+stieltjes(sp.space,v[2:end],z)
stieltjes(f::Fun)=Fun([0;f.coefficients],StieltjesSpace(space(f)))


union_rule(A::ConstantSpace,B::StieltjesSpace)=B
Base.ones{T<:Number}(::Type{T},S::StieltjesSpace)=Fun(ones(T,1),S)
Base.ones(S::StieltjesSpace)=Fun(ones(1),S)


# construct spaces for complement
# we assume boundedness

Space{DD<:Interval}(Γ::Complement{ComplexPlane,DD})=StieltjesSpace(JacobiWeight(0.5,0.5,Ultraspherical{1}(Γ.del)))
Space{DD<:Circle}(Γ::Complement{ComplexPlane,DD})=StieltjesSpace(Laurent(Γ.del))
Space{DD<:UnionDomain}(Γ::Complement{ComplexPlane,DD})=StieltjesSpace(PiecewiseSpace(map(d->Space(d).space,map(x->ComplexPlane()\x,Γ.del))))



## Interval

checkpoints{DD<:Interval}(C::Complement{ComplexPlane,DD})=[fromcanonical(C.del,1im),fromcanonical(C.del,2.)]

joukowsky(z)=(z+1./z)./2

points{DD<:Interval}(S::StieltjesSpace{JacobiWeight{Ultraspherical{1,DD},DD},DD},n)=
    fromcanonical(domain(S).del,joukowsky(points((1-1/n)*Circle(),n)))

plan_transform{DD<:Interval}(S::StieltjesSpace{JacobiWeight{Ultraspherical{1,DD},DD},DD},vals::Vector)=
        plan_transform(Taylor((1-1/length(vals))*Circle()),vals)

function transform{DD<:Interval}(S::StieltjesSpace{JacobiWeight{Ultraspherical{1,DD},DD},DD},vals,plan)
    @assert S.space.α==S.space.β==0.5
    n=length(vals)
    cfs=transform(Taylor((1-1/n)*Circle()),vals,plan)
    for k=2:length(cfs)
        cfs[k]*=1/((1-1/n)^(k-1)*π)
    end
    cfs
end


## ℂ\Circle()





# represents strip between a < ℑ(z*exp(-im*θ)) < b
immutable Strip{T} <: Domain{Complex{T},2}
    a::T
    b::T
    θ::T
end

Strip(a,b,θ) = Strip(promote(a,b,θ)...)
Strip(a,b) = Strip(a,b,0)
Strip(b) = Strip(-b,b)

Base.show(io::IO,d::ComplexPlane)=print(io,"ℂ")
Base.show(io::IO,d::Complement)=print(io,"$(d.full) \\$(d.del)")



end #module




