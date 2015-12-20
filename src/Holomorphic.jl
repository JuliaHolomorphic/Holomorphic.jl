__precompile__()
module Holomorphic
    using Base, ApproxFun, SingularIntegralEquations


import ApproxFun: UnivariateSpace, domain, evaluate, ComplexBasis, spacescompatible,
                    Space, defaultFun
import SingularIntegralEquations: stieltjes

# represents the whole Planess
immutable ComplexPlane <: Domain{Complex128,2} end

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

domain(sp::StieltjesSpace)=Complement(ComplexPlane(),domain(sp.space))


evaluate(v::AbstractVector,sp::StieltjesSpace,z)=stieltjes(sp.space,v,z)

stieltjes(f::Fun)=Fun(f.coefficients,StieltjesSpace(space(f)))

# construct spaces for complement
# we assume boundedness

Space{DD<:Interval}(Γ::Complement{ComplexPlane,DD})=StieltjesSpace(JacobiWeight(0.5,0.5,Γ.del))
Space{DD<:UnionDomain}(Γ::Complement{ComplexPlane,DD})=mapreduce(d->StieltjesSpace(d).space,PiecewiseSpace,Γ.del.domains)




# represents strip between a < im z < b
# TODO: rotations
immutable Strip <: Domain{Complex128,2}
    a::Float64
    b::Float64
end




end #module



