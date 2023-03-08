using AlgebraicPetri

using OrdinaryDiffEq
using Plots

using Catlab

@theory Group begin
    @op begin
        (⋅) := multiply
    end

    El::TYPE

    inv(A::El)::El
    id::El

    id ⋅ a == a ⊣ (a::El)
    a ⋅ id == a ⊣ (a::El)
    a ⋅ inv(a) == id ⊣ (a::El)
    inv(a) ⋅ a == id ⊣ (a::El)
    a ⋅ (b ⋅ c) == (a ⋅ b) ⋅ c ⊣ (a::El, B::El, C::El)
end

############ SIR MODEL ############





############ Lotka-Volterra ############
