# Zelnik et al. equations
# Integrate equation using forward Euler

@parallel function update_b_clonal!(b2, b, w, η_c::Data.Number, γ_c::Data.Number, dt::Data.Number, dx::Data.Number, dy::Data.Number)
    # γ_cwbc(1-bc)(1+η_cbc)^2 - bc + 6*s_c*∇^2(bc)
     @inn(b2) = @inn(b) + dt*( γ_c*@inn(w)*@inn(b)*(1-@inn(b))*(1+η_c*@inn(b))^2 - @inn(b) + (@d2_xi(b)/dx^2 + @d2_yi(b)/dy^2) );
     return
 end

 @parallel function update_b_seeds!(b2, b, w, b_diff, η_s::Data.Number, γ_s::Data.Number, λ::Data.Number, μ::Data.Number, k::Data.Number, s_s::Data.Number, dt::Data.Number, dx::Data.Number, dy::Data.Number)
    # γ_s*λ*w*bs*(1-bs)(1+η_sbs)^2 - μ*bs + 6*s_s*∇^2(bs)
     @inn(b2) = @inn(b) + dt*( γ_s*λ*@inn(w)*@inn(b)*(1-@inn(b)/k)*(1+η_s*@inn(b))^2 - μ*@inn(b) + s_s*(@inn(b_diff)-@inn(b)));
     return
 end
 
 @parallel function update_w!(w2::Data.Array, b_c::Data.Array, b_s::Data.Array, w::Data.Array, p::Data.Number, γ_c::Data.Number, γ_s::Data.Number, η_c::Data.Number, η_s::Data.Number, ρ::Data.Number, ν::Data.Number, dw::Data.Number, dt::Data.Number, dx::Data.Number, dy::Data.Number)
     # p - νw/(1+ρ(bs+bc))-γ_cwbc(1+η_cbc)^2 -γ_swbs(1+η_sbs)^2 + dw*∇^2(w)
     @inn(w2) = @inn(w) + dt*( p - ν*@inn(w)/(1+ρ*(@inn(b_c)+@inn(b_s))) - γ_c*@inn(w)*@inn(b_c)*(1+η_c*@inn(b_c))^2 - γ_s*@inn(w)*@inn(b_s)*(1+η_s*@inn(b_s))^2 + dw*(@d2_xi(w)/dx^2 + @d2_yi(w)/dy^2) );
     return
 end