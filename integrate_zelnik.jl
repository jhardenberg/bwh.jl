# Zelnik et al. equations
# Integrate equation using forward Euler

@parallel function update_b!(b2::Data.Array, b::Data.Array, w::Data.Array, p::Data.Number, η::Data.Number, λ::Data.Number, ρ::Data.Number, ν::Data.Number, db::Data.Number, dw::Data.Number, dt::Data.Number, dx::Data.Number, dy::Data.Number)
    # λwb(1-b)(1+ηb)^2 - b + db*∇^2(b)
     @inn(b2) = @inn(b) + dt*( λ*@inn(w)*@inn(b)*(1-@inn(b))*(1+η*@inn(b))^2 - @inn(b) + db*(@d2_xi(b)/dx^2 + @d2_yi(b)/dy^2) );
     return
 end
 
 @parallel function update_w!(w2::Data.Array, b::Data.Array, w::Data.Array, p::Data.Number, η::Data.Number, λ::Data.Number, ρ::Data.Number, ν::Data.Number, db::Data.Number, dw::Data.Number, dt::Data.Number, dx::Data.Number, dy::Data.Number)
     # p - νw(1-ρb)-λwb(1+ηb)^2 + dw*∇^2(w)
     @inn(w2) = @inn(w) + dt*( p - ν*@inn(w)*(1-ρ*@inn(b)) - λ*@inn(w)*@inn(b)*(1+η*@inn(b))^2 + dw*(@d2_xi(w)/dx^2 + @d2_yi(w)/dy^2) );
     return
 end
