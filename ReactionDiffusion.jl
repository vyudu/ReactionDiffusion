### Model with three reactions
# B -> [] 
# [] -> A
# A + 2B <-> 3B

# 2F <-> B_2
# B_2 + F <-> B_3
# B_3 + F <-> B_4
# 2F <-> A_2
# A_2 + F <-> A_3
# A_3 + F <-> A_4

using Plots
using FFMPEG

struct ReactionDiffusion
    D_A::Float64
    D_B::Float64
    dt::Float64 #timestep

    f::Float64 #feed rate
    k::Float64 #kill rate
    r::Float64 #reaction rate
    size::Int64 #grid size

    grid_A::Matrix{Float64}
    grid_B::Matrix{Float64}
    grid_T::Matrix{Float64}
end

function evolve(model::ReactionDiffusion) 
    ΔA_d = model.dt * diffuse(model.grid_A, model.D_A)
    ΔB_d = model.dt * diffuse(model.grid_B, model.D_B)
    ΔA_r, ΔB_r = react(model)

    model.grid_A .= model.grid_A + ΔA_d + ΔA_r
    model.grid_B .= model.grid_B + ΔB_d + ΔB_r
end

# finite-differences: [n(x+1) - 2n(x) + n(x-1)]  
function diffuse(grid, d) 
    uxx = circshift(grid,(0,-1))+circshift(grid,(0,1))-2*grid
    uyy = circshift(grid,(-1,0))+circshift(grid,(1,0))-2*grid
    d*(uxx + uyy)
end

function soretDiffusion(χ_grid, T_grid, d, d_T)
    uxx = (circshift(χ_grid,(0,-1))+circshift(χ_grid,(0,1))-2*χ_grid)
    uyy = (circshift(χ_grid,(-1,0))+circshift(χ_grid,(1,0))-2*χ_grid)

    Txx = (circshift(T_grid,(0,-1))+circshift(T_grid,(0,1))-2*T_grid)
    Tyy = (circshift(T_grid,(-1,0))+circshift(T_grid,(1,0))-2*T_grid)

    d*(uxx + uyy) + d_T * (χ_grid .* (1 - χ_grid))*(Txx + Tyy)
    # d_T / d is derived in Shiliang's paper
end
    

function react(model) 
    Δ_f = model.f * (1 .- model.grid_A) 
    Δ_k = model.k * (model.grid_B)
    Δ_r = model.r * model.grid_A .* (model.grid_B).^2

    Δ_f - Δ_r, Δ_r - Δ_k
end

function simulate(D, dt, f, k, grid_size, steps; random_seed=false, seed_size=5) 
    mid = Int(ceil(grid_size/2))
    grid_A = ones(grid_size, grid_size)
    grid_B = zeros(grid_size, grid_size)

    if !random_seed
        grid_B[mid-seed_size:mid+seed_size,mid-seed_size:mid+seed_size] .= 1.0
    else 
        for _ in 1:10
            idx = rand(1:grid_size,2)
            grid_B[idx...]=1.0
        end
    end

    anim = Animation()
    model = ReactionDiffusion(2*D, D, dt, f, k, 1, grid_size, grid_A, grid_B)
    
    for i in 1:steps
        evolve(model)
        heatdata = model.grid_B ./ (model.grid_A + model.grid_B) 
        hm = heatmap(x=1:grid_size, y=1:grid_size, heatdata, c = cgrad(:thermal))
        frame(anim, hm)
    end
    anim
end

function diffusionAnim(D, dt, size, steps)
    anim = Animation()
    cells = zeros(size, size)
    mid = Int(ceil(size/2))
    cells[mid, mid] = 1.0
    for i in 1:steps
        cells += dt * diffuse(cells, D)
        hm = heatmap(x=1:size, y=1:size, cells, c=cgrad(:thermal))
        frame(anim, hm)
    end
    anim
end

anim = simulate(0.5, 0.2, 0.034, 0.095, 101, 500)
gif(anim, "react1.gif", fps=60)
