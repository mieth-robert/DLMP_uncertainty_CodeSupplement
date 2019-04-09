   
function build_model(feeder::FeederTopo, settings::Dict)

    # Paramters for thermal line flow approximation
    # See "Distributed Generation Hosting Capacicty Evaluation for Distribution Systems Considering the Robust Optimal Operation of OLTC and SVC"
    a1 = [1, 1, 0.2679, -0.2679, -1, -1, -1, -1, -0.2679, 0.2679, 1, 1]
    a2 = [0.2679, 1, 1, 1, 1, 0.2679, -0.2679, -1, -1, -1, -1, -0.2679]
    a3 = -1 .* [-1, -1.366, -1, -1, -1.366, -1, -1, -1.366, -1, -1, -1.366, -1]

    # Get and prepare feeder data 
    buses = feeder.buses
    lines = feeder.lines
    generators = feeder.generators
    n_buses = feeder.n_buses
    root_bus = feeder.root_bus
    gen_buses = feeder.gen_buses
    line_to = feeder.line_to    
    v_root = 1

    Rd = feeder.R
    A = feeder.A[1:end, 2:end]
    A_check = A^-1
    R = A'*Rd*A
    R_check = R^(-1)

    e = ones(n_buses-1)

    bus_set = collect(1:n_buses)
    non_root_buses = setdiff(bus_set, [root_bus])

    # Get and prepare settings
    var_vec = settings["var_vec"]
    Σ = settings["Σ"]
    z_g = settings["z_g"]
    z_v = settings["z_v"]
    toggle_volt_cc = settings["toggle_volt_cc"]
    toggle_gen_cc = settings["toggle_gen_cc"]
    toggle_thermal_cc = settings["toggle_thermal_cc"]
    thermal_const_method = settings["thermal_const_method"]
    vfac = settings["vfac"]
    qcfac = settings["qcfac"]
    output_level = settings["output_level"]
    Ψ = settings["Ψ"]
    if "loadfac" in keys(settings)
        loadfac = settings["loadfac"]
        buses = change_load_same_pf(buses, loadfac)
    end


    # Build cost vector and matrices
    c = [] # linear cost vector
    C = zeros(n_buses,n_buses) # quadratic cost matrix
    for i in 1:n_buses
        if i in gen_buses
            c_i = buses[i].generator.cost
            C[i,i] = buses[i].generator.quad_cost * qcfac
        else
            c_i = 0
        end
        push!(c, c_i)
    end
    F = C^(1/2)

    # Prepare cc matrices
    Σ_rt = Σ[2:end, 2:end]^(1/2)
    s = sqrt(sum(Σ_rt))

    # Build model
    any_cc = toggle_volt_cc || toggle_gen_cc || toggle_thermal_cc

    m = Model(with_optimizer(Mosek.Optimizer, MSK_IPAR_LOG=output_level))

    # Standard Constraints 
    @variable(m, v[bus_set] >=0) # voltage square
    @variable(m, fp[bus_set]) # active power flow
    @variable(m, fq[bus_set]) # reactive power flow
    @variable(m, gp[bus_set]) # active power generation
    @variable(m, gq[bus_set]) # reactive power generation
    @variable(m, r_sched >= 0) # quadratic part of cost function 

    if any_cc
        @variable(m, α[bus_set] >=0) # Balancing Participation factor
        @variable(m, r_bal >= 0) #   
    end
    if toggle_volt_cc
        # additional variables needed for voltage soc reformulation
        n = n_buses-1
        @variable(m, ρ[1:n]) 
        @variable(m, t[1:n] >=0)
    end
    if thermal_const_method == 2
        # additional variable for thermal line flow soc reformulation
        n = n_buses-1
        @variable(m, ρ_f[1:n]) 
        @variable(m, t_f[1:n] >=0)
    end


    # Energy Balances
    @constraint(m, λ[b=bus_set], buses[b].d_P - gp[b] + sum(fp[k] for k in buses[b].children) == fp[b])
    @constraint(m, π[b=bus_set], buses[b].d_Q - gq[b] + sum(fq[k] for k in buses[b].children) == fq[b])

    non_root_buses = setdiff(bus_set, [root_bus])

    # Deterministic voltage equations
    @constraint(m, β[b=non_root_buses], v[b] == v[buses[b].ancestor[1]] - 2*(line_to[b].r * fp[b] + line_to[b].x * fq[b]))

    # Substation Constraints
    @constraint(m, v[root_bus] == v_root)
    @constraint(m, fp[root_bus] == 0)
    @constraint(m, fq[root_bus] == 0)

    # Non-generating bus constraints
    buses_without_generation = setdiff(bus_set, gen_buses)
    # buses_without_generation = non_root_buses
    @constraint(m, [b=buses_without_generation], gp[b] == 0)
    @constraint(m, [b=buses_without_generation], gq[b] == 0)

    if any_cc
        # Balancing of participation
        @constraint(m, γ, sum(α) == 1)
        @constraint(m, [b=buses_without_generation], α[b] == 0)
    end

    # Deterministic Voltage Constraints
    if  toggle_volt_cc
        # Voltage Chance Constraints
        eΣ_rt = Array(e'*Σ_rt)
        RΣ_rt = Array(R*Σ_rt)
        soc_vectors = []
        idx_to_bus = Dict()
        bus_to_idx = Dict()
        for (i, b) in enumerate(non_root_buses)
            idx_to_bus[i] = b
            bus_to_idx[b] = i
            y = RΣ_rt[i,:] - ρ[i].*eΣ_rt'
            soc = vcat(t[i], y)
            soc = vec(soc)
            push!(soc_vectors, soc)
        end
        # NOTE: indices of constraints refer to non-root indices 
        @constraint(m, ζ[i=1:n], soc_vectors[i] in SecondOrderCone())
        @constraint(m, η[i=1:n], sum(R_check[i,ii] * ρ[ii] for ii in 1:n) == α[i])
        # @constraint(m, η[i=1:n], sum(R[i,ii] * α[ii] for ii in 1:n) == ρ[i])
        @constraint(m, μp[b=non_root_buses], v[b] + 2*z_v*t[bus_to_idx[b]] <= (vfac > 0 ? (1+vfac)^2 : buses[b].v_max))
        @constraint(m, μm[b=non_root_buses], -v[b] + 2*z_v*t[bus_to_idx[b]] <= -(vfac > 0 ? (1-vfac)^2 : buses[b].v_min))
    else    
        idx_to_bus = collect(1:n_buses)
        bus_to_idx = collect(1:n_buses)
        @constraint(m, μp[b=non_root_buses], v[b] <= (vfac > 0 ? (1+vfac)^2 : buses[b].v_max))
        @constraint(m, μm[b=non_root_buses], v[b] >= (vfac > 0 ? (1-vfac)^2 : buses[b].v_min))
    end


    if toggle_gen_cc
        # Generation Chance Constraints    
        @constraint(m, δp[b=setdiff(gen_buses,[root_bus])], gp[b] + α[b]*z_g*s <= buses[b].generator.g_P_max)
        @constraint(m, δm[b=setdiff(gen_buses,[root_bus])], gp[b] - α[b]*z_g*s >= 0)
    #     @constraint(m, gp[root_bus] - α[root_bus]*z*s >= 0)
    else
        # Deterministic constraints on active power
        @constraint(m, δp[b=setdiff(gen_buses,[root_bus])], gp[b] <= buses[b].generator.g_P_max)
        @constraint(m, δm[b=setdiff(gen_buses,[root_bus])], gp[b] >= 0)
    #     @constraint(m, gp[root_bus] >= 0)
    end
    # Reactive Power Constraints
    @constraint(m, θp[b=gen_buses], gq[b] <= buses[b].generator.g_Q_max)
    @constraint(m, θm[b=gen_buses], gq[b] >= -buses[b].generator.g_Q_max)

    # constaints on thermal line capacity 
    if toggle_thermal_cc
        # Line Limit Chance Concstraints
        eΣ_rt = Array(e'*Σ_rt)
        AΣ_rt = Array(A*Σ_rt)
        soc_vectors_f = []
        idx_to_bus = Dict()
        bus_to_idx = Dict()
        for (i, b) in enumerate(non_root_buses)
            idx_to_bus[i] = b
            bus_to_idx[b] = i
            y_f = AΣ_rt[i,:] - ρ_f[i].*eΣ_rt'
            soc = vcat(t_f[i], y_f)
            soc = vec(soc)
            push!(soc_vectors_f, soc)
        end
        @constraint(m, ζ_f[i=1:n], soc_vectors_f[i] in SecondOrderCone())
        @constraint(m, η_f[i=1:n], sum(A_check[i,ii] * ρ_f[ii] for ii in 1:n) == α[i])
        @constraint(m, τ[b=setdiff(bus_set,[root_bus]), c=1:12], a1[c]*(fp[b] + t_f[bus_to_idx[b]]) + a2[c]*fq[b] <= a3[c]*line_to[b].s_max)
    else
        if thermal_const_method == 1
            @constraint(m, τ[b=setdiff(bus_set,[root_bus])], [line_to[b].s_max, fp[b], fq[b]] in SecondOrderCone())
        elseif thermal_const_method == 2
            @constraint(m, τ[b=setdiff(bus_set,[root_bus]), c=1:12], a1[c]*fp[b] + a2[c]*fq[b] <= a3[c]*line_to[b].s_max)
        else
            @warn("Thermal constraint method $(thermal_const_method) unknown. Proceeding with unconstrained lines.")
        end
    end


    # Linear Part of objective
    @expression(m, linear_cost, sum(gp[b]*c[b] for b in gen_buses))
    
    # Quadratic part of objective in as soc (see note below)
    # No quadratic cost on substation
    # quadratic cost term of substation is penalty for alpha to reduce deviation from schedule
    F_schedule = copy(F)
    F_schedule[root_bus, root_bus] = 0
    if any_cc
        cost_soc = vcat(r_sched, F_schedule'*gp)
        bal_soc = vcat(r_bal, F'*α)
        @constraint(m, bal_soc in SecondOrderCone())
        @expression(m, quad_cost, r_sched + r_bal)
    else
        cost_soc = vcat(r_sched, s.*F_schedule'*gp)
        @expression(m, quad_cost, r_sched)
    end
    @constraint(m, cost_soc in SecondOrderCone())

    # Variance penalty
    if toggle_volt_cc
        @expression(m, variance_penalty, Ψ*sum(t))
    else
        @expression(m, variance_penalty, 0)
    end

    @objective(m, Min, linear_cost + quad_cost + variance_penalty)

    # Some info to send to results handling, not super elegant but I leave it for now
    meta = Dict(
        "idx_to_bus" => idx_to_bus,
        "bus_to_idx" => bus_to_idx,
        "toggle_volt_cc" => toggle_volt_cc,
        "toggle_gen_cc" => toggle_gen_cc,
        "any_cc" => any_cc,
        "z_v" => z_v,
        "s" => s,
        )

    return m, meta

end




# Note on quadratic objective:
# The program is a conic program with quadratic objective. Usually solvers (like Mosek) do not allow a 
# mixture of those two approaches. Therefore the quadratic objective has been reformulated as a second order 
# cone constraint with linear objective rendering the whole problem in a SOC constrained problem with linear objective.
# The original quadratic objective value can be recovered from r as:
# corrected objective = linear_cost + 