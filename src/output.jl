
function results_to_df(m, meta, feeder)

    n_buses = feeder.n_buses
    root_bus = feeder.root_bus
    gen_buses = feeder.gen_buses
    bus_set = collect(1:n_buses)

    any_cc = meta["any_cc"]
    toggle_volt_cc = meta["toggle_volt_cc"]
    idx_to_bus = meta["idx_to_bus"]
    bus_to_idx = meta["bus_to_idx"]


    # Results from calculations
    b_idx = []
    gp_res = []
    gq_res = []
    voltages = []
    alphas = []
    lambdas = []
    pies = []
    delta_plus, delta_minus = [], []
    mu_plus, mu_minus = [], []
    gamma = []
    etas = []
    voltvar = []


    num_approx = 1e-7

    for b in bus_set
        push!(b_idx, b) 
        push!(gp_res, (abs(value(m[:gp][b])) - num_approx) > 0 ? value(m[:gp][b]) : 0 )
        push!(gq_res, (abs(value(m[:gq][b])) - num_approx) > 0 ? value(m[:gq][b]) : 0 )
        push!(voltages, sqrt(value(m[:v][b])))
        if any_cc 
            push!(alphas, (abs(value(m[:α][b])) - num_approx) > 0 ? value(m[:α][b]) : 0 ) 
        else
            push!(alphas, 0)
        end
        push!(lambdas, shadow_price(m[:λ][b]))
        push!(pies, shadow_price(m[:π][b]))
        if b in setdiff(gen_buses,[root_bus])
            push!(delta_plus, shadow_price(m[:δp][b]))
            push!(delta_minus, shadow_price(m[:δm][b]))
        else
            push!(delta_plus, 0)
            push!(delta_minus, 0)
        end
        
        if b == root_bus
            push!(mu_plus, 0)
            push!(mu_minus, 0)
        else
            push!(mu_plus, shadow_price(m[:μp][b]))
            push!(mu_minus, shadow_price(m[:μm][b]))
        end
        if !toggle_volt_cc || b == root_bus
            push!(voltvar, 0)
            push!(etas, 0)
        else
            push!(voltvar, value(m[:t][bus_to_idx[b]]))
            push!(etas, shadow_price(m[:η][bus_to_idx[b]]))
        end
    end

    objective = zeros(n_buses)
    objective[root_bus] = objective_value(m)
    objective_corrected = objective[root_bus]
    gamma = zeros(n_buses)
    if any_cc
        gamma[root_bus] = shadow_price(m[:γ])
    end


    results_df = DataFrame(
        objective = objective,
        bus = b_idx, 
        gp = gp_res,
        gq = gq_res,
        voltage = voltages,
        alpha = alphas,
        lambda = lambdas,
        pi = pies,
        gamma = gamma,
        eta = etas,
        delta_plus = delta_plus,
        delta_minus = delta_minus,
        mu_plus = mu_plus,
        mu_minus = mu_minus,
        voltvar = voltvar,
    )


    # DLMP Decomposition
    lambda_a = []
    rx_pi = []
    rx_pi_a = []
    r_sum_mu_d = []
    for bus in feeder.buses
        if bus.is_root
            push!(lambda_a, 0)
            push!(rx_pi, 0)
            push!(rx_pi_a, 0)
            push!(r_sum_mu_d, 0)
        else
            anc = bus.ancestor
            cs = bus.children
            push!(lambda_a, results_df[results_df[:bus] .== anc, :lambda][1])
            
            pi_i =  results_df[results_df[:bus] .== bus.index, :pi][1]
            push!(rx_pi, pi_i * (feeder.line_to[bus.index].r / feeder.line_to[bus.index].x))
            
            pi_a = pi_i =  results_df[results_df[:bus] .== anc, :pi][1]
            push!(rx_pi_a, pi_a * (feeder.line_to[bus.index].r / feeder.line_to[bus.index].x))

            downstream = traverse(feeder, bus.index)
            v = 2 * feeder.line_to[bus.index].r * sum((results_df[results_df[:bus] .== d, :mu_plus][1] - results_df[results_df[:bus] .== d, :mu_minus][1]) for d in downstream)
            push!(r_sum_mu_d, v)
        end
    end
  
    results_df[:lambda_anc] = lambda_a
    results_df[:rx_pi_i] = rx_pi
    results_df[:rx_pi_a] = rx_pi_a
    results_df[:r_sum_mu_d] = r_sum_mu_d


    return results_df
end