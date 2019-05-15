# Copyright (c) 2019 Robert Mieth
# Code supplement to the paper "Distribution Electricity Pricing under Uncertainty" by Robert Mieth and Yury Dvorkin
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# +++++
# output.jl
#
# 
# Provides a functin to process the results of the solved and returns a comprehesive DataFrame
# +++++
# Devnotes: ---



function results_to_df(m, meta, feeder)

    a1 = [1, 1, 0.2679, -0.2679, -1, -1, -1, -1, -0.2679, 0.2679, 1, 1]
    a2 = [0.2679, 1, 1, 1, 1, 0.2679, -0.2679, -1, -1, -1, -1, -0.2679]
    a3 = -1 .* [-1, -1.366, -1, -1, -1.366, -1, -1, -1.366, -1, -1, -1.366, -1]

    n_buses = feeder.n_buses
    root_bus = feeder.root_bus
    gen_buses = feeder.gen_buses
    bus_set = collect(1:n_buses)
    non_root_buses = setdiff(bus_set,[root_bus])

    any_cc = meta["any_cc"]
    toggle_volt_cc = meta["toggle_volt_cc"]
    idx_to_bus = meta["idx_to_bus"]
    bus_to_idx = meta["bus_to_idx"]
    toggle_thermal_cc = meta["toggle_thermal_cc"]
    thermal_const_method = meta["thermal_const_method"]

    Rd = feeder.R
    A = feeder.A[1:end, 2:end]
    R = A'*Rd*A
    R_check = R^(-1)

    # Results from calculations
    b_idx = []
    gp_res = []
    gq_res = []
    fp_res = []
    fq_res = []
    voltages = []
    alphas = []
    lambdas = []
    pies = []
    delta_plus, delta_minus = [], []
    mu_plus, mu_minus = [], []
    eta_plus, eta_minus = [], []
    eta_aP = []
    eta_aQ = []
    gamma = []
    rhos = []
    nus = []
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
            push!(fp_res, 0)
            push!(fq_res, 0)
            push!(eta_plus, 0)
            push!(eta_minus, 0)
            push!(eta_aP, 0)
            push!(eta_aQ, 0)
        else
            push!(mu_plus, shadow_price(m[:μp][b]))
            push!(mu_minus, shadow_price(m[:μm][b]))
            push!(fp_res, value(m[:fp][b]))
            push!(fq_res, value(m[:fq][b]))

            if toggle_thermal_cc
                push!(eta_plus, sum(shadow_price(m[:ηp][b, c]) for c in 1:12))
                push!(eta_minus, sum(shadow_price(m[:ηm][b, c]) for c in 1:12))
                push!(eta_aP, 0)
                push!(eta_aQ, 0)
            else
                if thermal_const_method == 2
                    push!(eta_plus, sum(shadow_price(m[:η][b, c]) for c in 1:12))
                    push!(eta_aP, sum(shadow_price(m[:η][b, c])*a1[c] for c in 1:12))
                    push!(eta_aQ, sum(shadow_price(m[:η][b, c])*a2[c] for c in 1:12))
                    push!(eta_minus, 0) 
                elseif thermal_const_method == 1
                    push!(eta_plus, "na")
                    push!(eta_minus, 0) 
                    push!(eta_aP, 0)
                    push!(eta_aQ, 0)
                else
                    push!(eta_plus, 0)
                    push!(eta_minus, 0)
                    push!(eta_aP, 0)
                    push!(eta_aQ, 0)
                end
            end

        end

        if !toggle_volt_cc || b == root_bus
            push!(voltvar, 0)
            push!(nus, 0)
            push!(rhos, 0)
        else
            push!(voltvar, value(m[:t][bus_to_idx[b]])^2)
            push!(nus, shadow_price(m[:ν][bus_to_idx[b]]))
            push!(rhos, value(m[:ρ][bus_to_idx[b]]))
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
        fp = fp_res,
        fq = fq_res,
        alpha = alphas,
        lambda = lambdas,
        pi = pies,
        gamma = gamma,
        rho = rhos,
        nu = nus,
        delta_plus = delta_plus,
        delta_minus = delta_minus,
        mu_plus = mu_plus,
        mu_minus = mu_minus,
        eta_plus = eta_plus,
        eta_minus = eta_minus,
        eta_aP = eta_aP, 
        eta_aQ = eta_aQ, 
        voltvar = voltvar,
    )


    # DLMP Decomposition
    lambda_a = []
    rx_pi = []
    rx_pi_a = []
    r_sum_mu_d = []
    nu_calc = []
    gamma_calc = []
    rx_etaQ = []
    for bus in feeder.buses
        if bus.is_root
            push!(lambda_a, 0)
            push!(rx_pi, 0)
            push!(rx_pi_a, 0)
            push!(r_sum_mu_d, 0)
            push!(nu_calc, 0)
            push!(rx_etaQ, 0)
        else
            anc = bus.ancestor
            cs = bus.children
            push!(lambda_a, results_df[results_df[:bus] .== anc, :lambda][1])
            
            pi_i =  results_df[results_df[:bus] .== bus.index, :pi][1]
            push!(rx_pi, pi_i * (feeder.line_to[bus.index].r / feeder.line_to[bus.index].x))
            
            pi_a = results_df[results_df[:bus] .== anc, :pi][1]
            push!(rx_pi_a, pi_a * (feeder.line_to[bus.index].r / feeder.line_to[bus.index].x))

            downstream = traverse(feeder, bus.index)
            v = 2 * feeder.line_to[bus.index].r * sum((results_df[results_df[:bus] .== d, :mu_plus][1] - results_df[results_df[:bus] .== d, :mu_minus][1]) for d in downstream)
            push!(r_sum_mu_d, v)

            rxeQ = results_df[results_df[:bus] .== bus.index, :eta_aQ][1] * (feeder.line_to[bus.index].r / feeder.line_to[bus.index].x)
            push!(rx_etaQ, rxeQ)

            if toggle_volt_cc
                nu = 2 * meta["z_v"] * sum(R[bus_to_idx[j], bus_to_idx[bus.index]] * (results_df[results_df[:bus] .== j, :mu_plus][1] + results_df[results_df[:bus] .== j, :mu_minus][1]) * ((sum(R[bus_to_idx[j], bus_to_idx[k]] * (meta["Σ"][k,k] + results_df[results_df[:bus] .== k, :alpha][1] * meta["s"]^2)  for k in non_root_buses)) / sqrt(results_df[results_df[:bus] .== j, :voltvar][1]))  for j in non_root_buses)
            else
                nu = 0
            end
            push!(nu_calc, nu)
        end
    end
  
    results_df[:lambda_anc] = lambda_a
    results_df[:rx_pi_i] = rx_pi
    results_df[:rx_pi_a] = rx_pi_a
    results_df[:rx_etaQ] = rx_etaQ
    results_df[:r_sum_mu_d] = r_sum_mu_d
    results_df[:nu_calc] = nu_calc


    return results_df
end