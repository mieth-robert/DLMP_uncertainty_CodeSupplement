function run_single_case(feeder, casesettings)
    println(">>>> Building Model")
    m, meta = build_model(feeder, casesettings)
    println(">>>> Running Model")
    solvetime = @elapsed optimize!(m)
    status = termination_status(m)
    println(">>>> Model finished with status $status in $solvetime seconds")

    println(">>>> Post-Processing")
    println("")
    return results_to_df(m, meta, feeder)
end


function run_experiment(experiment) 
    for case in keys(experiment)
        println("Running $(experiment[case]["verbose"])")
        println(">>>> Building Model")
        m, meta = build_model(experiment[case]["feeder"], experiment[case]["settings"])
        println(">>>> Running Model")
        solvetime = @elapsed optimize!(m)
        status = termination_status(m)
        println(">>>> Model finished with status $status in $solvetime seconds")

        println(">>>> Post-Processing")
        # experiment[case]["model"] = m
        experiment[case]["results"] = results_to_df(m, meta, experiment[case]["feeder"])
        println("")
    end 
    return experiment
end
   

function save_experiment(folder, experiment)
    timestamp = Dates.format(Dates.now(), "yymmdd_HHMM")
    for case in keys(experiment)
        println("Saving $(experiment[case]["verbose"])")
        
        mkpath("results/$(folder)/all/$(timestamp)/")
        mkpath("results/$(folder)/latest/")
        CSV.write("results/$(folder)/all/$(timestamp)/$(experiment[case]["verbose"]).csv", experiment[case]["results"])
        CSV.write("results/$(folder)/latest/$(experiment[case]["verbose"]).csv", experiment[case]["results"])

        save("results/$(folder)/all/$(timestamp)/experiment_dict.jld", "experiment", experiment)
    end
    println(">>>> Saved with timestamp $timestamp")
end

# Get all downstream nodes of bus_idx
function traverse(feeder, bus_idx; nodes=[])
    push!(nodes, bus_idx)
    if length(feeder.buses[bus_idx].children) > 0
        for c in feeder.buses[bus_idx].children
            nodes = traverse(feeder, c, nodes=nodes)
        end
    end
    return nodes 
end

# Recursively get the loss active and reactive demand loss factors 
# ∂l_i/∂d_k 
function get_lQi_dk(feeder, results, i, k)
    D = traverse(feeder, i)
    if k in D
        if feeder.buses[k].children == []
            fPi_dk = 1
            fQi_dk = 1
        else        
            fPi_dk = 1 + sum(get_lPi_dk(feeder, results, ii, k) for ii in setdiff(D, [i]))
            fQi_dk = 1 + sum(get_lQi_dk(feeder, results, ii, k) for ii in setdiff(D, [i]))
        end
    else
        fPi_dk = 0
        fQi_dk = 0
    end
    return 2*(results[:fp][i]*fPi_dk + results[:fq][i]*fQi_dk) * (feeder.line_to[i].x / (results[:voltage][i])^2)
end

function get_lPi_dk(feeder, results, i, k)
    D = traverse(feeder, i)
    if k in D
        if feeder.buses[k].children == []
            fPi_dk = 1
            fQi_dk = 1
        else        
            fPi_dk = 1 + sum(get_lPi_dk(feeder, results, ii, k) for ii in setdiff(D, [i]))
            fQi_dk = 1 + sum(get_lQi_dk(feeder, results, ii, k) for ii in setdiff(D, [i]))
        end
    else
        fPi_dk = 0
        fQi_dk = 0
    end
    return 2*(results[:fp][i] * fPi_dk + results[:fq][i] * fQi_dk) * (feeder.line_to[i].r / (results[:voltage][i])^2)
end
