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
        experiment[case]["model"] = m
        experiment[case]["results"] = results_to_df(m, meta, experiment[case]["feeder"])
        println("")
    end 
    return experiment
end
   

function save_experiment(folder, experiment)
    timestamp = Dates.format(Dates.now(), "yymmdd_HHMM")
    for case in keys(experiment)
        println("Saving $(experiment[case]["verbose"])")
        
        mkpath("results/$(folder)/$(timestamp)/")
        CSV.write("results/$(folder)/$(timestamp)/$(experiment[case]["verbose"]).csv", experiment[case]["results"])
        # save("results/$(folder)/$(timestamp)/experiment_dict.jld", "experiment", experiment)
    end
    println(">>>> Saved with timestamp $timestamp")
end

function traverse(feeder, bus_idx; nodes=[])
    push!(nodes, bus_idx)
    if length(feeder.buses[bus_idx].children) > 0
        for c in feeder.buses[bus_idx].children
            nodes = traverse(feeder, c, nodes=nodes)
        end
    end
    return nodes 
end