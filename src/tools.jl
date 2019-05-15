# Copyright (c) 2019 Robert Mieth
# Code supplement to the paper "Distribution Electricity Pricing under Uncertainty" by Robert Mieth and Yury Dvorkin
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# +++++
# tools.jl
#
# Provides auxiallary functions
# +++++
# Devnotes: ---



# Build, run and process a single model run
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


# Build, run and process a multiple runs from the experiment dict
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
   

# Save all results of the experiment
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

