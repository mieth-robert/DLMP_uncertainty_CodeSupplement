# Copyright (c) 2019 Robert Mieth
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# +++++
# input.jl
#
# Definition of typed to handle network data
# Functions to read and prepare data and settings
# +++++
# Devnotes: ---


# TYPE DEFINITIONS
mutable struct Generator
   index::Any
   bus_idx::Int
   g_P_max::Float64
   g_Q_max::Float64
   cost::Float64
   quad_cost::Float64
   function Generator(index, bus_idx, g_P_max, g_Q_max)
      g = new()
      g.index  = index
      g.bus_idx = bus_idx
      g.cost = 1.00
      g.quad_cost = 0.00
      g.g_P_max = g_P_max
      g.g_Q_max = g_Q_max
      return g
   end
end

mutable struct Bus
   index::Any
   is_root::Bool
   d_P::Float64
   d_Q::Float64
   cosphi::Float64
   tanphi::Float64
   v_max::Float64
   v_min::Float64
   children::Vector{Int}
   ancestor::Vector{Int}
   generator::Generator
   beta_0::Float64
   beta_1::Float64
   function Bus(index, d_P, d_Q, v_max, v_min)
      b = new()
      b.index = index
      b.is_root = false
      b.d_P = d_P
      b.d_Q = d_Q
      b.v_max = v_max
      b.v_min = v_min
      b.children = Int[]
      b.ancestor = Int[]
      cosphi = d_P/(sqrt(d_P^2 + d_Q^2))
      tanphi = d_Q/d_P
      if isnan(cosphi)
        b.cosphi = 0
        b.tanphi = 0
      else
        b.cosphi = cosphi
        b.tanphi = tan(acos(cosphi))
      end
      return b
   end
end
function set_bus_active_load(b::Bus, dP; auto_set_dQ=false)
    b.d_P = dP
    if auto_set_dQ
        b.d_Q = dP * b.tanphi
    end
end

mutable struct Line
   index::Any
   to_node::Int # the "to" node
   from_node::Int # the "from" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   b::Float64 # the susceptance value
   s_max::Float64 # the capacity of the line
   function Line(index, to_node, from_node, r, x, s_max)
      l = new()
      l.index = index
      l.to_node = to_node
      l.from_node = from_node
      l.r = r
      l.x = x
      l.b = (x/(r^2 + x^2))
      l.s_max = s_max
      return l
   end
end

mutable struct FeederTopo
  label::Any
  buses::Array{Bus}
  lines::Array{Line}
  generators::Array{Generator}
  n_buses::Int
  gen_buses::Array{Int}
  dr_buses::Array{Int}
  line_to::Dict{}
  A::Array{Float64, 2}
  R::Array{Float64, 2}
  X::Array{Float64, 2}
  γ::Array{Float64, 2}
  root_bus::Int
  function FeederTopo(buses, lines, generators; label="RadFeeder")
    f  = new()
    f.label = label
    f.buses = buses
    f.lines = lines
    f.generators = generators
    f.n_buses = length(buses)
    f.gen_buses = [g.bus_idx for g in generators]
    for (i,b) in enumerate(buses)
      if b.is_root
        f.root_bus = i
        break
      end
    end
    line_to = Dict()
    for l in lines
      line_to[l.to_node] = l
    end
    f.line_to = line_to
    f.A = zeros(length(lines), length(buses))
    for b in 1:f.n_buses
      a = b
      while a != f.root_bus
        f.A[a-1, b] = 1
        a = buses[a].ancestor[1]
      end
    end
    γ_vec = [buses[b].tanphi for b in 1:f.n_buses]
    f.γ = diagm(0 => γ_vec)
    R_vec = [f.line_to[b].r for b in setdiff(collect(1:f.n_buses), [f.root_bus])]
    X_vec = [f.line_to[b].x for b in setdiff(collect(1:f.n_buses), [f.root_bus])]
    f.R = diagm(0 => R_vec)
    f.X = diagm(0 => X_vec)
    return f
  end
end


function set_root_feeder_cost(feeder, cost)
  r = feeder.root_bus
  if r in feeder.gen_buses
    feeder.buses[r].generator.cost = cost
  else
    println("!!! Could not set root feeder cost. No generator at root")
  end
end


function load_feeder(datadir)
  # READ RAW DATA
  println(">>>>> Reading feeder data from $(datadir)")

  nodes_raw = CSV.read("$datadir/nodes.csv")
  sum(nonunique(nodes_raw, :index)) != 0 ? @warn("Ambiguous Node Indices") : nothing

  lines_raw = CSV.read("$datadir/lines.csv")
  sum(nonunique(lines_raw, :index)) != 0  ? @warn("Ambiguous Line Indices") : nothing

  generators_raw = CSV.read("$datadir/generators.csv")
  sum(nonunique(generators_raw, :index)) != 0 ? @warn("Ambiguous Generator Indices") : nothing

  # PREPARING MODEL DATA
  buses = []
  for n in 1:nrow(nodes_raw)
      index = nodes_raw[n, :index]
      d_P = nodes_raw[n, :d_P]
      d_Q = nodes_raw[n, :d_Q]
      v_max = nodes_raw[n, :v_max]
      v_min = nodes_raw[n, :v_min]
      newb = Bus(index, d_P, d_Q, v_max, v_min)
      in(:beta_0, names(nodes_raw)) ? newb.beta_0 = nodes_raw[n, :beta_0] : nothing
      in(:beta_1, names(nodes_raw)) ? newb.beta_1 = nodes_raw[n, :beta_1] : nothing
      push!(buses, newb)
  end

  lines = []
  for l in 1:nrow(lines_raw)
      index = lines_raw[l, :index]
      from_node = lines_raw[l, :from_node]
      to_node = lines_raw[l, :to_node]
      r = lines_raw[l, :r]
      x = lines_raw[l, :x]
      b = lines_raw[l, :b]
      s_max = lines_raw[l, :s_max]
      newl = Line(index, to_node, from_node, r, x, s_max)
      push!(buses[newl.from_node].children, newl.to_node)
      push!(buses[newl.to_node].ancestor, newl.from_node)
      push!(lines, newl)
  end

  # Check topology
  r = 0
  root_bus = 0
  for (i,b) in enumerate(buses)
      l = length(b.ancestor)
      if l > 1
          @warn("Network not Radial (Bus $(b.index))")
      elseif l == 0
          b.is_root = true
          root_bus = i
          r += 1
      end
  end
  if r == 0
      @warn("No root detected")
      root_bus = 0
  elseif r > 1
      @warn("More than one root detected")
  end

  generators = []
  for g in 1:nrow(generators_raw)
      index = generators_raw[g, :index]
      bus_idx = generators_raw[g, :node]
      g_P_max = generators_raw[g, :p_max]
      g_Q_max = generators_raw[g, :q_max]
      cost = generators_raw[g, :cost]
      quad_cost = generators_raw[g, :quad_cost]
      newg = Generator(index, bus_idx, g_P_max, g_Q_max)
      newg.cost = cost
      newg.quad_cost = quad_cost
      buses[newg.bus_idx].generator = newg
      push!(generators, newg)
  end

  feeder = FeederTopo(buses, lines, generators)

  return feeder
end

function read_price_data(filename)
   data_df = CSV.read(filename, header=3)
   prices = data_df[:price]
   timestamps = data_df[:t]
   return prices, timestamps
end


function change_load_same_pf(buses_in, α)
# always uses the GLOBAL BUSES as reference
   buses = deepcopy(buses_in)
   for i in keys(buses)
      P = buses[i].d_P
      Q = buses[i].d_Q
      PF = P/sqrt(P^2 + Q^2)
      if !(isnan(PF))
         buses[i].d_P = P*α
         buses[i].d_Q = P*α*sqrt((1-PF^2)/PF^2)
      end
   end

   return buses
end

function add_generator(index, bus, g_P_max, g_Q_max, cost)
    newg = Generator(index, bus, g_P_max, g_Q_max)
    newg.cost = cost
    BUSES[newg.bus_idx].generator = newg
    GENERATORS[newg.index] = newg
end


function load_timeseries(;datafolder="")
  info("Reading Timeseries Data from $(datafolder)")

  rel_demand_hourly_raw = CSV.read("data/$datafolder/forecast_hourly.csv", datarow=1)
  rel_demand_5min_raw = CSV.read("data/$datafolder/load_5min.csv", datarow=1)
  rel_lbmp_hourly_raw = CSV.read("data/$datafolder/prices_hourly.csv", datarow=1)

  rel_demand_hourly = rel_demand_hourly_raw[:Column1]
  rel_demand_5min = rel_demand_5min_raw[:Column1]
  rel_lbmp_hourly = rel_lbmp_hourly_raw[:Column1]


  return rel_demand_hourly, rel_demand_5min, rel_lbmp_hourly
end
