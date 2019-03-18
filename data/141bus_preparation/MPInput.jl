using MatpowerCases
using DataFrames

type Bus
   nodeID::Int
   root::Int
   Pd::Float64
   Qd::Float64
   Vmax::Float64
   Vmin::Float64
   B::Float64
   R::Float64
   X::Float64
   children::Vector{Int}
   ancestor::Vector{Int}
   genids::Int
   function Bus(nodeID, root, Pd, Qd, B, R, X, Vmax, Vmin)
      b = new(nodeID, root, Pd, Qd)
      b.Vmax = Vmax
      b.Vmin = Vmin
      b.R = R
      b.X = X
      b.B = B
      b.children = Int[]
      b.ancestor = Int[]
      b.genids = 0
      return b
   end
end

function setg(b::Bus, genidx)
   b.genids = genidx
end


#################################################################
type Generator
   genID::Int
   busidx::Int
   Pgmax::Float64
   Pgmin::Float64
   Qgmax::Float64
   Qgmin::Float64
   cost::Float64
   function Generator(genID, busidx, Pgmax, Pgmin, Qgmax, Qgmin)
      g = new()
      g.genID  = genID
      g.busidx = busidx
      g.cost = 1.00
      g.Pgmax = Pgmax
      g.Pgmin = Pgmin
      g.Qgmax = Qgmax
      g.Qgmin = Qgmin
      return g
   end
end
##################################################################
type Line
   arcID::Int
   tail::Int # the "to" node
   head::Int # the "from" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   u::Float64 # the capacity of the line
   function Line(arcID, tail, head, r, x, u)
      line = new(arcID, tail, head, r, x)
      line.u = u
      return line
   end
end
#########################################

function mpc_to_mat(mpc)
   # Branch Data
   branchdat = zeros(countnz(mpc["branch"][:, 11]), size(mpc["branch"], 1))
   for row in 1:size(mpc["branch"], 1)
      if mpc["branch"][row, 11] > 0
         for col in 1:size(mpc["branch"], 2)
            branchdat[row, col] = mpc["branch"][row, col]
         end
      end
   end
   branchmat = DataFrame()
   branchmat[:to] = Int.(branchdat[:, 1])   # 'To' and 'from' are defined in the sense of
   branchmat[:from] = Int.(branchdat[:, 2]) # towards the substation
   branchmat[:r] = branchdat[:, 3]
   branchmat[:x] = branchdat[:, 4]
   branchmat[:b] = branchdat[:, 5]
   # branchmat[:s] = branchdat[:, 6]
   branchmat[:s] = fill(2.0, nrow(branchmat))

   # Bus Data
   busdat = mpc["bus"]
   busmat = DataFrame()
   busmat[:Node] = Int.(busdat[:, 1])
   busmat[:Pd] = busdat[:, 3]
   busmat[:Qd] = busdat[:, 4]
   busmat[:Vmax] = busdat[:, 12]
   busmat[:Vmin] =busdat[:, 13]
   busmat[:B] =busdat[:, 6]
   rename!(branchmat, :from, :Node)
   busmat = join(busmat, branchmat[:, [:Node, :r, :x]], on=:Node, kind=:left)
   rename!(branchmat, :Node, :from)
   rename!(busmat, :x, :X)
   rename!(busmat, :r, :R)
   # busmat[isna.(busmat[:R]),:R] = 0.0
   # busmat[isna.(busmat[:X]),:X] = 0.0
   sort!(busmat, cols=:Node)

   # Generator Data
   gendat = mpc["gen"]
   genmat = DataFrame()
   genmat[:Node] = Int.(gendat[:, 1])
   genmat[:Pmax] = gendat[:, 9]
   genmat[:Qmax] = gendat[:, 4]
   costdat = mpc["gencost"]
   pw = []
   sort!(genmat, cols=:Node)
   genmat[:Cost] = zeros(nrow(genmat))
   for i in 1:nrow(genmat)
       pw_cost_funct = costdat[i, 5:end]
       filter!(x -> x!=0, pw_cost_funct)
       genmat[i, :Cost] = middle(pw_cost_funct)
   end

   # Fill gendata for every bus
   genmat = join(busmat[:, [:Node]], genmat, on=:Node, kind=:left)
   # genmat[isna.(genmat[:Pmax]),:Pmax] = 0.0
   # genmat[isna.(genmat[:Qmax]),:Qmax] = 0.0
   # genmat[isna.(genmat[:Cost]),:Cost] = 0.0

   # to p.u.
   # mpc already in PU!!
   # Vbase = mpc["bus"][1, 10] * 1e3
   # Sbase = mpc["baseMVA"] * 1e6

   # branchmat[:, :r] = branchmat[:, :r] ./ (Vbase^2 / Sbase)
   # branchmat[:, :x] = branchmat[:, :x] ./ (Vbase^2 / Sbase)
   # busmat[:, :R] = busmat[:, :R] ./ (Vbase^2 / Sbase)
   # busmat[:, :X] = busmat[:, :X] ./ (Vbase^2 / Sbase)
   # busmat[:, :Pd] = busmat[:, :Pd] ./ 1e3
   # busmat[:, :Qd] = busmat[:, :Qd] ./ 1e3

   return branchmat, busmat, genmat
end

function load_data(;remove_dg = false, folder = "", mpcasename = "", Sbase = 1, Vbase = 1)

# mpcasename = "case33bw"
# remove_dg = true
# folder = ""

# READ RAW
   if mpcasename != ""
      mpc = loadcase(mpcasename)
      branchmat, busmat, genmat = mpc_to_mat(mpc)
   else
      busmat = readtable("$folder/Node.csv")
      genmat = readtable("$folder/Generator.csv")
      branchmat = readtable("$folder/Line.csv")
   end

Zbase = (Vbase^2)/Sbase

busmat[:Pd] = busmat[:Pd] ./ Sbase
busmat[:Qd] = busmat[:Qd] ./ Sbase
busmat[:r]  = busmat[:r]  ./ Zbase
busmat[:x]  = busmat[:x]  ./ Zbase

branchmat[:r] = branchmat[:r] ./ Zbase
branchmat[:x] = branchmat[:x] ./ Zbase
branchmat[:s] = branchmat[:s] ./ Sbase   


# TO TYPELISTS
   # Buses
   buses = Bus[]
   busIDmap = Dict()
   for i in 1:nrow(busmat)
   	   nodeID = i
         busIDmap[busmat[i, :Node]] = i

         if i == 1
            root = busIDmap[busmat[i, :Node]]
         else
            root = 0
         end

         Pd = busmat[i, :Pd]
         Qd = busmat[i, :Qd]
         R = busmat[i, :r]
         X = busmat[i, :x]
         Vmax = busmat[i, :Vmax]
         Vmin = busmat[i, :Vmin]
         B = busmat[i, :b]

         b = Bus(nodeID, root, Pd, Qd, B, R, X, Vmax, Vmin)
         push!(buses, b)
   end
   numbuses = length(buses)

   ## Generators
   generatorlist = Int[]
   generators = Generator[]
   for i in 1:nrow(genmat)
      genID = i
      busidx = busIDmap[genmat[i, :Node]]
      Pgmax = genmat[i, :Pmax]
      Pgmin = 0.0
      Qgmax = genmat[i, :Qmax]
      Qgmin = 0.0
      g = Generator(genID, busidx, Pgmax, Pgmin, Qgmax, Qgmin)
      push!(generators, g)
      setg(buses[busidx], i)
   end
   for g in 1:length(generators)
      generators[g].cost = genmat[g, :Cost]
   end

   #Remove DGs
   if remove_dg & (length(generators) > 1)
      for g in 2:length(generators)
         generators[g].Pgmax = 0
         generators[g].Qgmax = 0
      end
   end

   # Branches
   lines = Line[]
   for i in 1:nrow(branchmat)
      fbus = busIDmap[branchmat[i, :from]]
      tbus = busIDmap[branchmat[i, :to]]
      abus = busIDmap[branchmat[i, :to]]
      x = branchmat[i, :x]
      r = branchmat[i, :r]
      u = branchmat[i, :s]

      push!(buses[tbus].children, fbus)#children
      push!(buses[fbus].ancestor, abus)#ancestor

      l = Line(i, tbus, fbus, r, x, u)
      push!(lines,l)
   end

   return buses, generators, lines
end
