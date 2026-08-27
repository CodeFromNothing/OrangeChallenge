### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 31854e29-621f-4a87-80aa-61043915041b
using JSON3

# ╔═╡ da17332c-082e-4f30-bfcc-77fa4ed970c8
using HiGHS

# ╔═╡ d3242b57-8215-44ff-bc53-ee3659fee836
using JuMP

# ╔═╡ 5ec53aae-42b4-4d18-bd71-f4386fb3b776
using Graphs

# ╔═╡ db10464d-96c0-4d1f-9f8f-9fc377ea0f59
using DataFrames

# ╔═╡ 289f0d9d-c1ad-4eaa-97c9-1b0acf2fbc60
using CSV

# ╔═╡ 395b32b3-b141-4e3e-ba92-559819bce91f
# Load network topology from JSON into a graph (NetworkGraph struct)

# ╔═╡ 0631f47f-47f9-4ff3-8db1-71caa488492d
struct NetworkGraph{G<:AbstractGraph}
    graph::G
    node_data::Vector{NamedTuple{(:json_id, :name), Tuple{Int, String}}}
    json_to_vertex::Dict{Int, Int}
    edge_data::Dict{
        Tuple{Int, Int},
        NamedTuple{(:id, :metric, :capacity), Tuple{Int, Float64, Float64}}
    }
end

# ╔═╡ 0dd31498-70b1-4e50-a899-c089bed4dc25
function json_to_network(json_text::AbstractString)
    data = JSON3.read(json_text)

    hasproperty(data, :directed) || error("Missing JSON field: directed")
    hasproperty(data, :multigraph) || error("Missing JSON field: multigraph")
    hasproperty(data, :nodes) || error("Missing JSON field: nodes")
    hasproperty(data, :links) || error("Missing JSON field: links")

    Bool(data.multigraph) &&
        error("This notebook uses SimpleGraph/SimpleDiGraph and therefore does not support multigraph=true.")

    n = length(data.nodes)

    json_to_vertex = Dict{Int, Int}()
    node_data = Vector{NamedTuple{(:json_id, :name), Tuple{Int, String}}}(undef, n)

    for (v, node) in enumerate(data.nodes)
        id = Int(node.id)
        haskey(json_to_vertex, id) && error("Duplicate node id in JSON: $id")
        json_to_vertex[id] = v
        node_data[v] = (json_id = id, name = String(node.name))
    end

    g = Bool(data.directed) ? SimpleDiGraph(n) : SimpleGraph(n)

    edge_data = Dict{
        Tuple{Int, Int},
        NamedTuple{(:id, :metric, :capacity), Tuple{Int, Float64, Float64}}
    }()

    for link in data.links
        from_id = Int(link.from)
        to_id   = Int(link.to)

        haskey(json_to_vertex, from_id) || error("Unknown node id in link: $from_id")
        haskey(json_to_vertex, to_id) || error("Unknown node id in link: $to_id")

        u = json_to_vertex[from_id]
        v = json_to_vertex[to_id]

        key = Bool(data.directed) ? (u, v) : minmax(u, v)

        haskey(edge_data, key) &&
            error("Multiple links between the same endpoints are not supported when multigraph=false.")

        added = add_edge!(g, u, v)
        added || error("Could not add edge ($from_id,$to_id). Check the JSON for duplicate links.")

        edge_data[key] = (id = Int(link.id), metric = Float64(link.metric), capacity = Float64(link.capacity))
    end

    return NetworkGraph(g, node_data, json_to_vertex, edge_data)
end

# ╔═╡ d283a9eb-f143-4a63-8b6a-95a6d5979050
function load_network_json(filename::AbstractString)
    return json_to_network(read(filename, String))
end

# ╔═╡ 09498a00-ce7b-410b-89bb-738a19c18a8b
json_id(net::NetworkGraph, v::Integer) = net.node_data[v].json_id

# ╔═╡ 36ba962c-7140-488b-9617-e7a4fe0476f9
node_name(net::NetworkGraph, v::Integer) = net.node_data[v].name

# ╔═╡ 225c4004-a90e-4cae-a834-75e93564d2a1
vertex_from_json_id(net::NetworkGraph, id::Integer) = net.json_to_vertex[Int(id)]

# ╔═╡ 63d4ffa9-048f-423f-ac3c-a2e8893b788f
function edge_attributes(net::NetworkGraph, u::Integer, v::Integer)
    key = is_directed(net.graph) ? (Int(u), Int(v)) : minmax(Int(u), Int(v))
    return net.edge_data[key]
end

# ╔═╡ 46809a87-4979-4dcd-ac4d-c01af752adfe
function metric_matrix(net::NetworkGraph)
    g = net.graph
    n = nv(g)
    D = fill(Inf, n, n)
    for v in vertices(g)
        D[v, v] = 0.0
    end
    for e in edges(g)
        u, v = src(e), dst(e)
        w = edge_attributes(net, u, v).metric
        w < 0 && error("Dijkstra's algorithm requires non-negative metric values.")
        D[u, v] = w
        if !is_directed(g)
            D[v, u] = w
        end
    end
    return D
end

# ╔═╡ 4c3c53f6-ab6a-4712-9235-f377b55c0224
network = load_network_json("setA-01-net.json")

# ╔═╡ 2ef5c08c-75cd-4fab-8e7c-392506b17886
nv(network.graph)

# ╔═╡ 21da7235-b934-424b-a588-70bf7c00155d
ne(network.graph)

# ╔═╡ 68d54122-0d7a-4e2b-94b9-f203e7015bc2
# Run Dijkstra from every vertex to get all-pairs distances (d) and shortest-path counts (sigma)

# ╔═╡ 146f861f-3813-4ff7-b20c-c8eeac9e9611
function all_pairs_shortest_data(g, distmx)
	n = nv(g)
# d[i,j] = shortest-path distance from i to j
	d = fill(Inf, n, n)
# sigma[i,j] = number of shortest paths from i to j
	sigma = zeros(Float64, n, n)
	for i in vertices(g)
		state = dijkstra_shortest_paths(
		g,
		i,
		distmx;
		allpaths = true
	)
	d[i, :] .= state.dists
	sigma[i, :] .= state.pathcounts
	end
	return d, sigma
end

# ╔═╡ 939ea361-6fab-4d3c-8f89-f279c93d31e7
# For a fixed segment (i,j), computes r(i,j,a,0) for every arc a — the fraction
# of i->j shortest-path traffic that flows through each arc

# ╔═╡ 29b9c6c6-53b5-4820-9314-7acfef8a2d43
function r_coefficients(net::NetworkGraph, d, sigma, i::Integer, j::Integer)
    n = nv(net.graph)
    r = zeros(Float64, n, n)   # r[u,v] = r(i,j,(u,v),0)

    sigma_ij = sigma[i, j]
    if sigma_ij == 0
        return r   # no path from i to j at all -> everything stays 0
    end

    for e in edges(net.graph)
        u, v = src(e), dst(e)
        ω = edge_attributes(net, u, v).metric
        ω == Inf && continue

        if d[i, u] + ω + d[v, j] == d[i, j]
            r[u, v] = (sigma[i, u] * sigma[v, j]) / sigma_ij
        end
    end

    return r
end

# ╔═╡ 99e38a49-3eed-40ca-aaa4-e3659ee64828
# Uses distance/path-count matrices instead of enumerating shortest paths explicitly (avoids combinatorial blowup)

# ╔═╡ b1733e2f-7b0b-4621-8bf7-58902304b111
begin
    # Reads demand list + traffic volumes ν(d,t) from the traffic matrix JSON
    struct Demand
        s::Int          # source (Julia vertex, after remapping)
        t::Int          # target (Julia vertex, after remapping)
        volumes::Vector{Float64}   # ν(d, time_step) for each time step
    end
    
    function load_traffic_matrix(filename::AbstractString, net::NetworkGraph)
        data = JSON3.read(read(filename, String))
        demands = Demand[]
        for dem in data.demands
            s = vertex_from_json_id(net, Int(dem.s))
            t = vertex_from_json_id(net, Int(dem.t))
            volumes = Float64.(dem.v)
            push!(demands, Demand(s, t, volumes))
        end
        return demands
    end
end

# ╔═╡ e9f57963-0df1-4a4e-bcf0-87250e7abbe6
# Reads maxSeg, budget κ(t), and intervention scenario q(t) from the scenario JSON
function load_scenario(filename::AbstractString)
    data = JSON3.read(read(filename, String))
    max_seg = Int(data.max_segments)

    budget = Dict{Int,Int}()
    for b in data.budget
        budget[Int(b.t)] = Int(b.value)
    end

    interventions = Dict{Int,Vector{Int}}()
    for iv in data.interventions
        interventions[Int(iv.t)] = Int.(iv.links)
    end

    return max_seg, budget, interventions
end

# ╔═╡ b3433f28-0ad2-4025-8931-da86bd1c68b2
# Step 3 — Data: assembles r(i,j,a,0), c(a), ν(d,0), maxSeg for the model
function build_data(net_file, tm_file, scenario_file)
    net = load_network_json(net_file)
    demands = load_traffic_matrix(tm_file, net)
    max_seg, _, _ = load_scenario(scenario_file)   # budget/interventions unused at t=0

    d, sigma = all_pairs_shortest_data(net.graph, metric_matrix(net))

    n = nv(net.graph)
    r = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    for i in 1:n, j in 1:n
        r[(i,j)] = r_coefficients(net, d, sigma, i, j)
    end

    return net, demands, max_seg, r
end

# ╔═╡ 0145572b-8d63-4de6-9852-a5f7e7537b2c
instances = ["setA-01", "setA-02", "setA-03", "setA-04", "setA-05",
             "setA-06", "setA-07", "setA-08", "setA-09", "setA-10", "setA-11", "setA-12", "setA-13", "setA-14", "setA-15", "setA-16", "setA-17", "setA-18", "setA-19", "setA-20"]

# ╔═╡ dc2e6e08-8028-4002-a3d6-74ea8f3d2b9e
begin
    function run_instance(instance_prefix::AbstractString)
        net_file      = instance_prefix * "-net.json"
        tm_file       = instance_prefix * "-tm.json"
        scenario_file = instance_prefix * "-scenario.json"
    
        net, demands, max_seg, r = build_data(net_file, tm_file, scenario_file)
        model, x, λmax = build_period0_model(net, demands, max_seg, r)
    
        t_start = time()
        optimize!(model)
        cpu_time = time() - t_start
    
        status = termination_status(model)
    
        if status == MOI.OPTIMAL
            obj = objective_value(model)
            gap = 0.0
        elseif status == MOI.INFEASIBLE
            obj = missing
            gap = missing
        elseif has_values(model)
            obj = objective_value(model)
            lb  = objective_bound(model)
            gap = 1.0 - lb / obj
        else
            obj, gap = missing, missing
        end
    
        return (
            instance    = instance_prefix,
            n_vertices  = nv(net.graph),
            n_links     = ne(net.graph),
            n_demands   = length(demands),
            status      = string(status),
            objective   = obj,
            cpu_time    = cpu_time,
            gap         = gap,
        )
    end

    function run_instance_safe(inst)
        try
            return run_instance(inst)
        catch err
            return (instance=inst, n_vertices=missing, n_links=missing, n_demands=missing,
                    status="ERROR: $err", objective=missing, cpu_time=missing, gap=missing)
        end
    end

    results = [run_instance_safe(inst) for inst in instances]
    df = DataFrame(results)
end

# ╔═╡ 61186540-cf7c-4cc8-80ae-0b1b9c41c293
CSV.write("results_period0.csv", df)

# ╔═╡ 3f9ae4cd-7504-4696-931a-f144fbebb95e
function metric_matrix_with_interventions(net::NetworkGraph, down_link_ids::Vector{Int})
    D = metric_matrix(net)
    down_set = Set(down_link_ids)
    for ((u, v), attrs) in net.edge_data
        if attrs.id in down_set
            D[u, v] = Inf
        end
    end
    return D
end

# ╔═╡ a0beb51a-866b-4d77-a544-59e31ee20cf0
function build_two_period_data(net_file, tm_file, scenario_file)
    net = load_network_json(net_file)
    demands = load_traffic_matrix(tm_file, net)
    max_seg, budget, interventions = load_scenario(scenario_file)
    n = nv(net.graph)

    D0 = metric_matrix(net)
    d0, sigma0 = all_pairs_shortest_data(net.graph, D0)

    down_ids = get(interventions, 1, Int[])
    D1 = metric_matrix_with_interventions(net, down_ids)
    d1, sigma1 = all_pairs_shortest_data(net.graph, D1)

    r0 = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    r1 = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    for i in 1:n, j in 1:n
        r0[(i,j)] = r_coefficients(net, d0, sigma0, D0, i, j)
        r1[(i,j)] = r_coefficients(net, d1, sigma1, D1, i, j)
    end

    return net, demands, max_seg, budget, Set(down_ids), r0, r1
end

# ╔═╡ 7b93da35-51eb-4c89-9bf7-8dd5943e3de6
function build_two_period_model(net::NetworkGraph, demands::Vector{Demand}, max_seg::Int,
                                 budget::Dict{Int,Int}, down_ids::Set{Int},
                                 r0::Dict{Tuple{Int,Int}, Matrix{Float64}},
                                 r1::Dict{Tuple{Int,Int}, Matrix{Float64}};
                                 time_limit::Float64 = 1800.0)
    n  = nv(net.graph)
    Ds = 1:length(demands)
    T  = 0:1

    model = Model(HiGHS.Optimizer)
    set_time_limit_sec(model, time_limit)

    @variable(model, x[d in Ds, t in T, i in 1:n, j in 1:n; i != j], Bin)
    @variable(model, λmax >= 0)
    @variable(model, z[d in Ds, i in 1:n, j in 1:n; i != j] >= 0)

    # (1) flow conservation, each period
    for t in T, d in Ds
        s, tt = demands[d].s, demands[d].t
        for i in 1:n
            rhs = i == s ? 1 : (i == tt ? -1 : 0)
            @constraint(model,
                sum(x[d,t,i,j] for j in 1:n if j != i) -
                sum(x[d,t,j,i] for j in 1:n if j != i) == rhs)
        end
    end

    # (2) segment cap, each period
    for t in T, d in Ds
        @constraint(model, sum(x[d,t,i,j] for i in 1:n, j in 1:n if i != j) <= max_seg)
    end

    # (3) load: period 0 all arcs, period 1 skips down arcs
    for e in edges(net.graph)
        u, v = src(e), dst(e)
        link_id = edge_attributes(net, u, v).id
        cap = edge_attributes(net, u, v).capacity

        @constraint(model,
            sum(r0[(i,j)][u,v] * demands[d].volumes[1] * x[d,0,i,j]
                for d in Ds, i in 1:n, j in 1:n if i != j) <= λmax * cap)

        if !(link_id in down_ids)
            @constraint(model,
                sum(r1[(i,j)][u,v] * demands[d].volumes[2] * x[d,1,i,j]
                    for d in Ds, i in 1:n, j in 1:n if i != j) <= λmax * cap)
        end
    end

    # (4) budget constraint between period 0 and period 1
    κ1 = get(budget, 1, typemax(Int))
    for d in Ds, i in 1:n, j in 1:n
        i == j && continue
        @constraint(model, z[d,i,j] >= x[d,1,i,j] - x[d,0,i,j])
        @constraint(model, z[d,i,j] >= x[d,0,i,j] - x[d,1,i,j])
    end
    @constraint(model, sum(z[d,i,j] for d in Ds, i in 1:n, j in 1:n if i != j) <= κ1)

    @objective(model, Min, λmax)

    return model, x, λmax
end

# ╔═╡ 2386dda0-2be8-41d9-abcf-bb6303c8b100
begin
    function run_instance_two_period(instance_prefix::AbstractString)
        net_file      = instance_prefix * "-net.json"
        tm_file       = instance_prefix * "-tm.json"
        scenario_file = instance_prefix * "-scenario.json"

        net, demands, max_seg, budget, down_ids, r0, r1 = build_two_period_data(net_file, tm_file, scenario_file)
        model, x, λmax = build_two_period_model(net, demands, max_seg, budget, down_ids, r0, r1; time_limit=1800.0)

        t_start = time()
        optimize!(model)
        cpu_time = time() - t_start

        status = termination_status(model)

        if status == MOI.OPTIMAL
            obj, gap = objective_value(model), 0.0
        elseif status == MOI.INFEASIBLE
            obj, gap = missing, missing
        elseif has_values(model)
            obj = objective_value(model)
            lb  = objective_bound(model)
            gap = 1.0 - lb / obj
        else
            obj, gap = missing, missing
        end

        return (
            instance   = instance_prefix,
            n_vertices = nv(net.graph),
            n_links    = ne(net.graph),
            n_demands  = length(demands),
            status     = string(status),
            objective  = obj,
            cpu_time   = cpu_time,
            gap        = gap,
        )
    end

    function run_instance_two_period_safe(inst)
        try
            return run_instance_two_period(inst)
        catch err
            return (instance=inst, n_vertices=missing, n_links=missing, n_demands=missing,
                    status="ERROR: $err", objective=missing, cpu_time=missing, gap=missing)
        end
    end

    results8 = [run_instance_two_period_safe(inst) for inst in instances]
    df8 = DataFrame(results8)
end

# ╔═╡ fb75dcce-90b0-4b0f-be24-c070dd8bdf65
CSV.write("results_period0_1.csv", df8)

# ╔═╡ 245111c8-c16d-4584-b7f1-f9b41594baec
begin
	net, demands, max_seg, r = build_data("setA-01-net.json", "setA-01-tm.json", "setA-01-scenario.json")
	model, x, λmax = build_period0_model(net, demands, max_seg, r)
	optimize!(model)
	objective_value(model)   # the minimized max link utilization at t=0
end

# ╔═╡ c5c1b4c9-a317-46f3-b41b-af3c60f5c181
begin
    function build_period0_model(net::NetworkGraph, demands::Vector{Demand}, max_seg::Int,
                                  r::Dict{Tuple{Int,Int}, Matrix{Float64}})
        n  = nv(net.graph)
        Ds = 1:length(demands)
    
        model = Model(HiGHS.Optimizer)
        set_time_limit_sec(model, 900.0)
    
        # x[d,i,j] only makes sense for i != j
        @variable(model, x[d in Ds, i in 1:n, j in 1:n; i != j], Bin)
        @variable(model, λmax >= 0)
    
        # (1) flow conservation
        for d in Ds
            s, t = demands[d].s, demands[d].t
            for i in 1:n
                rhs = i == s ? 1 : (i == t ? -1 : 0)
                @constraint(model,
                    sum(x[d,i,j] for j in 1:n if j != i) -
                    sum(x[d,j,i] for j in 1:n if j != i) == rhs)
            end
        end
    
        # (2) segment cap
        for d in Ds
            @constraint(model, sum(x[d,i,j] for i in 1:n, j in 1:n if i != j) <= max_seg)
        end
    
        # (3) load, only over real arcs of the graph (q(0) = ∅, so all arcs are available)
        for e in edges(net.graph)
            u, v = src(e), dst(e)
            cap  = edge_attributes(net, u, v).capacity
            @constraint(model,
                sum(r[(i,j)][u,v] * demands[d].volumes[1] * x[d,i,j]
                    for d in Ds, i in 1:n, j in 1:n if i != j) <= λmax * cap)
        end
    
        # (5), top level of the lex-min: minimize the worst-case load
        @objective(model, Min, λmax)
    
        return model, x, λmax
    end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
HiGHS = "87dc4568-4c63-4d18-b0c0-bb2238e4078b"
JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"

[compat]
CSV = "~0.10.16"
DataFrames = "~1.8.2"
Graphs = "~1.14.0"
HiGHS = "~1.24.1"
JSON3 = "~1.14.3"
JuMP = "~1.31.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.7"
manifest_format = "2.0"
project_hash = "60d490b4f6b0bb108b4560fcb6ccc8b568e43245"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "8d8e0b0f350b8e1c91420b5e64e5de774c2f0f4d"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.16"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "TranscodingStreams"]
git-tree-sha1 = "84990fa864b7f2b4901901ca12736e45ee79068c"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.8.5"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "970758a3d591a2a5c2a907c53f2e2f8c1b1d3537"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.9"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.1+2"

[[deps.Crayons]]
git-tree-sha1 = "54b76cbb40d9a0f5368c880725b2f141da77c94f"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.2.0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "5fab31e2e01e70ad66e3e24c968c264d1cf166d6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.2"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "79a2aca180a85c690c58a020d47b426954b590f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.16.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "1b86cca764a61dcac4fef4c5e16e378e5ed6953c"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.4.5"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Inflate", "LinearAlgebra", "Random", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7eb45fe833a5b7c51cf6d89c5a841d5967e44be3"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.14.0"

    [deps.Graphs.extensions]
    GraphsSharedArraysExt = "SharedArrays"

    [deps.Graphs.weakdeps]
    Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"
    SharedArrays = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.HiGHS]]
deps = ["HiGHS_jll", "LinearAlgebra", "MathOptIIS", "MathOptInterface", "OpenBLAS32_jll", "PrecompileTools", "SparseArrays"]
git-tree-sha1 = "01a5241985559c08a5baadbcebd6d87daaf84a84"
uuid = "87dc4568-4c63-4d18-b0c0-bb2238e4078b"
version = "1.24.1"

[[deps.HiGHS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Zlib_jll", "libblastrampoline_jll"]
git-tree-sha1 = "2d9747b79d17c4320fe48048a3a768fe6d6d82de"
uuid = "8fd58aa0-07eb-5a78-9b36-339c94fd15ea"
version = "1.15.1+1"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "8f3d257792a522b4601c24a577954b0a8cd7334d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.5"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "c7345ab1a7ca4dc8a02c9f6510da0d9857bbe513"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.7.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "411eccfe8aba0814ffa0fdf4860913ed09c34975"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.3"

    [deps.JSON3.extensions]
    JSON3ArrowExt = ["ArrowTypes"]

    [deps.JSON3.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JuMP]]
deps = ["LinearAlgebra", "MacroTools", "MathOptInterface", "MutableArithmetics", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays"]
git-tree-sha1 = "614b22ff014355192982b1f9a12c61298ce6a908"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "1.31.1"

    [deps.JuMP.extensions]
    JuMPDimensionalDataExt = "DimensionalData"

    [deps.JuMP.weakdeps]
    DimensionalData = "0703355e-b756-11e9-17c0-8b28908087d0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "bba2d9aa057d8f126415de240573e86a8f39d2a1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.1"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathOptIIS]]
deps = ["MathOptInterface"]
git-tree-sha1 = "3b3d69130d8ab8c39d5fa4d30e20a8e6428c9d37"
uuid = "8c4f8055-bd93-4160-a86b-a0c04941dbff"
version = "0.2.0"

[[deps.MathOptInterface]]
deps = ["CodecBzip2", "CodecZlib", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays", "SpecialFunctions", "Test"]
git-tree-sha1 = "f1ccd9ffcb8577e207deb9aaebeb3f961de70380"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.52.0"

    [deps.MathOptInterface.extensions]
    MathOptInterfaceBenchmarkToolsExt = "BenchmarkTools"
    MathOptInterfaceCliqueTreesExt = "CliqueTrees"

    [deps.MathOptInterface.weakdeps]
    BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "dc5b2c4c111c46bc79ac4405eeb563523b39c004"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.8.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "libblastrampoline_jll"]
git-tree-sha1 = "30870d0f2dc0b2dba76b10df1c58c7f018413e56"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.34+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "3de8f5e6e90ebfa8d6d1f86997d6cdcd6a912ff3"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.7"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "4ac881f5432bd93463a41767a814a45245be22b6"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.4.6"

    [deps.PrettyTables.extensions]
    PrettyTablesExcelExt = "XLSX"
    PrettyTablesTypstryExt = "Typstry"

    [deps.PrettyTables.weakdeps]
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"
    XLSX = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "084c47c7c5ce5cfecefa0a98dff69eb3646b5a80"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.10"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "429071b23f4c9a13fb6582f807cc2ef454082408"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.9.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "fac51faf3bb96e8bc0bf6f9f39ca4955652776bb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.19"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "8a90c1d77c3277a5d43b83927b3cbe2c70a37484"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.7"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "2d0fc55c61321ba245c47be599570d11bac50303"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.5"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "0f38a06c83f0007bbab3cf911262841c9a0f07e0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.13.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "0716e01c3b40413de5dedbc9c5c69f27cddfddfc"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.3"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"
"""

# ╔═╡ Cell order:
# ╠═395b32b3-b141-4e3e-ba92-559819bce91f
# ╠═31854e29-621f-4a87-80aa-61043915041b
# ╠═da17332c-082e-4f30-bfcc-77fa4ed970c8
# ╠═d3242b57-8215-44ff-bc53-ee3659fee836
# ╠═5ec53aae-42b4-4d18-bd71-f4386fb3b776
# ╠═db10464d-96c0-4d1f-9f8f-9fc377ea0f59
# ╠═289f0d9d-c1ad-4eaa-97c9-1b0acf2fbc60
# ╠═0631f47f-47f9-4ff3-8db1-71caa488492d
# ╠═0dd31498-70b1-4e50-a899-c089bed4dc25
# ╠═d283a9eb-f143-4a63-8b6a-95a6d5979050
# ╠═09498a00-ce7b-410b-89bb-738a19c18a8b
# ╠═36ba962c-7140-488b-9617-e7a4fe0476f9
# ╠═225c4004-a90e-4cae-a834-75e93564d2a1
# ╠═63d4ffa9-048f-423f-ac3c-a2e8893b788f
# ╠═46809a87-4979-4dcd-ac4d-c01af752adfe
# ╠═4c3c53f6-ab6a-4712-9235-f377b55c0224
# ╠═2ef5c08c-75cd-4fab-8e7c-392506b17886
# ╠═21da7235-b934-424b-a588-70bf7c00155d
# ╠═68d54122-0d7a-4e2b-94b9-f203e7015bc2
# ╠═146f861f-3813-4ff7-b20c-c8eeac9e9611
# ╠═939ea361-6fab-4d3c-8f89-f279c93d31e7
# ╠═29b9c6c6-53b5-4820-9314-7acfef8a2d43
# ╠═99e38a49-3eed-40ca-aaa4-e3659ee64828
# ╠═b1733e2f-7b0b-4621-8bf7-58902304b111
# ╠═e9f57963-0df1-4a4e-bcf0-87250e7abbe6
# ╠═b3433f28-0ad2-4025-8931-da86bd1c68b2
# ╠═c5c1b4c9-a317-46f3-b41b-af3c60f5c181
# ╠═245111c8-c16d-4584-b7f1-f9b41594baec
# ╠═0145572b-8d63-4de6-9852-a5f7e7537b2c
# ╠═dc2e6e08-8028-4002-a3d6-74ea8f3d2b9e
# ╠═61186540-cf7c-4cc8-80ae-0b1b9c41c293
# ╠═3f9ae4cd-7504-4696-931a-f144fbebb95e
# ╠═a0beb51a-866b-4d77-a544-59e31ee20cf0
# ╠═7b93da35-51eb-4c89-9bf7-8dd5943e3de6
# ╠═2386dda0-2be8-41d9-abcf-bb6303c8b100
# ╠═fb75dcce-90b0-4b0f-be24-c070dd8bdf65
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
