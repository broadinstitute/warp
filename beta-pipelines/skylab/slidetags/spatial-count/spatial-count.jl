using CSV
using HDF5
using FASTX
using CodecZlib
using IterTools: product
using StatsBase: countmap, sample
using DataFrames
using StringViews
using LinearAlgebra: dot
using Combinatorics: combinations

# fastqpath is the path to a directory containing all fastq files to process
# puckpath is the path to a directory containing all puck files to reference
# Running this script will process all files in both directories and produce one output
# - so be sure that these directories only contain the files for a specific run

# Read the command-line arguments
if length(ARGS) != 2
    error("Usage: julia spatial-count.jl fastq_path puck_path")
end
fastqpath = ARGS[1]
println("FASTQ path: "*fastqpath)
@assert isdir(fastqpath) "FASTQ path not found"
@assert !isempty(readdir(fastqpath)) "FASTQ path is empty"
puckpath = ARGS[2]
println("Puck path: "*puckpath)
@assert isdir(puckpath) "Puck path not found"
@assert !isempty(readdir(puckpath)) "Puck path is empty"

# Recognized bead types:
# JJJJJJJJ  TCTTCAGCGTTCCCGAGA JJJJJJJ  NNNNNNNVV (V10)
# JJJJJJJJJ TCTTCAGCGTTCCCGAGA JJJJJJJJ NNNNNNNNN (V17)
# JJJJJJJJJJJJJJJ   CTGTTTCCTG NNNNNNNNN          (V15)
# JJJJJJJJJJJJJJJJJ CTGTTTCCTG NNNNNNNNN          (V16)

##### Load the puck data #######################################################

# Load the pucks
puck_paths = filter(x -> endswith(x, ".csv"), readdir(puckpath, join=true)) ; println("Pucks: ", basename.(puck_paths))
puckdfs = [rename(CSV.read(puck, DataFrame, header=false, types=[String, Float64, Float64]), [:sb,:x,:y]) for puck in puck_paths]

# Get the spatial barcode length
function get_sb_length(puckdfs)
    sb_sizes = union([Set(length(s) for s in puckdf.sb) for puckdf in puckdfs]...)
    @assert length(sb_sizes) == 1 "ERROR: Not all puck spatial barcodes have the same length"
    return(first(sb_sizes))
end
const sb_len = get_sb_length(puckdfs) ; println("Spatial barcode length: $sb_len")
@assert sb_len in [14, 15, 17] "ERROR: Unrecognized spatial barcode length"

# Validate the pucks
for (puck, puckdf) in zip(puck_paths, puckdfs)
    println("Loaded $(basename(puck)): $(nrow(puckdf)) spatial barcodes found")
    @assert all(!ismissing, puckdf[!, :sb]) "ERROR: Puck column <sb> contains missing values"
    @assert all(!ismissing, puckdf[!, :x]) "ERROR: Puck column <x> contains missing values"
    @assert all(!ismissing, puckdf[!, :y]) "ERROR: Puck column <y> contains missing values"
    @assert all(length(s) == sb_len for s in puckdf.sb) "ERROR: Some spatial barcodes in $puck do not have $(sb_len)bp"
    @assert length(puckdf.sb) == length(Set(puckdf.sb)) "ERROR: Some spatial barcodes in $puck are repeated"
end

# Check for low quality spatial barcodes
const num_lowQbeads = [sum(count(c -> c == 'N', sb) > 1 for sb in puckdf.sb) for puckdf in puckdfs]
println("Removed $(sum(num_lowQbeads)) bead barcode(s) for having 2+ N bases")    

# Concatenate the puck dataframes
for (i, puckdf) in enumerate(puckdfs)
    puckdf[!, :puck_index] = fill(UInt8(i), nrow(puckdf))
end
const puckdf = vcat(puckdfs...)
empty!(puckdfs)

# Create the sb_whitelist
const sb_whitelist = [str for str in sort(collect(Set(puckdf.sb))) if count(c -> c == 'N', str) < 2]
println("Total number of unique spatial barcodes: $(length(sb_whitelist))")
@assert length(sb_whitelist) < 2^32-1 "Must change type of sb_i from U32 to U64"

println("") ; flush(stdout) ; GC.gc()

##### Load the FASTQ data ######################################################

# Load the FASTQ paths
fastq_paths = readdir(fastqpath, join=true)
fastq_paths = filter(fastq -> endswith(fastq, ".fastq.gz"), fastq_paths)
@assert length(fastq_paths) > 1 "ERROR: No FASTQ pairs found"
const R1s = filter(s -> occursin("_R1_", s), fastq_paths) ; println("R1s: ", basename.(R1s))
const R2s = filter(s -> occursin("_R2_", s), fastq_paths) ; println("R2s: ", basename.(R2s))
@assert length(R1s) > 0 && length(R2s) > 0 "ERROR: No FASTQ pairs found"
@assert length(R1s) == length(R2s) "ERROR: R1s and R2s are not all paired"
@assert [replace(R1, "_R1_"=>"", count=1) for R1 in R1s] == [replace(R2, "_R2_"=>"", count=1) for R2 in R2s]
println("$(length(R1s)) pair(s) of FASTQs found\n")

# Read structure methods
const SeqView = StringView{SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true}}
@inline function get_V10(seq::SeqView)
    @inbounds sb_1 = seq[1:8]
    @inbounds up = seq[9:26]
    @inbounds sb_2 = seq[27:33]
    return sb_1, sb_2, up
end
@inline function get_V17(seq::SeqView)
    @inbounds sb_1 = seq[1:9]
    @inbounds up = seq[10:27]
    @inbounds sb_2 = seq[28:35]
    return sb_1, sb_2, up
end
@inline function get_V15(seq::SeqView)
    @inbounds sb_1 = seq[1:8]
    @inbounds sb_2 = seq[9:15]
    @inbounds up = seq[16:25]
    return sb_1, sb_2, up
end
@inline function get_V16(seq::SeqView)
    @inbounds sb_1 = seq[1:9]
    @inbounds sb_2 = seq[10:17]
    @inbounds up = seq[18:27]
    return sb_1, sb_2, up
end
const UP1 = "TCTTCAGCGTTCCCGAGA"
const UP2 = "CTGTTTCCTG"

# UMI compressing (between 0x00000000 and 0x00ffffff for a 12bp UMI)
const px12 = [convert(UInt32, 4^i) for i in 0:11]
function UMItoindex(UMI::SeqView)::UInt32
    return dot(px12, (codeunits(UMI).>>1).&3)
end
# Convert compressed UMIs back into strings
const bases = ['A','C','T','G'] # MUST NOT change this order
function indextoUMI(i::UInt32)::String
    return String([bases[(i>>n)&3+1] for n in 0:2:22])
end

# Learn the read structure
function learn_bead(R)
    iter = R |> open |> GzipDecompressorStream |> FASTQ.Reader
    counts = Dict("V10"=>0, "V17"=>0, "V15"=>0, "V16"=>0)
    for (i, record) in enumerate(iter)
        i > 100000 ? break : nothing
        seq = FASTQ.sequence(record)
        length(seq) < 25 ? continue : nothing
        counts["V15"] += get_V15(seq)[3] == UP2
        length(seq) < 27 ? continue : nothing
        counts["V16"] += get_V16(seq)[3] == UP2
        length(seq) < 33 ? continue : nothing
        counts["V10"] += get_V10(seq)[3] == UP1
        length(seq) < 35 ? continue : nothing
        counts["V17"] += get_V17(seq)[3] == UP1
    end
    @assert sum(values(counts)) > 0
    println(counts)
    return findmax(counts)
end

R1_types = [learn_bead(R1) for R1 in R1s]
R2_types = [learn_bead(R2) for R2 in R2s]

switch_list = [r1[1] > r2[1] for (r1, r2) in zip(R1_types, R2_types)]
@assert all(switch_list) || !any(switch_list) "ERROR: the R1/R2 read assignment is not consistent"
const switch = all(switch_list)
println("Switch R1/R2: $switch")
if switch == true
    println("Switching R1 and R2")
    temp = R1s ; R1s = R2s ; R2s = temp
    temp = R1_types ; R1_types = R2_types ; R2_types = temp
end

bead_list = [r[2] for r in R2_types]
@assert length(Set(bead_list)) == 1 "ERROR: the bead type is not consistent"
const bead = first(Set(bead_list))
println("Bead type: $bead")

if bead == "V10"
    get_SB = get_V10
    const UP = UP1
    const R2_len = 33
    if sb_len == 14 # for backwards compatibility
        println("NOTE: Puck has 14-bp barcodes, but V10 beads have 15-bp barcodes - only using first 14 bases")
        const sb1_len = 8
        const sb2_len = 6
    else
        @assert sb_len == 15 "ERROR: Puck has $sb_len-bp barcodes, but V10 beads have 15-bp barcodes"
        const sb1_len = 8
        const sb2_len = 7
    end
elseif bead == "V17"
    get_SB = get_V17
    const UP = UP1
    const R2_len = 35
    @assert sb_len == 17 "ERROR: Puck has $sb_len-bp barcodes, but V17 beads have 17-bp barcodes"
    const sb1_len = 9
    const sb2_len = 8
elseif bead == "V15"
    get_SB = get_V15
    const UP = UP2
    const R2_len = 25
    @assert sb_len == 15 "ERROR: Puck has $sb_len-bp barcodes, but V15 beads have 15-bp barcodes"
    const sb1_len = 8
    const sb2_len = 7
elseif bead == "V16"
    get_SB = get_V16
    const UP = UP2
    const R2_len = 27
    @assert sb_len == 17 "ERROR: Puck has $sb_len-bp barcodes, but V16 beads have 17-bp barcodes"
    const sb1_len = 9
    const sb2_len = 8
else
    error("ERROR: Unhandled bead type ($bead)")
end

##### Helper methods ###########################################################

# Returns a set of all strings within a certain hamming distance of the input
function listHDneighbors(str, hd, charlist = ['A','C','G','T','N'])::Set{String}
    res = Set{String}()
    for inds in combinations(1:length(str), hd)
        chars = [str[i] for i in inds]
        pools = [setdiff(charlist, [char]) for char in chars]
        prods = product(pools...)
        for prod in prods
            s = str
            for (i, c) in zip(inds, prod)
                s = s[1:i-1]*string(c)*s[i+1:end]
            end
            push!(res,s)
        end
    end
    return(res)
end

# Returns a (string, index) tuple for all strings withing 1 HD of the input
# the index stores the 1-indexed position of the mismatched base
function listHD1neighbors(str, charlist = ['A','C','G','T','N'])::Vector{Tuple{String, UInt8}}
    res = Vector{Tuple{String, UInt8}}()
    for i in 1:length(str)
        for char in setdiff(charlist, str[i])
            s = String(str[1:i-1]*string(char)*str[i+1:end])
            push!(res, (s,i))
        end
    end
    return(res)
end

# Given a string with N, return a list of strings with all possible replacements
function expand_N(s, charlist = ['A','C','G','T'])::Vector{String}
    if !occursin('N', s)
        return [s]
    end
    combins = String[]
    for nucleotide in charlist
        new_str = replace(s, 'N' => nucleotide, count=1)
        append!(combins, expand_N(new_str))
    end
    return combins
end

# Create spatial barcode matching dictionary (sb1, sb2) -> (sb_i, pos)
# sb_i: ==0 is ambiguous, >0 is the unique index into sb_whitelist
# pos: the position of the fuzzy matched base, 0 for exact/ambig matches
function create_SBtoindex(sb_whitelist)
    SBtoindex = Dict{Tuple{String15, String15}, Tuple{UInt32, Int8}}()

    # Fuzzy matches
    for (i, sb) in enumerate(sb_whitelist)
        numN = count(c -> c == 'N', sb)
        if numN == 0
            # Add Hamming distance = 1 strings
            for (sb_f, ind) in listHD1neighbors(sb)
                sbtuple = (sb_f[1:sb1_len], sb_f[sb1_len+1:sb1_len+sb2_len])
                haskey(SBtoindex, sbtuple) ? SBtoindex[sbtuple] = (0, 0) : SBtoindex[sbtuple] = (i, ind)
            end
        elseif numN == 1
            ind = findfirst(isequal('N'), sb)
            sbtuple = (sb[1:sb1_len], sb[sb1_len+1:sb1_len+sb2_len])
            haskey(SBtoindex, sbtuple) ? SBtoindex[sbtuple] = (0, 0) : SBtoindex[sbtuple] = (i, ind)
            for sb_f in expand_N(sb)
                sbtuple = (sb_f[1:sb1_len], sb_f[sb1_len+1:sb1_len+sb2_len])
                haskey(SBtoindex, sbtuple) ? SBtoindex[sbtuple] = (0, 0) : SBtoindex[sbtuple] = (i, ind)
            end
        end
    end

    # Exact matches
    for (i, sb) in enumerate(sb_whitelist)
        if !occursin('N', sb)
            sbtuple = (sb[1:sb1_len], sb[sb1_len+1:sb1_len+sb2_len])
            SBtoindex[sbtuple] = (i, 0)
        end
    end
        
    return(SBtoindex)
end

print("Creating matching dictionaries... ") ; flush(stdout)

const SBtoindex = create_SBtoindex(sb_whitelist)

const UP_TOL = round(Int, length(UP)/6)
const UP_HD_whitelist = reduce(union, [listHDneighbors(UP, i) for i in 0:UP_TOL])
const UP_GG_whitelist = reduce(union, [listHDneighbors("G"^length(UP), i) for i in 0:UP_TOL])

const UMI_TOL = round(Int, 12/6)
const umi_homopolymer_whitelist = reduce(union, [listHDneighbors(c^12, i) for c in bases for i in 0:UMI_TOL])

println("done") ; flush(stdout) ; GC.gc()

##### Read the FASTQS ##########################################################

println("Reading FASTQs... ") ; flush(stdout)

function process_fastqs(R1s, R2s)
    # Create data structures
    mat = Dict{Tuple{UInt32, UInt32, UInt32}, UInt32}() # (cb_i, umi_i, sb_i) -> reads
    cb_dictionary = Dict{String31, UInt32}() # cb -> cb_i
    metadata = Dict("reads"=>0, "reads_filtered"=>0,
                    "R1_tooshort"=>0, "R2_tooshort"=>0,
                    "UMI"=>Dict("N"=>0, "homopolymer"=>0),
                    "UP"=>Dict("exact"=>0, "fuzzy"=>0, "GG"=>0, "none"=>0),
                    "SB"=>Dict("exact"=>0, "HD1"=>0, "HD1ambig"=>0, "none"=>0),
                    "SB_HD"=>Dict(i=>0 for i in collect(1:sb_len)))
    
    for fastqpair in zip(R1s, R2s)
        println(fastqpair) ; flush(stdout)
        it1 = fastqpair[1] |> open |> GzipDecompressorStream |> FASTQ.Reader;
        it2 = fastqpair[2] |> open |> GzipDecompressorStream |> FASTQ.Reader;
        for record in zip(it1, it2)
            metadata["reads"] += 1

            seq1 = FASTQ.sequence(record[1])
            seq2 = FASTQ.sequence(record[2])

            # Validate the sequence length
            skip = false
            if length(seq1) < 28
                metadata["R1_tooshort"] += 1
                skip = true
            end
            if length(seq2) < R2_len
                metadata["R2_tooshort"] += 1
                skip = true
            end
            if skip
                continue
            end

            # Load R1
            cb  = seq1[1:16]
            umi = seq1[17:28]
            # Filter UMI
            if occursin('N', umi)
                metadata["UMI"]["N"] += 1
                continue
            elseif in(umi, umi_homopolymer_whitelist)
                metadata["UMI"]["homopolymer"] += 1
                continue
            end

            # Load R2
            sb1, sb2, up = get_SB(seq2)
            # Filter UP
            if up == UP
                metadata["UP"]["exact"] += 1
            elseif in(up, UP_HD_whitelist)
                metadata["UP"]["fuzzy"] += 1
            elseif in(up, UP_GG_whitelist)
                metadata["UP"]["GG"] += 1
                continue
            else
                metadata["UP"]["none"] += 1
                continue
            end
            # Filter SB
            sb_i, ind = get(SBtoindex, (sb1, sb2), (-1, 0))
            if sb_i > 0 && ind == 0
                metadata["SB"]["exact"] += 1
            elseif sb_i > 0 && ind > 0
                metadata["SB"]["HD1"] += 1
                metadata["SB_HD"][ind] += 1
            elseif sb_i == 0
                metadata["SB"]["HD1ambig"] += 1
                continue
            else
                metadata["SB"]["none"] += 1
                continue
            end

            # update counts
            cb_i = get!(cb_dictionary, cb, length(cb_dictionary) + 1)
            umi_i = UMItoindex(umi)
            key = (cb_i, umi_i, sb_i)
            mat[key] = get(mat, key, 0) + 1
            metadata["reads_filtered"] += 1

            if metadata["reads_filtered"] % 10_000_000 == 0
                println(metadata["reads_filtered"]) ; flush(stdout)
            end
        end
    end
    
    # Turn matrix from dictionary into dataframe
    df = DataFrame(cb_i = UInt32[], umi_i = UInt32[], sb_i = UInt32[], reads = UInt32[])
    for key in keys(mat)
        value = pop!(mat, key)
        push!(df, (key[1], key[2], key[3], value))
    end

    # Turn cb_dictionary into cb_whitelist
    cb_whitelist = DataFrame(cb = collect(String31, keys(cb_dictionary)), cb_i = collect(UInt32, values(cb_dictionary)))
    sort!(cb_whitelist, :cb_i)
    @assert cb_whitelist.cb_i == 1:size(cb_whitelist, 1)

    return(df, cb_whitelist.cb, metadata)
end

df, cb_whitelist, metadata = process_fastqs(R1s, R2s)

@assert metadata["UP"]["exact"] + metadata["UP"]["fuzzy"] == sum(values(metadata["SB"]))
@assert metadata["SB"]["exact"] + metadata["SB"]["HD1"] == metadata["reads_filtered"] == sum(df.reads)
@assert sum(values(metadata["SB_HD"])) == metadata["SB"]["HD1"]

println("...done") ; flush(stdout) ; GC.gc()

# Create a downsampling curve
downsampling = UInt32[]
table = countmap(df.reads)
for prob in 0:0.05:1
    s = [length(unique(floor.(sample(0:k*v-1, round(Int,k*v*prob), replace=false)/k))) for (k,v) in zip(keys(table),values(table))]
    append!(downsampling, sum(s))
    GC.gc()
end

##### Save results #############################################################

print("Saving results... ") ; flush(stdout)

h5open("SBcounts.h5", "w") do file
    create_group(file, "lists")
    file["lists/cb_list", compress=9] = String31.(cb_whitelist) # Vector{String}
    file["lists/sb_list", compress=9] = String31.(sb_whitelist) # Vector{String}
    file["lists/puck_list"] = String255.(basename.(puck_paths)) # Vector{String}
    file["lists/R1_list"] = String255.(basename.(R1s))          # Vector{String}
    file["lists/R2_list"] = String255.(basename.(R2s))          # Vector{String}

    create_group(file, "matrix")
    file["matrix/cb_index", compress=9] = df.cb_i # Vector{UInt32}
    file["matrix/umi", compress=9] = df.umi_i     # Vector{UInt32}
    file["matrix/sb_index", compress=9] = df.sb_i # Vector{UInt32}
    file["matrix/reads", compress=9] = df.reads   # Vector{UInt32}

    create_group(file, "puck")
    file["puck/sb", compress=9] = String31.(puckdf.sb)      # Vector{String}
    file["puck/x", compress=9] = puckdf.x                   # Vector{Float64}
    file["puck/y", compress=9] = puckdf.y                   # Vector{Float64}
    file["puck/puck_index", compress=9] = puckdf.puck_index # Vector{UInt8}
    
    create_group(file, "metadata")
    
    file["metadata/switch"] = convert(Int8, switch)
    file["metadata/bead"] = String7(bead)
    file["metadata/num_lowQbeads"] = num_lowQbeads
    
    file["metadata/reads"] = metadata["reads"]
    file["metadata/reads_filtered"] = metadata["reads_filtered"]
    file["metadata/R1_tooshort"] = metadata["R1_tooshort"]
    file["metadata/R2_tooshort"] = metadata["R2_tooshort"]

    create_group(file, "metadata/UMI")
    file["metadata/UMI/type"] = keys(metadata["UMI"]) |> x -> String15.(x) |> collect
    file["metadata/UMI/count"] = values(metadata["UMI"]) |> collect
    
    create_group(file, "metadata/UP")
    file["metadata/UP/type"] = keys(metadata["UP"]) |> x -> String15.(x) |> collect
    file["metadata/UP/count"] = values(metadata["UP"]) |> collect

    create_group(file, "metadata/SB")
    file["metadata/SB/type"] = keys(metadata["SB"]) |> x -> String15.(x) |> collect
    file["metadata/SB/count"] = values(metadata["SB"]) |> collect

    create_group(file, "metadata/SB_HD")
    file["metadata/SB_HD/type"] = keys(metadata["SB_HD"]) |> collect
    file["metadata/SB_HD/count"] = values(metadata["SB_HD"]) |> collect
    
    file["metadata/downsampling"] = downsampling # Vector{UInt32}
end;

println("done") ; flush(stdout) ; GC.gc()

##### Documentation ############################################################

# The workflow for processing the FASTQs is as follows:
#   1) UMI filter: discard reads that have an N in the umi or a homopolymer umi
#   2) UP filter: discard reads that don't have a detectable UP site
#      - UP site does fuzzy matching - allow 1 mismatch for every 6 bases
#   3) SB filter: discard reads whose spatial barcode (sb) doesn't match the puck whitelist
#      - does fuzzy matching with 1 mismatch
#   4) matrix update
#      - for reads that made it this far, the (cb_i, umi_i, sb_i) key of the dictionary will be incremented by 1
#      - cb_i, sb_i represent an index into the cb_whitelist, sb_whitelist vectors respectively
#      - umi_i is the 2-bit encoded umi