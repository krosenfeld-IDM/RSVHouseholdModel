using JLD

files = readdir("input")

# # remove DailyIFATbyAgeCat.jld from files
# files = filter(x -> x != "DailyIFATbyAgeCat.jld", files)

# # remove daily_RSV_hosp.jld from files
# files = filter(x -> x != "daily_RSV_hosp.jld", files)

# loop through files and load them
for file in files
    println("Loading $file")
    data = load("input/$file")
    # println(typeof(data)) # everythign is a dict 
    for (key, value) in data
        println("Key: $key, Type: $(typeof(value))")
    end
    println()
end

using Roots
f(x) = exp(x) - x^4;
α₀, α₁, α₂ = -0.8155534188089607, 1.4296118247255556, 8.6131694564414;
fzero(f, (8,9), Roots.Bisection()) ≈ α₂ # a bisection method has the bracket specified