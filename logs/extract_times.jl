using Statistics

# Get the directory of the file
cd(dirname(Base.source_path()))

# Open the log file
file = open("20240407-fit.log")

# Initialize an array to store the time values
times = Float64[]

# Read the file line by line
for line in eachline(file)
    # Use a regular expression to match lines that contain the time information
    matched = match(r"(\d+\.\d+) secs", line)
    if matched !== nothing
        # Extract the time value and convert it to a float
        time = parse(Float64, matched[1])
        push!(times, time)
    end
end

# Calculate the differences between consecutive time values
differences = diff(times)

# remove any negative differences
differences = filter(x -> x >= 0, differences)

mean_diff = mean(differences)
println("The mean bewteen iterations is: $mean_diff")
println("The standard deviation bewteen iterations is: $(std(differences))")

# Close the file
close(file)