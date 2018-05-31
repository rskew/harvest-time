########
# Heuristic model:
# Relax the problem by searching for a time to harvest each property once,
# then fix the harvests chosen for the next year. Slide the window forward,
# updating the initial ages, and repeat.
########

using JuMP
using Gurobi

include("read_data.jl")

minimum_yield_per_year = 500_000
maximum_yield_per_year = 1_250_000
n_properties = length(initialAge)
properties = 1:n_properties
replanting_cost = 1000
time_limit = 30 # seconds

property_age = 1:length(yields[1,:])
initial_ages = initialAge + 1

# Parameters to tune
n_years_window = 15
n_years_to_plan = 50
discount_rate = 0.0
yield_decay_rate = 0.05
decay_start = 10
years_window = 1:n_years_window
years_to_plan = 1:n_years_to_plan

# Extend plateau section. Some properties may not be harvested until later.
n_yields_to_fill = (n_years_window + n_years_to_plan)
if n_yields_to_fill > 0
    plateau = repeat(yields[:,end],outer=[1,n_yields_to_fill])
    yields_extended_plateau = hcat(yields, plateau)
end

# Add yield decay
yields_decayed = zeros(size(yields_extended_plateau))
for y in 1:size(yields_extended_plateau,2)
    if y < decay_start
        yields_decayed[:,y] = yields_extended_plateau[:,y]
    else
        yields_decayed[:,y] = yields_extended_plateau[:,y] * (1 - yield_decay_rate)^(y-decay_start)
    end
end


function print_harvests(harvests)
    println("harvests[p,y]:")
    for p in properties
        for y in years_window
            print("$(Int(round(getvalue(harvests[p,y])))) ")
        end
        println()
    end
end
function print_harvest_choices(harvest_choices)
    println("harvest_choices:")
    for p in properties
        for year in years_to_plan
            print("$(Int(round(harvest_choices[p,year]))) ")
        end
        println()
    end
end
function print_yields(yields)
    println("yields[p,y]:")
    for p in properties
        for y in years_window
            print("$(yields[p,y]) ")
        end
        println()
    end
end

#yields = calc_yields(initial_ages, slopes, plateau_years)

harvest_choices = zeros(length(properties),length(years_to_plan))


for year in years_to_plan
    println("Year: $(year)")
    m = Model(solver=GurobiSolver(TimeLimit=time_limit))

    @variable(m, harvest[properties,years_window], Bin)

    @objective(m, Max, sum(harvest[p,y] * yields_decayed[p,y+initial_ages[p]] * (1 - discount_rate)^(y-1) for p in properties, y in years_window))

    # Only harvest each property once
    @constraint(m, [p in properties],
                sum(harvest[p,y] for y in years_window) <= 1)

    # Must provide a minimum quantity per year.
    # Otherwise properties would be harvested at the latest possible date.
    @constraint(m, [y in years_window],
                sum(harvest[p,y] * yields_decayed[p,y+initial_ages[p]] for p in properties) >= minimum_yield_per_year)

    # The amount of product that can be processed is limited.
    @constraint(m, [y in years_window],
                sum(harvest[p,y] * yields_decayed[p,y+initial_ages[p]] for p in properties) <= maximum_yield_per_year)

    solve(m)

    println()
    println()
    println("Next Iteration:")
    println()
    println("initial_ages:")
    println("$initial_ages")
    println()
    print_harvests(harvest)
    println()
    println("Total cost $(getobjectivevalue(m))")
    println("Bound is $(getobjectivebound(m))")
    harvest_choices[:,year] = getvalue(harvest[:,1])


    # Update initial ages
    for p in properties
        if getvalue(harvest[p,1]) > 0.5
            initial_ages[p] = 0
        else
            initial_ages[p] = initial_ages[p] + 1
        end
    end

end

println()
println()
println()
println()
print_harvest_choices(harvest_choices)

ages = zeros(size(harvest_choices))
ages[:,1] = initialAge[properties,:] + 1
for y in years_to_plan[1:end-1]
    for p in properties
        if harvest_choices[p,y] > 0.5
            ages[p,y+1] = 1
        else
            ages[p,y+1] = ages[p,y] + 1
        end
    end
end

harvest_ages = ages .* harvest_choices

println("harvest_choices")
print_harvest_choices(harvest_choices)
println("ages")
print_harvest_choices(ages)
println("harvest_ages")
print_harvest_choices(harvest_ages)

heuristic_harvest_choices = round.(Int,harvest_choices)
heuristic_ages = round.(Int,ages)
heuristic_harvest_ages = round.(Int,harvest_ages)

writecsv("heuristic_harvest_choices.csv",heuristic_harvest_choices)
writecsv("heuristic_ages.csv",heuristic_ages)
writecsv("heuristic_harvest_ages.csv",heuristic_harvest_ages)

include("plot_test_plan_next_year.jl")
