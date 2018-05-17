########
# Heuristic model:
# Relax the problem by searching for a time to harvest each property once,
# then fix the harvests chosen for the next year. Slide the window forward,
# updating the initial ages, and repeat.
########

using JuMP
using Cbc

#include("readerFunction.jl")
#include("writeSolutionToFile.jl")
#include("validateSolution.jl")

#filename = ARGS[1]


# Dummy data for model development
properties = 1:23
initial_ages = zeros(length(properties),1)
slopes = [1,2,3,2,1,5,5,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]
plateau_years = 4*ones(length(properties),1)
minimum_yield_per_year = 3
maximum_yield_per_year = 1000000

# Parameters to tune
years_window = 1:10
years_to_plan = 1:20

tie_breaker_slope = -0.001

function calc_yields(initial_ages, slopes, plateau_years)
    yields = zeros(length(properties),length(years_window))
    for p in properties
        for y in years_window
            current_age = initial_ages[p] + y
            if current_age < plateau_years[y]
                yields[p,y] = slopes[p] * current_age
            else
                yields[p,y] = slopes[p] * plateau_years[p] +
                    tie_breaker_slope * (y - plateau_years[p])
            end
        end
    end
    return yields
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

yields = calc_yields(initial_ages, slopes, plateau_years)

harvest_choices = zeros(length(properties),length(years_to_plan))


for year in years_to_plan
    m = Model(solver=CbcSolver(log=3, Sec=30))

    @variable(m, harvest[properties,years_window], Bin)

    @objective(m, Max, sum(harvest[p,y] * yields[p,y] for p in properties, y in years_window))

    # Only havest each property once
    @constraint(m, [p in properties],
                sum(harvest[p,y] for y in years_window) <= 1)

    # Must provide a minimum quantity per year.
    # Otherwise properties would be harvested at the latest possible date.
    @constraint(m, [y in years_window],
                sum(harvest[p,y] * yields[p,y] for p in properties) >= minimum_yield_per_year)

    # The amount of product that can be processed is limited.
    @constraint(m, [y in years_window],
                sum(harvest[p,y] * yields[p,y] for p in properties) <= maximum_yield_per_year)

    solve(m)

    println()
    println()
    println("Next Iteration:")
    println()
    println("initial_ages:")
    println("$initial_ages")
    println()
    print_yields(yields)
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

    yields = calc_yields(initial_ages, slopes, plateau_years)
end

println()
println()
println()
println()
print_harvest_choices(harvest_choices)
