using JuMP
using Cbc

#include("readerFunction.jl")
#include("writeSolutionToFile.jl")
#include("validateSolution.jl")

#filename = ARGS[1]

# Dummy data for model development
properties = 1:20
initial_ages = ones(length(properties),1)
slopes = [1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]
plateau_years = 4*ones(length(properties),1)
minimum_yield_per_year = 5
maximum_yield_per_year = 100
years = 1:5
M = 100000000

yields = zeros(length(properties),length(years))
for p in properties
    for y in years
        current_age = initial_ages[p] + y
        if current_age < plateau_years[y]
            yields[p,y] = slopes[p] * current_age
        else
            yields[p,y] = slopes[p] * plateau_years[p]
        end
    end
end


function print_harvests(harvests)
    for p in properties
        for y in years
            print("$(Int(round(getvalue(harvests[p,y])))) ")
        end
        println()
    end
end

### Modelling

m = Model(solver=CbcSolver(log=3, Sec=30))

# 'harvest' is 1 for a property for a given year if that property is
# harvested that year
@variable(m, harvest[properties,years], Bin)

# 'age' gives the number of years since the previous harvest of a property
@variable(m, age[properties,years], Int)

# harvest_age stores the age of a property on the years it is harvested.
@variable(m, harvest_age[properties,years], Int)

@objective(m, Max, sum(harvest_age[p,y] * yields[p,y] for p in properties, y in years))

# All variables greater than zero
@constraint(m, [p in properties, y in years],
            harvest[p,y] >= 0)
@constraint(m, [p in properties, y in years],
            age[p,y] >= 0)
@constraint(m, [p in properties, y in years],
            harvest_age[p,y] >= 0)

# harvest_age is zero for years where there is no harvest
@constraint(m, [p in properties, y in years],
            harvest_age[p,y] <= M * harvest[p,y])
# harvest_age is the age of the property when it is harvested.
# To maximise the objective, the solver will automatically push
# harvest_age up to the age of the property on harvest years.
@constraint(m, [p in properties, y in years],
            harvest_age[p,y] <= age[p,y])
# Force 'harvest_age[p,y]' to be 'age[p,y]' when harvest[p,y]
@constraint(m, [p in properties, y in years],
            harvest_age[p,y] >= age[p,y] - M * (1 - harvest[p,y]))

# Only harvest each property once
@constraint(m, [p in properties],
            sum(harvest[p,y] for y in years) <= 1)


# assign initial ages
@constraint(m, [p in properties],
            age[p,1] == initial_ages[1])
# 'age' must be 1 for the year after a harvest. For any other year, 'age' must
# be the previous age +1.
@constraint(m, [p in properties, y in years],
            age[p, y] >= 1)
@constraint(m, [p in properties, prev_y in years, next_y in years; prev_y+1 == next_y],
            age[p, next_y] <= age[p, prev_y] + 1)

@constraint(m, [p in properties, prev_y in years, next_y in years; prev_y+1 == next_y],
            age[p, next_y] <= 1 + M * (1 - harvest[p,prev_y]))

@constraint(m, [p in properties, prev_y in years, next_y in years; prev_y+1 == next_y],
            age[p, next_y] >= age[p, prev_y] + 1 - M * harvest[p,prev_y])


# Must provide a minimum quantity per year.
@constraint(m, [y in years],
            sum(harvest_age[p,y] * slopes[p]
                for p in properties)
            >= minimum_yield_per_year)

# The amount of product that can be processed is limited.
@constraint(m, [y in years],
            sum(harvest_age[p,y] * slopes[p]
                for p in properties)
            <= maximum_yield_per_year)

solve(m)

println("harvest[p,y]:")
print_harvests(harvest)
println("age[p,y]:")
print_harvests(age)
println("harvest_ages[p,y]:")
print_harvests(harvest_age)

println("Total cost $(getobjectivevalue(m))")
println("Bound is $(getobjectivebound(m))")
