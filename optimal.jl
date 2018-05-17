using JuMP
using Cbc

#include("readerFunction.jl")
#include("writeSolutionToFile.jl")
#include("validateSolution.jl")

#filename = ARGS[1]


# dummy data for model development
properties = 1:5
initial_ages = [1,1,1,1,1]
nYears = 5
years = 1:nYears
yield_age_slope = 1
yield_age_function(age) = age * yield_age_slope

### Modelling

m = Model(solver=CbcSolver(log=3, Sec=30))

# 'harvest' is 1 for a property for a given year if that property is
# harvested that year
@variable(m, harvest[properties, years], Bin)
# 'age' gives the number of years since the previous harvest of a property
#@variable(m, age[properties, years], Int)

# Maximise profit = yield*value_per_unit - costs
@objective(m, Max, sum(harvest_sequence_values[hs] * harvest_sequence_choice[p]
                       for hs in harvest_sequences, p in properties))

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
            age[p, next_y] >= age[p, prev_y] + 1 - M * harvest[p,y])

#solve(m)
#
#for p in properties
#    for y in years
#  	    println("harvest[$p,$y]: ",getvalue(harvest[p,y]))
#  	    println("age[$p,$y]: ",getvalue(age[p,y]))
#    end
#end
