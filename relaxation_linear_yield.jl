using JuMP
#using Cbc
using Gurobi

include("read_data.jl")

#filename = ARGS[1]

# Dummy data for model development
minimum_yield_per_year = 500
maximum_yield_per_year = 2000
n_properties = 199
n_years = 21
replanting_cost = 10
time_limit = 60 # seconds


properties = 1:n_properties
initial_ages = initialAge[1:n_properties] + 1
property_age = 1:length(density[1,:])
years_to_plan = 1:n_years
M = 100000000


function print_harvests(harvests)
    for p in properties
        for t in years_to_plan
            print("$(Int(round(getvalue(harvests[p,t])))) ")
        end
        println()
    end
end
function print_yields(yields)
    for p in properties
        for t in 1:length(density[1,:])
            print("$(yields[p,t]) ")
        end
        println()
    end
end

### Modelling

#m = Model(solver=CbcSolver(log=1, Sec=time_limit))
m = Model(solver=GurobiSolver(TimeLimit=time_limit))

# 'harvest' is 1 for a propertt for a given year if that property is
# harvested that year
@variable(m, harvest[properties,years_to_plan], Bin)

# 'age' gives the number of years_to_plan since the previous harvest of a property
@variable(m, age[properties,years_to_plan], Int)

# HARVEST_age stores the age of a property on the years_to_plan it is harvested.
@variable(m, harvest_age[properties,years_to_plan], Int)

# Linearise by parts the age vs yield function.
# Each 'harvest_age[p,t]' has its own linearised curve variables.
# The 'y's from the lecture notes:
@variable(m, y[properties, years_to_plan, property_age-1], Bin)
# The lambdas from the lecture notes:
@variable(m, 0 <= lambda[properties, years_to_plan, property_age] <= 1)

# Maximise profit: balance yield against replanting cost
@objective(m, Max,
           sum(harvest_age[p,t] * growthRate[p]
               - harvest[p,t] * replanting_cost
               for p in properties,
               t in years_to_plan))

# Integer variables >= 0
@constraint(m, [p in properties, t in years_to_plan],
            age[p,t] >= 0)
@constraint(m, [p in properties, t in years_to_plan],
            harvest_age[p,t] >= 0)

# harvest_age is zero for years_to_plan where there is no harvest
@constraint(m, [p in properties, t in years_to_plan],
            harvest_age[p,t] <= M * harvest[p,t])
# harvest_age is the age of the property when it is harvested.
# To maximise the objective, the solver will automatically push
# harvest_age up to the age of the property on harvest years_to_plan.
@constraint(m, [p in properties, t in years_to_plan],
            harvest_age[p,t] <= age[p,t])
# Force 'harvest_age[p,t]' to be 'age[p,t]' when harvest[p,t]
@constraint(m, [p in properties, t in years_to_plan],
            harvest_age[p,t] >= age[p,t] - M * (1 - harvest[p,t]))

# Assign initial ages
@constraint(m, [p in properties],
            age[p,1] == initial_ages[p])
# 'age' must be 1 for the year after a harvest. For any other year, 'age' must
# be the previous age +1.
@constraint(m, [p in properties, t in years_to_plan],
            age[p, t] >= 1)
@constraint(m, [p in properties, prev_t in years_to_plan, next_t in years_to_plan; prev_t+1 == next_t],
            age[p, next_t] <= age[p, prev_t] + 1)

@constraint(m, [p in properties, prev_t in years_to_plan, next_t in years_to_plan; prev_t+1 == next_t],
            age[p, next_t] <= 1 + M * (1 - harvest[p,prev_t]))

@constraint(m, [p in properties, prev_t in years_to_plan, next_t in years_to_plan; prev_t+1 == next_t],
            age[p, next_t] >= age[p, prev_t] + 1 - M * harvest[p,prev_t])


# Must provide a minimum quantity per year.
@constraint(m, [t in years_to_plan],
            sum(harvest_age[p,t] * growthRate[p]
                for p in properties)
            >= minimum_yield_per_year)

# The amount of product that can be processed is limited.
@constraint(m, [t in years_to_plan],
            sum(harvest_age[p,t] * growthRate[p]
                for p in properties)
            <= maximum_yield_per_year)

solve(m)

println("harvest[p,t]:")
print_harvests(harvest)
println("age[p,t]:")
print_harvests(age)
println("harvest_ages[p,t]:")
print_harvests(harvest_age)
println("growth rates:")
println("$(growthRate[1:n_properties])")
println("yields:")
for t in years_to_plan
    print("$(getvalue(sum(harvest_age[p,t] * growthRate[p]
            for p in properties))) ")
end
println()
println("Total cost $(getobjectivevalue(m))")
println("Bound is $(getobjectivebound(m))")
