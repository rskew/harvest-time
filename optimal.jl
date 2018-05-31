using JuMP
using Gurobi

include("read_data.jl")

properties = 1:length(initialAge)
minimum_yield_per_year = 500_000
maximum_yield_per_year = 1_500_000
n_properties = 199
n_years = 30
replanting_cost = 1000
discount_rate = 0.0
time_limit = 600 # seconds


properties = 1:n_properties
initial_ages = initialAge
property_age = 1:length(yields[1,:])
years_to_plan = 1:n_years
M = 10000


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
        for t in 1:length(yields[1,:])
            print("$(yields[p,t]) ")
        end
        println()
    end
end

### Modelling

m = Model(solver=GurobiSolver(TimeLimit=time_limit))

# 'harvest' is 1 for a property for a given year if that property is
# harvested that year
heuristic_harvest_choices = readcsv("heuristic_harvest_choices.csv")
@variable(m, harvest[p=properties,y=years_to_plan], Bin, start=heuristic_harvest_choices[p,y])

# 'age' gives the number of years_to_plan since the previous harvest of a property
heuristic_ages = readcsv("heuristic_ages.csv")
@variable(m, age[p=properties,y=years_to_plan], Int, start=heuristic_ages[p,y])

# HARVEST_age stores the age of a property on the years_to_plan it is harvested.
heuristic_harvest_ages = readcsv("heuristic_harvest_ages.csv")
@variable(m, harvest_age[p=properties,y=years_to_plan], Int, start=heuristic_harvest_ages[p,y])

# Linearise by parts the age vs yield function.
# Each 'harvest_age[p,t]' has its own linearised curve variables.
# The 'y's from the lecture notes:
@variable(m, y[properties, years_to_plan, property_age-1], Bin)
# The lambdas from the lecture notes:
@variable(m, 0 <= lambda[properties, years_to_plan, property_age] <= 1)

# Maximise profit: balance yield against replanting cost
@objective(m, Max,
           sum(lambda[p,t,pa] * yields[p,pa] * netReturns[p]
               * (1 - discount_rate)^(t-1)
               for p in properties,
               t in years_to_plan,
               pa in property_age)
           - sum(harvest[p,t] * replanting_cost * area[p]
                 * (1 - discount_rate)^(t-1)
                 for p in properties,
                 t in years_to_plan))

# Integer variables >= 0
@constraint(m, [p in properties, t in years_to_plan],
            age[p,t] >= 0)
@constraint(m, [p in properties, t in years_to_plan],
            harvest_age[p,t] >= 0)

# Linearise by parts constraints
@constraint(m, [p in properties, t in years_to_plan],
            harvest_age[p,t] == sum(lambda[p,t,pa] * pa
                                    for pa in property_age))
@constraint(m, [p in properties, t in years_to_plan],
            sum(lambda[p,t,pa] for pa in property_age) == harvest[p,t])
@constraint(m, [p in properties, t in years_to_plan],
            sum(y[p,t,pa] for pa in property_age[1:end-1]) == 1)
@constraint(m, [p in properties, t in years_to_plan],
            lambda[p,t,1] <= y[p,t,1])
@constraint(m, [p in properties, t in years_to_plan, pa in property_age[1:end-2]],
            lambda[p,t,pa+1] <= y[p,t,pa] + y[p,t,pa+1])
@constraint(m, [p in properties, t in years_to_plan],
            lambda[p,t,length(property_age)] <= y[p,t,length(property_age)-1])

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
            age[p,1] == initial_ages[p] + 1)
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
            sum(lambda[p,t,pa] * yields[p,pa]
                for p in properties,
                pa in property_age)
            >= minimum_yield_per_year)

# The amount of product that can be processed is limited.
@constraint(m, [t in years_to_plan],
            sum(lambda[p,t,pa] * yields[p,pa]
                for p in properties,
                pa in property_age)
            <= maximum_yield_per_year)

solve(m)

#println("harvest[p,t]:")
#print_harvests(harvest)
#println("age[p,t]:")
#print_harvests(age)


#println("harvest_ages[p,t]:")
#print_harvests(harvest_age)
#println("yields[p,t]:")
#print_yields(yields)
#println("yields per year:")
#for t in years_to_plan
#    print("$(sum(getvalue(lambda[p,t,pa]*yields[p,pa]) for p in properties, pa in property_age)) ")
#end
#println()

#println("lambda:")
#for pa in property_age
#    for t in years_to_plan
#        print("$(getvalue(lambda[length(properties),t,pa])) ")
#    end
#    println()
#end
#
#println("y:")
#for pa in property_age[1:end-1]
#    for t in years_to_plan
#        print("$(Int(round(getvalue(y[length(properties),t,pa])))) ")
#    end
#    println()
#end


println("Total cost $(getobjectivevalue(m))")
println("Bound is $(getobjectivebound(m))")

years = years_to_plan
include("plot_test.jl")
