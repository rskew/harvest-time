using JuMP
using Gurobi

include("read_data.jl")

#filename = ARGS[1]

# Dummy data for model development

properties = 1:length(initialAge)
minimum_yield_per_year = 500_000
maximum_yield_per_year = 1_500_000
n_properties = 199
n_years = 8
n_future_chunks = 4
n_chunk_spacing = 3
replanting_cost = 1000
time_limit = 120 # seconds


properties = 1:n_properties
initial_ages = initialAge
property_age = 1:length(density[1,:])
years = 1:n_years
years_plus_chunks = 1:(n_years+n_future_chunks)
future_chunks = 1:n_future_chunks
M = 100000000

years = 1:1
minimum_yield_per_year = 5
maximum_yield_per_year = 100


function print_harvests(harvests)
    for p in properties
        for t in years
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

#m = Model(solver=CbcSolver(log=1, Sec=300))
m = Model(solver=GurobiSolver(TimeLimit=time_limit))
# 'harvest' is 1 for a property for a given year if that property is
# harvested that year
@variable(m, harvest[properties,years], Bin)
@variable(m, future_harvest[properties,future_chunks], Bin)

# 'age' gives the number of years since the previous harvest of a property
@variable(m, age[properties,years], Int)
@variable(m, future_age[properties,future_chunks], Int)

# harvest_age stores the age of a property on the years it is harvested.
@variable(m, harvest_age[properties,years], Int)
@variable(m, future_harvest_age[properties,future_chunks], Int)

# Linearise by parts the age vs yield function.
# Each 'harvest_age[p,t]' has its own linearised curve variables.
# The 'y's from the lecture notes:
@variable(m, y[properties, years_plus_chunks, property_age-1], Bin)
# The lambdas from the lecture notes:
@variable(m, 0 <= lambda[properties, years_plus_chunks, property_age] <= 1)

# Maximise profit: balance yield against replanting cost
@objective(m, Max,
           sum(lambda[p,t,pa] * yields[p,pa] * netReturns[p]
               for p in properties,
               t in years_plus_chunks,
               pa in property_age)
           - sum(harvest[p,t] * replanting_cost * area[p]
                 for p in properties,
                 t in years)
           - sum(future_harvest[p,c] * replanting_cost * area[p]
                 for p in properties,
                 t in future_chunks))

# Integer variables >= 0
@constraint(m, [p in properties, t in years],
            age[p,t] >= 0)
@constraint(m, [p in properties, t in years],
            future_age[p,t] >= 0)
@constraint(m, [p in properties, t in years],
            harvest_age[p,t] >= 0)
@constraint(m, [p in properties, t in years],
            future_harvest_age[p,t] >= 0)

# Linearise by parts constraints
@constraint(m, [p in properties, t in years],
            harvest_age[p,t] == sum(lambda[p,t,pa] * pa
                                    for pa in property_age))
@constraint(m, [p in properties, c in future_chunks],
            future_harvest_age[p,c] == sum(lambda[p,c+n_years,pa] * pa
                                           for pa in property_age))

@constraint(m, [p in properties, t in years],
            sum(lambda[p,t,pa] for pa in property_age) == harvest[p,t])
@constraint(m, [p in properties, c in future_chunks],
            sum(lambda[p,c+n_years,pa] for pa in property_age) == future_harvest[p,c])

@constraint(m, [p in properties, t in years_plus_chunks],
            sum(y[p,t,pa] for pa in property_age[1:end-1]) == 1)

@constraint(m, [p in properties, t in years_plus_chunks],
            lambda[p,t,1] <= y[p,t,1])

@constraint(m, [p in properties, t in years_plus_chunks, pa in property_age[1:end-2]],
            lambda[p,t,pa+1] <= y[p,t,pa] + y[p,t,pa+1])

@constraint(m, [p in properties, t in years_plus_chunks],
            lambda[p,t,length(property_age)] <= y[p,t,length(property_age)-1])

# harvest_age is zero for years where there is no harvest
@constraint(m, [p in properties, t in years],
            harvest_age[p,t] <= M * harvest[p,t])
@constraint(m, [p in properties, c in future_chunks],
            future_harvest_age[p,c] <= M * future_harvest[p,c])
# harvest_age is the age of the property when it is harvested.
# To maximise the objective, the solver will automatically push
# harvest_age up to the age of the property on harvest years.
@constraint(m, [p in properties, t in years],
            harvest_age[p,t] <= age[p,t])
@constraint(m, [p in properties, c in future_chunks],
            future_harvest_age[p,c] <= future_age[p,c])
# Force 'harvest_age[p,t]' to be 'age[p,t]' when harvest[p,t]
@constraint(m, [p in properties, t in years],
            harvest_age[p,t] >= age[p,t] - M * (1 - harvest[p,t]))
@constraint(m, [p in properties, c in chunks],
            future_harvest_age[p,t] >= future_age[p,t] - M * (1 - harvest[p,t]))

# Assign initial ages
@constraint(m, [p in properties],
            age[p,1] == initial_ages[p] + 1)
# 'age' must be 1 for the year after a harvest. For any other year, 'age' must
# be the previous age +1.
@constraint(m, [p in properties, t in years],
            age[p, t] >= 1)
@constraint(m, [p in properties, prev_t in years, next_t in years; prev_t+1 == next_t],
            age[p, next_t] <= age[p, prev_t] + 1)

@constraint(m, [p in properties, prev_t in years, next_t in years; prev_t+1 == next_t],
            age[p, next_t] <= 1 + M * (1 - harvest[p,prev_t]))

@constraint(m, [p in properties, prev_t in years, next_t in years; prev_t+1 == next_t],
            age[p, next_t] >= age[p, prev_t] + 1 - M * harvest[p,prev_t])


# Must provide a minimum quantity per year.
@constraint(m, [t in years],
            sum(lambda[p,t,pa] * yields[p,pa]
                for p in properties,
                pa in property_age)
            >= minimum_yield_per_year)

# The amount of product that can be processed is limited.
@constraint(m, [t in years],
            sum(lambda[p,t,pa] * yields[p,pa]
                for p in properties,
                pa in property_age)
            <= maximum_yield_per_year)

solve(m)

println("harvest_ages[p,t]:")
print_harvests(harvest_age)
println("yields[p,t]:")
print_yields(yields)
println("yields per year:")
for t in years
    print("$(sum(getvalue(lambda[p,t,pa]*yields[p,pa]) for p in properties, pa in property_age)) ")
end
println()


println("Total cost $(getobjectivevalue(m))")
println("Bound is $(getobjectivebound(m))")

include("plot_test.jl")
