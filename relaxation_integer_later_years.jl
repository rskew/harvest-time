using JuMP
using Gurobi

include("read_data.jl")

properties = 1:length(initialAge)
#minimum_yield_per_year = 500_000
#maximum_yield_per_year = 1_500_000
#n_properties = 199
minimum_yield_per_year = 50_0
maximum_yield_per_year = 150_00000
n_properties = 20
n_years_int = 5
n_years_real = 5
#n_real_harvests = n_years_real
n_real_harvests = 7
n_years_total = n_years_int + n_years_real
replanting_cost = 1000
time_limit = 30 # seconds


properties = 1:n_properties
initial_ages = initialAge
property_age = 1:length(yields[1,:])
real_harvests = 1:n_real_harvests
years_int = 1:n_years_int
years_real = 1:n_years_real
years_total = 1:n_years_total
harvests = 1:(n_years_int+n_real_harvests)
M = 100000000


function print_harvests(harvests)
    for p in properties
        for t in years_int
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
@variable(m, harvest[properties,years_int], Bin)

# 'age' gives the number of years since the previous harvest of a property
@variable(m, age[properties,years_int], Int)

# HARVEST_age stores the age of a property on the years it is harvested.
@variable(m, harvest_age[properties,years_int], Int)

########################
# Relaxation: add a set of real valued harvest age variables for each property
# for years after the integer planning period
# eg integer harvest ages and years up to 5 years of planning,
# then relaxed real values for ages and harvest dates from that to 30 years

# Whether or not to use this real harvest variable
@variable(m, real_harvest[properties,real_harvests], Bin)

# Which real valued time in years this harvest happens
@variable(m, n_years_int <= real_harvest_year[properties,real_harvests] <= n_years_total)

# The age of the property at the harvest time in years
@variable(m, real_harvest_age[properties,real_harvests])


# Bunch the real harvests at the start of the array.
# When they stop being used, there can be no more harvests to the right
@constraint(m, [p in properties, r in real_harvests ; r < n_real_harvests],
            real_harvest[p,r+1] <= real_harvest[p,r])

# Real harvests are listed in consecutive order
@constraint(m, [p in properties, r in real_harvests ; r < n_real_harvests],
            real_harvest_year[p,r+1] >= real_harvest_year[p,r])

# First real harvest age is the age at the end of the integer period plus the time since then
@constraint(m, [p in properties],
            real_harvest_age[p,1] == age[p,n_years_int] + (real_harvest_year[p,1] - n_years_int))
# Calculate real_harvest_age from real_harvest_year
@constraint(m, [p in properties, r in real_harvests ; r >= 2],
            real_harvest_age[p,r] == real_harvest_year[p,r] - real_harvest_year[p,r-1])


########################
# Linearise by parts the age vs yield function.
# Each 'harvest_age[p,t]' has its own linearised curve variables.
# The 'y's from the lecture notes:
@variable(m, y[properties, harvests, property_age[1:end-1]], Bin)
# The lambdas from the lecture notes:
@variable(m, 0 <= lambda[properties, harvests, property_age] <= 1)

# (Int) harvest age interpolates two lambdas
@constraint(m, [p in properties, t in years_int],
            harvest_age[p,t] == sum(lambda[p,t,pa] * pa
                                    for pa in property_age))
# (Real) harvest age interpolated from lambdas
@constraint(m, [p in properties, r in real_harvests],
            real_harvest_age[p,r] == sum(lambda[p,r+n_years_int,pa] * pa
                                    for pa in property_age))

# (Int) Lambdas must sum to 1 to be valid interpolation
@constraint(m, [p in properties, t in years_int],
            sum(lambda[p,t,pa] for pa in property_age) == harvest[p,t])
# (Real) Lambdas must sum to 1 to be valid interpolation
@constraint(m, [p in properties, r in real_harvests],
            sum(lambda[p,r+n_years_int,pa] for pa in property_age) == real_harvest[p,r])

# (Int) Only choose 1 pair of consecutive lambdas using y variable to signal which pair
@constraint(m, [p in properties, t in harvests],
            sum(y[p,t,pa] for pa in property_age[1:end-1]) == 1)
# (Real) Only choose 1 pair of consecutive lambdas using y variable to signal which pair
@constraint(m, [p in properties, r in real_harvests],
            sum(y[p,r+n_years_int,pa] for pa in property_age[1:end-1]) == real_harvest[p,r])
# Left endpoint
@constraint(m, [p in properties, t in harvests],
            lambda[p,t,1] <= y[p,t,1])
# Midpoints
@constraint(m, [p in properties, t in harvests, pa in property_age[1:end-2]],
            lambda[p,t,pa+1] <= y[p,t,pa] + y[p,t,pa+1])
# Right endpoint
@constraint(m, [p in properties, t in harvests],
            lambda[p,t,length(property_age)] <= y[p,t,length(property_age)-1])


########################
# Need var to keep track of the discrete year that a property was harvested
# so that max/min harvest constraints can be enforced each year.
# This is done with another linearisation by parts type method, with lambda
# and y vars allowing a real harvest age to be allocated to a bucket of
# harvests that happen in a certain year.
@variable(m, real_harvest_year_bucket[properties,real_harvests,years_real[1:end-1]], Bin)
@variable(m, 0 <= real_harvest_year_lambda[properties,real_harvests,years_real] <= 1)

# harvest year interpolated from lambdas
@constraint(m, [p in properties, r in real_harvests],
            real_harvest_year[p,r] == sum(real_harvest_year_lambda[p,r,t] * (t+n_years_int)
                                          for t in years_real))

# Lambdas must sum to 1 to be valid interpolation, or 0 if not harvested
@constraint(m, [p in properties, r in real_harvests],
            sum(real_harvest_year_lambda[p,r,t] for t in years_real) == real_harvest[p,r])

# Only choose 1 pair of consecutive lambdas using bucket variable to signal which pair
@constraint(m, [p in properties, r in real_harvests],
            sum(real_harvest_year_bucket[p,r,t] for t in years_real[1:end-1]) == real_harvest[p,r])
# Left endpoint
@constraint(m, [p in properties, r in real_harvests],
            real_harvest_year_lambda[p,r,1] <= real_harvest_year_bucket[p,r,1])
# Midpoints
@constraint(m, [p in properties, r in real_harvests, t in years_real[1:end-2]],
            real_harvest_year_lambda[p,r,t+1] <= real_harvest_year_bucket[p,r,t] + real_harvest_year_bucket[p,r,t+1])
# Right endpoint
@constraint(m, [p in properties, r in real_harvests],
            real_harvest_year_lambda[p,r,n_years_real] <= real_harvest_year_bucket[p,r,n_years_real-1])

# Need a variable to hold the yield harvested from a property if it was
# harvested in a particular year bucket, so it can be summed for enforcing
# constraints on yield per whole year
@variable(m, real_yield_if_harvested_on_year[properties,real_harvests,years_real])

@constraint(m, [p in properties, r in real_harvests, t in years_real[1:end-1]],
            real_yield_if_harvested_on_year[p,r,t] <= M * real_harvest_year_bucket[p,r,t])

@constraint(m, [p in properties, r in real_harvests, t in years_real[1:end-1]],
            real_yield_if_harvested_on_year[p,r,t] <= sum(lambda[p,r+n_years_int,pa] * yields[p,pa]
                                                          for pa in property_age))

@constraint(m, [p in properties, r in real_harvests, t in years_real[1:end-1]],
            real_yield_if_harvested_on_year[p,r,t] >= sum(lambda[p,r+n_years_int,pa] * yields[p,pa]
                                                          for pa in property_age)
                                                      - M * (1 - real_harvest_year_bucket[p,r,t]))

##################################################
# (Int) age constraints

# Assign initial ages
@constraint(m, [p in properties],
            age[p,1] == initial_ages[p] + 1)
# 'age' must be 1 for the year after a harvest. For any other year, 'age' must
# be the previous age +1.
@constraint(m, [p in properties, t in years_int],
            age[p, t] >= 1)
@constraint(m, [p in properties, prev_t in years_int, next_t in years_int; prev_t+1 == next_t],
            age[p, next_t] <= age[p, prev_t] + 1)

@constraint(m, [p in properties, prev_t in years_int, next_t in years_int; prev_t+1 == next_t],
            age[p, next_t] <= 1 + M * (1 - harvest[p,prev_t]))

@constraint(m, [p in properties, prev_t in years_int, next_t in years_int; prev_t+1 == next_t],
            age[p, next_t] >= age[p, prev_t] + 1 - M * harvest[p,prev_t])

# harvest_age is zero for years where there is no harvest
@constraint(m, [p in properties, t in years_int],
            harvest_age[p,t] <= M * harvest[p,t])
# harvest_age is the age of the property when it is harvested.
# To maximise the objective, the solver will automatically push
# harvest_age up to the age of the property on harvest years.
@constraint(m, [p in properties, t in years_int],
            harvest_age[p,t] <= age[p,t])
# Force 'harvest_age[p,t]' to be 'age[p,t]' when harvest[p,t]
@constraint(m, [p in properties, t in years_int],
            harvest_age[p,t] >= age[p,t] - M * (1 - harvest[p,t]))


##################################################
# Must provide a minimum quantity per year.
@constraint(m, [t in years_int],
            sum(lambda[p,t,pa] * yields[p,pa]
                for p in properties,
                pa in property_age)
            >= minimum_yield_per_year)
# Real
@constraint(m, [t in years_real],
            sum(real_yield_if_harvested_on_year[p,r,t]
                for p in properties,
                r in real_harvests)
            >= minimum_yield_per_year)

# The amount of product that can be processed is limited.
@constraint(m, [t in years_int],
            sum(lambda[p,t,pa] * yields[p,pa]
                for p in properties,
                pa in property_age)
            <= maximum_yield_per_year)
# Real
@constraint(m, [t in years_real],
            sum(real_yield_if_harvested_on_year[p,r,t]
                for p in properties,
                r in real_harvests)
            <= maximum_yield_per_year)


# Maximise profit: balance yield against replanting cost
@objective(m, Max,
           sum(lambda[p,t,pa] * yields[p,pa] * netReturns[p]
               for p in properties,
               t in years_int,
               pa in property_age)
           + sum(lambda[p,r+n_years_int,pa] * yields[p,pa] * netReturns[p]
                 for p in properties,
                 r in real_harvests,
                 pa in property_age)
           - sum(harvest[p,t] * replanting_cost * area[p]
                 for p in properties,
                 t in years_int)
           - sum(real_harvest[p,r] * replanting_cost * area[p]
                 for p in properties,
                 r in real_harvests))

# Integer variables >= 0
@constraint(m, [p in properties, t in years_int],
            age[p,t] >= 0)
@constraint(m, [p in properties, t in years_int],
            harvest_age[p,t] >= 0)
@constraint(m, [p in properties, r in real_harvests],
            real_harvest_age[p,r] >= 0)
@constraint(m, [p in properties, r in real_harvests, t in years_real],
            real_yield_if_harvested_on_year[p,r,t] >= 0)


solve(m)

println("harvest_ages[p,t]:")
print_harvests(harvest_age)
println("yields[p,t]:")
print_yields(yields)
println("yields per year:")
for t in years_int
    print("$(sum(getvalue(lambda[p,t,pa]*yields[p,pa]) for p in properties, pa in property_age)) ")
end
println()

println("real_harvest:")
for p in properties
    for r in real_harvests
        print("$(Int(round(getvalue(real_harvest[p,r])))) ")
    end
    println()
end

println("real_harvest_year:")
for p in properties
    for r in real_harvests
        print("$(getvalue(real_harvest_year[p,r])) ")
    end
    println()
end

println("real_harvest_age:")
for p in properties
    for r in real_harvests
        print("$(getvalue(real_harvest_age[p,r])) ")
    end
    println()
end

println("real_harvest_year_lambda:")
for t in years_real
    for r in real_harvests
        print("$(getvalue(real_harvest_year_lambda[n_properties,r,t])) ")
    end
    println()
end

println("real_harvest_year_bucket:")
for t in years_real[1:end-1]
    for r in real_harvests
        print("$(Int(round(getvalue(real_harvest_year_bucket[n_properties,r,t])))) ")
    end
    println()
end

println("real_yield_if_harvested_on_year:")
#for p in properties
p = n_properties
    for t in years_real
        for r in real_harvests
            print("$(getvalue(real_yield_if_harvested_on_year[p,r,t])) ")
        end
        println()
    end
#    println()
#    println()
#end


println("Total cost $(getobjectivevalue(m))")
println("Bound is $(getobjectivebound(m))")
