using JuMP
using Plots
using StatPlots

# Yeild Per Year
yields_per_year = Array{Float64}(years_to_plan)
prop_yields = zeros(length(years_to_plan),length(properties))
age_tally = zeros(length(yields[1,:]))
year_profit = zeros(years_to_plan)

for t in years_to_plan
    year_total = 0
    for p in properties
        A = Int(round(harvest_ages[p,t]))
        if A != 0
            prop_yields[t,p] = round(yields[p,A])
            year_total = year_total + yields[p,A]
            age_tally[A] = age_tally[A]+1
            tempProfit = yields[p,A] * netReturns[p]
            tempCost = replanting_cost * area[p]
            year_profit[t] = round(year_profit[t] + tempProfit - tempCost)
        end
    end
    yields_per_year[t] = year_total
end

#PROFIT PER YEAR
display(Plots.bar(year_profit,
    title="Total Profit/Year",xlabel="Year in Plan",ylabel="Profit",
    legend = false))
# Plots.savefig("C:\Users\Beau\Google Drive\aaaUNI\Optimisation for Industry\Figures\TotalSProfitYear.png")

#YIELD PER YEAR
display(Plots.bar(years_to_plan,yields_per_year,
        title="Total Yield/Year",xlabel="Year in Plan",ylabel="Tonnes",
        legend = false))

#HARVEST AGES
display(Plots.bar(age_tally,
    title="Ages of Properties at Harvest",xlabel="Age",ylabel="Number of Properties",
    legend = false))

#YIELD PER YEAR STACKED
display(groupedbar(prop_yields, bar_position = :stack, bar_width=0.7,
    title="Total Yield/Year/Property",  xlabel="Years in Plan",ylabel="Tonnes",
    legend = false))


# Harvests Per Year
#

# x = [1 2 3 4 5 6]
# y = x*2
# p1 = plot(x,y,lw=5) # Make a line plot
# p2 = scatter(x,y) # Make a scatter plot
# p3 = plot(x,y,xlabel="This one is labelled",lw=3,title="Subtitle")
# p4 = histogram(x,y) # Four histograms each with 10 points? Why not!
# display(plot(p1,p2,p3,p4,layout=(2,2),legend=false,show=true))
