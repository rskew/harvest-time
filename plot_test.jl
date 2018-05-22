using JuMP
using Plots
using Plotly
using StatPlots

# Yeild Per Year
yields_per_year = Array{Float64}(years_to_plan)
prop_yields = Array{Int}(length(years_to_plan),length(properties))
age_tally = zeros(length(yields[1,:]))

for t in years_to_plan
    year_total = 0
    for p in properties
        A = Int(round(getvalue(harvest_age[p,t])))
        if A != 0
            prop_yields[t,p] = round(A*yields[p,A])
            year_total = year_total + A*yields[p,A]
            age_tally[A] = age_tally[A]+1
        else
            prop_yields[t,p] = 0
        end
    end
    yields_per_year[t] = year_total
end

display(Plots.bar(years_to_plan,yields_per_year,
    title="Total Yield Per Year",xlabel="Year in Plan",ylabel="Tonnes"))
display(Plots.bar(age_tally,
    title="Ages of Properties at Harvest",xlabel="Age",ylabel="Number of Properties"))
display(groupedbar(prop_yields, bar_position = :stack, bar_width=0.7,
    title="Total yields per year broken down into properties",
    xlabel="Years in Plan",ylabel="Tonnes"),legend=)


# Harvests Per Year
#

# x = [1 2 3 4 5 6]
# y = x*2
# p1 = plot(x,y,lw=5) # Make a line plot
# p2 = scatter(x,y) # Make a scatter plot
# p3 = plot(x,y,xlabel="This one is labelled",lw=3,title="Subtitle")
# p4 = histogram(x,y) # Four histograms each with 10 points? Why not!
# display(plot(p1,p2,p3,p4,layout=(2,2),legend=false,show=true))
