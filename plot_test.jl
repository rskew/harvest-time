using JuMP
using Plots
using Plotly
using StatPlots

# Yeild Per Year
yields_per_year = Array{Float64}(years_to_plan)

for t in years_to_plan
    year_total = 0
    for p in properties
        A = Int(round(getvalue(harvest_age[p,t])))
        if A != 0
            year_total = year_total + yields[p,A]
        end
    end
    yields_per_year[t] = year_total
end
println("Yeilds Per Year")
println(yields_per_year)

display(Plots.bar(years_to_plan,yields_per_year))

prop_yields = Array{Int}(length(years_to_plan),length(properties))
# age_tally = Array{Int}(length(25))
for t in years_to_plan
    for p in properties
         A = Int(round(getvalue(harvest_age[p,t])))
        if A != 0
            prop_yields[t,p] = round(A*yields[p,A])
            # age_tally[A] = age_tally[A]+1
        else
            prop_yields[t,p] = 0
        end
    end
end
display(groupedbar(prop_yields, bar_position = :stack, bar_width=0.7,
    title="Total yields per year broken down into properties",
    xlabel="Years",ylabel="Tonnes",legend = false))
# display(bar(age_tally))

#Bug in re setting the age
#Big in the total yields


# Harvests Per Year
#

# x = [1 2 3 4 5 6]
# y = x*2
# p1 = plot(x,y,lw=5) # Make a line plot
# p2 = scatter(x,y) # Make a scatter plot
# p3 = plot(x,y,xlabel="This one is labelled",lw=3,title="Subtitle")
# p4 = histogram(x,y) # Four histograms each with 10 points? Why not!
# display(plot(p1,p2,p3,p4,layout=(2,2),legend=false,show=true))
