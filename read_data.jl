using JuMP
using Cbc
using Clp

firstPropertyRow = 5
ageColumn = 3
areaColumn = 4
firstDensityColumn = 8
lastDensityColumn = 29
nProperties = 10
inputfile = "YieldsV2.csv"
tie_break_slope = 0.001
minGrowth = 1

function readData(filename)

    data = readcsv(filename)

    properties = 1:nProperties
    densityYears = 1:lastDensityColumn-firstDensityColumn

    initalAge = data[firstPropertyRow:firstPropertyRow+nProperties-1,ageColumn]
    area = data[firstPropertyRow:firstPropertyRow+nProperties-1,areaColumn]
    density = data[firstPropertyRow:firstPropertyRow+nProperties-1,firstDensityColumn:lastDensityColumn]
    #Times by age to get Yeild Loop up table

    yields = Array{Float64}(nProperties,lastDensityColumn-firstDensityColumn)

    for i in properties
        for j in densityYears
            yields[i,j] = density[i,j]*area[i]
            if j>1
                if yields[i,j] - yields[i,j-1] <= minGrowth
                    yields[i,j] = yields[i,j] - tie_break_slope*yields[i,j]
                end
            end
        end
    end

    growthRate =  Dict{Int, Float64}()
    for i in properties
        growthRate[i] = (density[i,lastDensityColumn-firstDensityColumn]-density[i,1])/
                                        (lastDensityColumn-firstDensityColumn)
    end

    #returns
    initalAge, area, density, growthRate, yields
end

initalAge, area, density, growthRate, yields = readData(inputfile)
