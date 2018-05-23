using JuMP
using Cbc
using Clp

firstPropertyRow = 5
ageColumn = 3
areaColumn = 4
net_return_column = 5
firstDensityColumn = 9
lastDensityColumn = 30
nProperties = 199
inputfile = "YieldsV3.csv"
tie_break_slope = 0.001
minGrowth = 1

function readData(filename)

    data = readcsv(filename)

    properties = 1:nProperties
    densityYears = 1:(lastDensityColumn-firstDensityColumn+1)

    netReturns = data[firstPropertyRow:firstPropertyRow+nProperties-1,net_return_column]
    initalAge = data[firstPropertyRow:firstPropertyRow+nProperties-1,ageColumn]
    area = data[firstPropertyRow:firstPropertyRow+nProperties-1,areaColumn]
    density = data[firstPropertyRow:firstPropertyRow+nProperties-1,firstDensityColumn:lastDensityColumn]
    #Times by age to get Yeild Loop up table

    yields = Array{Float64}(nProperties,lastDensityColumn-firstDensityColumn+1)

    for i in properties
        for j in densityYears
            yields[i,j] = density[i,j]*area[i]
            # if j>1
            #     if yields[i,j] - yields[i,j-1] <= minGrowth
            #         yields[i,j] = yields[i,j] - tie_break_slope*yields[i,j]
            #     end
            # end
        end
    end

    growthRate =  zeros(length(properties),1)
    for i in properties
        growthRate[i] = (density[i,lastDensityColumn-firstDensityColumn]-density[i,1])/
                                        (lastDensityColumn-firstDensityColumn)
    end

    #returns
    initalAge, area, density, growthRate, yields, netReturns
end

initialAge, area, density, growthRate, yields, netReturns = readData(inputfile)
