using PrecipSim
using Documenter

makedocs(;
    modules=[PrecipSim],
    authors="Patrick Eastham <peastham@math.fsu.edu> and contributors",
    repo="https://github.com/pseastham/PrecipSim.jl/blob/{commit}{path}#L{line}",
    sitename="PrecipSim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pseastham.github.io/PrecipSim.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pseastham/PrecipSim.jl",
)
