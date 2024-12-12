import DataStructures: OrderedDict
using Documenter
using Literate
using HolisticElectricityModel

if haskey(ENV, "GITHUB_ACTIONS")
    ENV["JULIA_DEBUG"] = "Documenter"
end

pages = OrderedDict(
        "Welcome Page" => "index.md",
        "Quick Start Guide" => "quickstartguide.md",
        "Developer Guide" => Any[
            "dev_guide/logging.md",
            "dev_guide/inputs.md"
        ],
        "API" => "HolisticElectricityModel.md"
)

# This code performs the automated addition of Literate - Generated Markdowns. The desired
# section name should be the name of the file for instance foo_bar.jl -> Foo Bar
julia_file_filter = x -> occursin(".jl", x)
folders = Dict(
    "Developer Guide" => filter(julia_file_filter, readdir("docs/src/dev_guide")),
)

for (section, folder) in folders
    for file in folder
        section_folder_name = lowercase(replace(section, " " => "_"))
        outputdir = joinpath(pwd(), "docs", "src", "$section_folder_name")
        inputfile = joinpath("$section_folder_name", "$file")
        outputfile = string("generated_", replace("$file", ".jl" => ""))
        Literate.markdown(joinpath(pwd(), "docs", "src", inputfile),
                          outputdir;
                          name = outputfile,
                          credit = false,
                          execute=true)
        subsection = titlecase(replace(split(file, ".")[1], "_" => " "))
        push!(pages[section], ("$subsection" =>  joinpath("$section_folder_name", "$(outputfile).md")))
    end
end

makedocs(
    modules = [HolisticElectricityModel],
    format = Documenter.HTML(prettyurls = haskey(ENV, "GITHUB_ACTIONS"),),
    sitename = "HolisticElectricityModel.jl",
    authors = "Elaine Hale, Nongchao Guo, Burcin Cakir Erdener, Haley Ross, Ellie Estreich, and Purboday Ghosh",
    pages = Any[p for p in pages]
)

deploydocs(
    repo = "https://github.nrel.gov/HEM/HolisticElectricityModel.jl",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#"],
    push_preview = true,
)
