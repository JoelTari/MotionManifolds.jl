using Documenter
using ManifoldExtras

# DocMeta.setdocmeta!(ManifoldExtras, :DocTestSetup, :(using ManifoldExtras); recursive=true)
makedocs(
    sitename = "ManifoldExtras",
    format = Documenter.HTML(),
    modules = [ManifoldExtras]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
