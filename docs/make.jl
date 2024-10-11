using Documenter
using MotionManifolds

# DocMeta.setdocmeta!(MotionManifolds, :DocTestSetup, :(using MotionManifolds); recursive=true)
makedocs(
    sitename = "MotionManifolds",
    format = Documenter.HTML(),
    modules = [MotionManifolds]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
