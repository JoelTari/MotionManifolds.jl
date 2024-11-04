using Documenter, MotionManifolds, StaticArrays

DocMeta.setdocmeta!(MotionManifolds, :DocTestSetup, :(using StaticArrays, MotionManifolds); recursive=true)

makedocs(sitename = "MotionManifolds", modules = [MotionManifolds], doctest = true)
