using Test, Documenter, Hubbard

DocMeta.setdocmeta!(Hubbard, :DocTestSetup, :(using Hubbard); recursive=true)
doctest(Hubbard)
