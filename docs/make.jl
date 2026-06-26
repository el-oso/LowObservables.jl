using Documenter
using DocumenterVitepress
using LowObservables
using WGLMakie, Bonito   # interactive figures inlined into static HTML (offline export)
WGLMakie.activate!()

makedocs(;
    modules = [LowObservables],
    sitename = "LowObservables.jl",
    authors = "el_oso",
    format = DocumenterVitepress.MarkdownVitepress(
        devbranch = "master",
        devurl = "dev",
        repo = "github.com/el-oso/LowObservables.jl",
        sidebar_drawer = true,
    ),
    pages = [
        "Home" => "index.md",
        "What is RCS?" => "what-is-rcs.md",
        "Scattering Regimes & BOR" => "scattering-regimes.md",
        "Meshes" => "meshes.md",
        "Physical Optics" => "physical-optics.md",
        "Edge Waves (PTD)" => "edge-waves.md",
        "Radar Equation & Stealth" => "radar-equation.md",
        "API Reference" => "api.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs],
    remotes = nothing,
    doctest = false,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/el-oso/LowObservables.jl",
    devbranch = "master",
    push_preview = true,
)
