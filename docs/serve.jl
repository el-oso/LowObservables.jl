# Local preview server for the built docs.
#
#   julia --project=docs docs/serve.jl [port]   # default port 8000
#
# Unlike `LiveServer.serve`, this maps VitePress "clean URLs" to the on-disk
# files, so `http://localhost:8000/physical-optics` works (no `.html` needed).
# The built site uses root-relative asset paths, so it MUST be served with the
# build output as the web root — which is exactly what this does.

using HTTP, Sockets

const ROOT = let b1 = joinpath(@__DIR__, "build", "1")
    isfile(joinpath(b1, "index.html")) ? b1 : joinpath(@__DIR__, "build")
end

const MIMES = Dict(
    ".html"=>"text/html", ".js"=>"text/javascript", ".mjs"=>"text/javascript",
    ".css"=>"text/css", ".json"=>"application/json", ".map"=>"application/json",
    ".woff2"=>"font/woff2", ".woff"=>"font/woff", ".svg"=>"image/svg+xml",
    ".png"=>"image/png", ".jpg"=>"image/jpeg", ".ico"=>"image/x-icon", ".txt"=>"text/plain",
)

# Map a request path to a file on disk, trying clean-URL fallbacks.
function resolve(target)
    path = HTTP.URIs.unescapeuri(first(split(target, '?')))
    path == "/" && (path = "/index.html")
    cands = [path]
    if endswith(path, "/")
        push!(cands, path * "index.html")
    elseif isempty(splitext(path)[2])          # extensionless → VitePress clean URL
        push!(cands, path * ".html", path * "/index.html")
    end
    for c in cands
        f = normpath(joinpath(ROOT, lstrip(c, '/')))
        startswith(f, ROOT) || continue        # block path traversal
        isfile(f) && return f
    end
    return nothing
end

function handler(req::HTTP.Request)
    f = resolve(req.target)
    f === nothing && return HTTP.Response(404, "404 Not Found: $(req.target)")
    ct = get(MIMES, lowercase(splitext(f)[2]), "application/octet-stream")
    return HTTP.Response(200, ["Content-Type" => ct], read(f))
end

port = isempty(ARGS) ? 8000 : parse(Int, ARGS[1])
@info "Serving $ROOT at http://localhost:$port  (clean URLs work, e.g. /physical-optics)"
HTTP.serve(handler, Sockets.localhost, port)
