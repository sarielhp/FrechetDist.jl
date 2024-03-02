#! /bin/julia

tryusing(pkgsym) = try
    @eval using $pkgsym
    return true
catch e
    return e
end

if tryusing(:UrlDownload) != true
    using  Pkg
    Pkg.add("UrlDownload")
end


