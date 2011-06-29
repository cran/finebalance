.First.lib =
function(lib, pkg)
{
    # Startup Mesage and Desription:
    MSG <- if(getRversion() >= "2.4.1") packageStartupMessage else message
    dsc <- packageDescription(pkg)
    if(interactive() || getOption("verbose")) {
        # not in test scripts
        MSG(sprintf("qlspack Package %s (%s) loaded.", pkg, dsc$Version))
    }

    # Load dll:
    # library.dynam("fPortfolio", pkg, lib)
    # use "Rdonlp2"
}

.onLoad <-
    function(libname, pkgname)
{

}
