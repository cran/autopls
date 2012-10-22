.onLoad <- function(lib, pkg)
{
	pkg.info <- utils::packageDescription('autopls')
	packageStartupMessage(paste("autopls", pkg.info[["Version"]]))
}
