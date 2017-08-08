.onAttach <- function(lib, pkg){
    pkgversion <- utils::packageVersion("DMCHMM")
    pkgdate <- utils::packageDescription("DMCHMM")$Date
    pkgdescrip <- utils::packageDescription("DMCHMM")$Description
    pkgBugReports <- utils::packageDescription("DMCHMM")$BugReports
    packageStartupMessage(
        paste('DMCHMM package, Version ',pkgversion,', Released ',pkgdate,'\n',
    pkgdescrip, '\nBugReports: ',pkgBugReports , sep = "")
    )
}

