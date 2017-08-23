.onAttach <- function(lib, pkg) {
    packageStartupMessage("DMCHMM package, Version ",
                                utils::packageVersion("DMCHMM"),
        ", Released ", utils::packageDescription("DMCHMM")$Date,
        "\n", utils::packageDescription("DMCHMM")$Description,
        "\nBugReports: ", utils::packageDescription("DMCHMM")$BugReports
    )
}

