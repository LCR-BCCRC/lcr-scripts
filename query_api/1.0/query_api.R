# USAGE:
# source("query_api.R")
# Requires that your GSC JIRA/wiki password is stored in an environment variable
# See https://www.bcgsc.ca/wiki/display/PIPE/API+interface+module+How-to for instructions

# Find merged bam file paths for some libraries
# bams <- query_api(endpoint = "merge", library = c("A68797", "A55760"))

# Find fastq files for some libraries
# fastqs <- query_api(endpoint = "merge_fastq", library = "F138635")

# Get comprehensive details about individual lanes of sequencing
# libcore <- query_api(endpoint = "aligned_libcore/info", library = c("A68797", "A55760"))
# libcore$bio_qc_comments

# For more details about different endpoints see https://www.bcgsc.ca/data/bioapps-docs/bioapps.api/docs/.

suppressWarnings(suppressPackageStartupMessages({
    message("Loading packages...")

    required_packages <- c(
        # Included to avoid mysterious error
        "xml2",
        # Password prompt
        "getPass",
        # Data handling
        "readr",
        "dplyr",
        "purrr",
        "tidyr",
        "lubridate",
        # Web requests
        "httr",
        "jsonlite"
    )

    installed_packages <- rownames(installed.packages())

    for (pkg in required_packages) {
        if (!pkg %in% installed_packages) {
            message(paste0("Installing ", pkg, "..."))
            install.packages(pkg, repos = "https://mirror.rcg.sfu.ca/mirror/CRAN/")
        }
        library(pkg, character.only = TRUE)
    }
}))

# Setup session -----------------------------------------------------------

message("Setting up session...")

base_url <- handle("https://sbs.bcgsc.ca:8100/")

if (Sys.getenv("API_PWD") != "") {
    message("Using API_PWD environment variable to login. ")
} else {
    stop("Please set your API_PWD following the instructions at https://www.bcgsc.ca/wiki/display/PIPE/API+interface+module+How-to")
}


session_response <- POST(
    handle = base_url,
    path = "session",
    encode = "json",
    body = list(
        username = Sys.getenv("USER"),
        password = Sys.getenv("API_PWD")
    )
)
session_content <- content(session_response)

if (status_code(session_response) != 200) {
    error <- pluck(session_content, "errors", 1, "description")

    if (is.null(error)) {
        print(session_content)
        stop(paste0("Failed to authenticate with the above error."))
    } else {
        stop(paste0("Failed to authenticate with the following error: \n", error))
    }
}

session_headers <- add_headers(
    "X-Token" = session_content$token,
    "Content-Type" = "application/json",
    "Accept" = "application/json"
)

query_api <- function(endpoint, ..., session = session_headers, url = base_url) {
    params <- list(...)

    # If using multiple queries, ensure that a nested list can be generated
    if (length(params) > 1) {
        stopifnot(
            "The second query must be a single element or the same length as the first query." = (length(params[[2]]) == 1 | length(params[[2]] == length(params[[1]])))
        )
    }

    query <- lapply(params, function(x) paste0(x, collapse = ","))

    # Query the API with the query values
    response <- GET(
        handle = url,
        path = endpoint,
        query = query,
        config = session
    )
    contents <- content(response)

    # Check if response is empty
    is_response_empty <- length(contents) == 0

    # Convert the response to a dataframe
    if (!is_response_empty) {
        contents <-
            response %>%
            content(as = "text", encoding = "UTF-8") %>%
            fromJSON() %>%
            jsonlite::flatten()

        contents
    } else {
        # If the response is empty, print the empty response
        message("The query produced an empty response. ")
        response
    }
}
