library(ini)
library(jsonlite)

#' load_rs_config
#'
#' Load ReadStore configuration from file or ENV variables
#' 
#' @param filename Path to configuration file
#' @param default_endpoint_url Default endpoint URL to set if not found in file or ENV
#' @param default_fastq_extensions Default fastq extensions to set if not found in file or ENV
#' @param default_output Default output format to set if not found in file or ENV
#' @return List containing username, token, endpoint URL, fastq extensions and output format
load_rs_config <- function(
    filename = NULL,
    default_endpoint_url = NULL,
    default_fastq_extensions = NULL,
    default_output = NULL
) {
    
    # Initialize variables
    username <- NULL
    token <- NULL
    endpoint_url <- NULL
    fastq_extensions <- NULL
    output <- NULL
    
    # If filename is provided, check if file exists
    if (!is.null(filename) && file.exists(filename)) {
        rs_config <- ini::read.ini(file = filename)
        
        username <- rs_config$credentials$username
        token <- rs_config$credentials$token
        endpoint_url <- rs_config$general$endpoint_url
        fastq_extensions <- rs_config$general$fastq_extensions
        fastq_extensions <- gsub("'", "\"", fastq_extensions)
        fastq_extensions <- jsonlite::fromJSON(fastq_extensions)
        output <- rs_config$general$output
    }
    
    # Check if ENV variables are set
    if (!(Sys.getenv("READSTORE_USERNAME") == "")) {
        username <- Sys.getenv("READSTORE_USERNAME")
    }
    if (!(Sys.getenv("READSTORE_TOKEN") == "")) {
        token <- Sys.getenv("READSTORE_TOKEN")
    }
    if (!(Sys.getenv("READSTORE_ENDPOINT_URL") == "")) {
        endpoint_url <- Sys.getenv("READSTORE_ENDPOINT_URL")
    }
    if (!(Sys.getenv("READSTORE_FASTQ_EXTENSIONS") == "")) {
        fastq_extensions <- strsplit(Sys.getenv("READSTORE_FASTQ_EXTENSIONS"), ",")[[1]]
    }
    if (!(Sys.getenv("READSTORE_DEFAULT_OUTPUT") == "")) {
        output <- Sys.getenv("READSTORE_DEFAULT_OUTPUT")
    }
    
    # If config parameters are not found in file or ENV, try to use default arguments if provided
    if (is.null(username)) {
        stop("Config: Username Not Found")
    }
    if (is.null(token)) {
        stop("Config: Token Not Found")
    }
    
    if (is.null(endpoint_url)) {
        if (is.null(default_endpoint_url)) {
            stop("Config: Endpoint URL Not Found")
        }

        endpoint_url <- default_endpoint_url
    }
    
    if (is.null(fastq_extensions)) {
        if (is.null(default_fastq_extensions)) {
            stop("Config: Fastq Extensions Not Found")
        }
        
        fastq_extensions <- default_fastq_extensions
    }
    
    if (is.null(output)) {
        if (is.null(default_output)) {
            stop("Config: Output Not Found")
        }
        
        output <- default_output
    }
    
    if (!(grepl('/$', endpoint_url))) {
        endpoint_url <- paste0(endpoint_url, '/')
    }

    return(list(username = username, token = token, endpoint_url = endpoint_url, fastq_extensions = fastq_extensions, output = output))
}

#' write_rs_config
#'
#' Write ReadStore configuration to file in INI format
#' Set permissions to 0600 (read/write only by owner)
#' 
#' @param filename Path to configuration file
#' @param username Username
#' @param token Token
#' @param endpoint_url Endpoint URL
#' @param fastq_extensions Fastq extensions
#' @param output Output format
write_rs_config <- function(filename,
                            username,
                            token,
                            endpoint_url,
                            fastq_extensions,
                            output) {

    if (!dir.exists(dirname(filename))) {
        stop("Directory Not Found")
    }

    # Create a list to hold the configuration
    config <- list(
        general = list(
            endpoint_url = endpoint_url,
            fastq_extensions = paste(fastq_extensions, collapse = ","),
            output = output
        ),
        credentials = list(
            username = username,
            token = token
        )
    )

    # Write the configuration to a file
    ini::write.ini(config, filename)
    Sys.chmod(filename, mode = "0600")
}