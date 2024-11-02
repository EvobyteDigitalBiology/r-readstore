library(ini)
library(assertthat)


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
    assert_that(!is.null(username), msg = "Username Not Found")
    assert_that(!is.null(token), msg = "Token Not Found")
    
    if (is.null(endpoint_url)) {
        assert_that(!is.null(default_endpoint_url), msg = "Config: Endpoint URL Not Found")
        endpoint_url <- default_endpoint_url
    }
    
    if (is.null(fastq_extensions)) {
        assert_that(!is.null(default_fastq_extensions), msg = "Config: Fastq Extensions Not Found")
        fastq_extensions <- default_fastq_extensions
    }
    
    if (is.null(output)) {
        assert_that(!is.null(default_output), msg = "Config: Output Format Not Found")
        output <- default_output
    }
    
    return(list(username = username, token = token, endpoint_url = endpoint_url, fastq_extensions = fastq_extensions, output = output))
}

write_rs_config <- function(filename, username, token, endpoint_url, fastq_extensions, output) {
            # Write configuration file
            #
            # Write configuration file containing
            # username, token, endpoint URL, fastq extensions and output format.
            #
            # Args:
            #     filename: Path to write configuration file
            #     username: Username
            #     token: Token
            #     endpoint_url: URL of the ReadStore API
            #     fastq_extensions (List[str]): List of fastq extensions
            #     output: Default output format

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