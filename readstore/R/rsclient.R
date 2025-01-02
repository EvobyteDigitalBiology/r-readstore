library(httr)
library(base64enc)
library(jsonlite)

REST_API_VERSION = "api_x_v1/"
USER_AUTH_TOKEN_ENDPOINT = "auth_token/"
FASTQ_UPLOAD_ENDPOINT = "fq_file_upload/"
FQ_DATASET_ENDPOINT = "fq_dataset/"
FQ_FILE_ENDPOINT = "fq_file/"
FQ_ATTACHMENT_ENDPOINT = "fq_attachment/"
PROJECT_ENDPOINT = "project/"
PROJECT_ATTACHMENT_ENDPOINT = "project_attachment/"
PRO_DATA_ENDPOINT = "pro_data/"

METADATA_RESERVED_KEYS = c('id',
                        'name',
                        'project',
                        'project_ids',
                        'project_names',
                        'owner_group_name',
                        'qc_passed',
                        'paired_end',
                        'index_read',
                        'created',
                        'description',
                        'owner_username',
                        'fq_file_r1',
                        'fq_file_r2',
                        'fq_file_i1',
                        'fq_file_i2',
                        'id_project',
                        'name_project',
                        'name_og',
                        'archived',
                        'collaborators',
                        'dataset_metadata_keys',
                        'data_type',
                        'version',
                        'valid_to',
                        'upload_path',
                        'owner_username',
                        'fq_dataset',
                        'id_fq_dataset',
                        'name_fq_dataset')

#' validate_charset
#'
#' Check if the query string contains only allowed characters
#' Allowed characters are: 0-9, a-z, A-Z, _, -, ., @
#' 
#' @param query_str String to validate 
#' @return TURE if valid charset, FALSE otherwise
validate_charset <- function(query_str) {
  # Validate charset for query string
  
  allowed <- c(0:9, letters, LETTERS, '_', '-', '.', '@')
  allowed <- unique(allowed)
  
  return(all(strsplit(query_str, NULL)[[1]] %in% allowed))
}

#' validate_metadata
#'
#' Check if metadata keys are valid
#' Allowed characters are: 0-9, a-z, A-Z, _, -, ., @
#' Must not be empty
#' Must not be a reserved keyword
#' 
#' @param metadata List if metadata entries
validate_metadata <- function(metadata) {
  # Validate metadata list
  
  # Ensure keys are non-empty, valid charset, and not reserved
  
  for (key in names(metadata)) {
    if (key == "") {
      stop("Empty Key")
    }
    
    if (!validate_charset(key)) {
      stop("Invalid character in key. Must be alphanumeric or _-.@")
    }
    
    if (key %in% METADATA_RESERVED_KEYS) {
      stop(paste("Reserved Keyword not allowed in metadata:", key))
    }
  }
}

#' test_server_connection
#'
#' Check if the server at the endpoint is reachable
#' 
#' @param endpoint_url Endpoint to check 
#' @return TURE if the server is reachable, FALSE otherwise
test_server_connection <- function(endpoint_url) {
    
    # Create a list of headers
    parsed_url <- httr::parse_url(endpoint_url)
    scheme <- parsed_url$scheme

    if (!(scheme %in% c("http", "https"))) {        
        return(FALSE)
    } else {
        tryCatch({
            res <- httr::HEAD(endpoint_url)

            if (res$status_code == 200) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        }, error = function(e) {
            return(FALSE)
        })
    }
}

#' auth_user_token
#'
#' Check if the user token is valid
#' 
#' @param username username
#' @param token token to authenticate
#' @param endpoint_url Endpoint to check
#' @return TRUE if the server token and user is valid, FALSE otherwise
auth_user_token <- function(username, token, endpoint_url) {

    tryCatch({
    # Construct the authentication endpoint URL
    auth_endpoint <- file.path(endpoint_url, USER_AUTH_TOKEN_ENDPOINT, fsep = "")
    
    # Send a POST request
    res <- httr::POST(auth_endpoint, httr::authenticate(username, token), encode = "json")
    
    # Check the status code
    if (httr::status_code(res) != 200) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }, error = function(e) {
    # Catch connection errors
    return(FALSE)
  })
}

#' get_rs_client
#'
#' Get a rs client object (client template)
#' 
#' @param username username
#' @param token token to authenticate
#' @param endpoint_url Endpoint to check
#' @return list with username, token and endpoint_url defining the ReadStore client
get_rs_client <- function(username, token, endpoint_url) {

    endpoint_url = file.path(endpoint_url, REST_API_VERSION, fsep = "")

    if (!test_server_connection(endpoint_url)) {
        stop("Connection Error")
    }
    if (!auth_user_token(username, token, endpoint_url)) {
        stop("Authentication Error")
    }

    return(list(username = username, token = token, endpoint_url = endpoint_url))
}

#' upload_fastq_rs
#'
#' Upload fastq files to the ReadStore API
#' 
#' @param client ReadStore client
#' @param fastq_file string or vector of fastq files to upload
#' @param fastq_name optional string or vector of fastq names to upload
#' @param read_type optional string or vector of read types to upload
upload_fastq_rs <- function(client,
                            fastq_file,
                            fastq_name = NULL,
                            read_type = NULL) {

    fq_upload_endpoint = file.path(client$endpoint_url, FASTQ_UPLOAD_ENDPOINT, fsep = "")
        
    fq_path = normalizePath(fastq_file)

    json_payload = list(fq_file_path = fq_path)

    if (!is.null(fastq_name)) {
        if (fastq_name == "") {
            stop("fastq_name cannot be empty")
        } 
        if (!validate_charset(fastq_name)) {
            stop("fastq_name contains invalid characters")
        }
        json_payload$fq_file_name = fastq_name
    }

    if (!is.null(read_type)) {
        if (!(read_type %in% c("R1","R2","I1","I2"))) {
            stop("Invalid Read Type")
        }
        json_payload$read_type = read_type
    }

    if (!(file.exists(fq_path))) {
        stop(paste("fastq file not found for path", fq_path))
    }
    if (file.access(fq_path, mode=4) != 0) {
        stop(paste("No read permissions for path", fq_path))
    }
    
    res <- httr::POST(fq_upload_endpoint,
                        body = json_payload,
                        config = httr::authenticate(client$username, client$token),
                        encode = "json")

    if (!(httr::status_code(res) %in% c(200, 204))) {
        stop("upload_fastq failed")
    }
}

#' get_fq_file_rs
#'
#' Get a fastq file from the ReadStore API
#' 
#' @param client ReadStore client
#' @param fq_file_id id of the fastq file
#' @return json object (list) with fastq file
get_fq_file_rs <- function(client, fq_file_id) {

    fq_file_endpoint = file.path(client$endpoint_url, FQ_FILE_ENDPOINT, fsep = "")
    fq_file_endpoint = file.path(fq_file_endpoint, fq_file_id, '/', fsep = "")

    res <- httr::GET(fq_file_endpoint,httr::authenticate(client$username, client$token), encode = "json")

    if (!(httr::status_code(res) %in% c(200, 204))) {
        stop("get_fq_file Failed")
    } else {
        json = httr::content(res, "parsed")
        if (length(json) == 0) {
            return(list())
        } else if (length(json) == 1) {
           return(json[[1]])
        } else {
            stop("get_fq_file: Multiple fq files found for fq file id")
        }
    }
}

#' get_fq_file_upload_path_rs
#'
#' Get upload path for a fastq file from the ReadStore API
#' 
#' @param client ReadStore client
#' @param fq_file_id id of the fastq file
#' @return character with the upload path
get_fq_file_upload_path_rs <- function(client, fq_file_id) {
    # Get the upload path for a fastq file from the ReadStore API
    #
    # Args:
    #     client: ReadStore client
    #     fq_file_id: Fastq file ID
    #
    # Returns:
    #     Upload path

    fq_file = get_fq_file_rs(client, fq_file_id)

    if (!("upload_path" %in% names(fq_file))) {
        stop("Upload Path Not Found")
    }

    return(fq_file$upload_path)
}

#' list_fq_files_rs
#'
#' List fastq files from the ReadStore API
#' 
#' @param client ReadStore client
#' @return (json) list of projects (list objects)
list_fq_files_rs <- function(client) {

    fq_file_endpoint = file.path(client$endpoint_url, FQ_FILE_ENDPOINT, fsep = "")

    res <- httr::GET(fq_file_endpoint,httr::authenticate(client$username, client$token), encode = "json")

    if (!(httr::status_code(res) %in% c(200, 204))) {
        stop("get_fq_file Failed")
    } else {
        json = httr::content(res, "parsed")
        
        return(json)
    }
}

#' create_fq_file_rs
#'
#' Create fastq files from the ReadStore API
#' 
#' @param client ReadStore client
#' @param name Fastq file name
#' @param read_type Read type (R1, R2, I1, I2)
#' @param qc_passed QC Pass
#' @param read_length Read length
#' @param num_reads Number of reads
#' @param size_mb Size in MB
#' @param qc_phred_mean QC Phred Mean
#' @param qc_phred QC Phred (must be a named list in R)
#' @param upload_path Upload Path
#' @param md5_checksum MD5 Checksum
#' @param staging Staging
#' @param pipeline_version Pipeline Version
#' @return (json) list of projects (list objects)
create_fq_file_rs <- function(client,
                            name, 
                            read_type, 
                            qc_passed, 
                            read_length, 
                            num_reads, 
                            size_mb, 
                            qc_phred_mean, 
                            qc_phred, 
                            upload_path, 
                            md5_checksum, 
                            staging, 
                            pipeline_version) {
    # Create Fastq File in ReadStore

    fq_file_endpoint = file.path(client$endpoint_url, FQ_FILE_ENDPOINT, fsep = "")
    
    if (name == "") {
        stop("create_fq_file failed: Empty Name")
    }
    if (!(validate_charset(name))) {
        stop("create_fq_file failed: Invalid Name")
    }

    # Define JSON for POST request
    json_data <- list(
        name = name,
        bucket = 'hello',
        key = 'hello',
        read_type = read_type,
        qc_passed = qc_passed,
        read_length = read_length,
        num_reads = num_reads,
        size_mb = size_mb,
        qc_phred_mean = qc_phred_mean,
        qc_phred = qc_phred,
        upload_path = upload_path,
        md5_checksum = md5_checksum,
        staging = staging,
        pipeline_version = pipeline_version
    )

    # Make the POST request
    response <- httr::POST(
        url = fq_file_endpoint,
        body = json_data,
        encode = "json",
        config = httr::authenticate(client$username, client$token)
    )
    
    # Check the response
    if (httr::status_code(response) != 201) {
        tryCatch({
        detail <- httr::content(response, "parsed")
        }, error = function(e) {
        detail <- "No Message"
        })
        stop(paste("create_fq_file failed:", detail))
    } else {
        return(httr::content(response, "parsed"))
    }
}

#' update_fq_file_rs
#'
#' Update fastq files from the ReadStore API
#' 
#' @param client ReadStore client
#' @param fq_file_id Fastq file ID
#' @param name Fastq file name
#' @param read_type Read type (R1, R2, I1, I2)
#' @param qc_passed QC Pass
#' @param read_length Read length
#' @param num_reads Number of reads
#' @param size_mb Size in MB
#' @param qc_phred_mean QC Phred Mean
#' @param qc_phred QC Phred (must be a named list in R)
#' @param upload_path Upload Path
#' @param md5_checksum MD5 Checksum
#' @param staging Staging
#' @param pipeline_version Pipeline Version
#' @return (json) list of projects (list objects)
update_fq_file_rs <- function(client,
                            fq_file_id,
                            name, 
                            read_type, 
                            qc_passed, 
                            read_length, 
                            num_reads, 
                            size_mb, 
                            qc_phred_mean, 
                            qc_phred, 
                            upload_path, 
                            md5_checksum, 
                            staging, 
                            pipeline_version) {
    # Create Fastq File in ReadStore

    fq_file_endpoint <- file.path(client$endpoint_url, FQ_FILE_ENDPOINT, fsep = "")
    fq_file_endpoint <- file.path(fq_file_endpoint, fq_file_id, '/', fsep = "")

    if (name == "") {
        stop("update_fq_file failed: Empty Name")
    }
    if (!(validate_charset(name))) {
        stop("update_fq_file failed: Invalid Name")
    }

    # Define JSON for POST request
    json_data <- list(
        name = name,
        read_type = read_type,
        qc_passed = qc_passed,
        read_length = read_length,
        num_reads = num_reads,
        size_mb = size_mb,
        qc_phred_mean = qc_phred_mean,
        qc_phred = qc_phred,
        upload_path = upload_path,
        md5_checksum = md5_checksum,
        staging = staging,
        pipeline_version = pipeline_version
    )
    
    # Make the POST request
    response <- httr::PUT(
        url = fq_file_endpoint,
        body = json_data,
        encode = "json",
        config = httr::authenticate(client$username, client$token)
    )
    
    # Check the response
    if (httr::status_code(response) != 200) {
        tryCatch({
        detail <- httr::content(response, "parsed")
        }, error = function(e) {
        detail <- "No Message"
        })
        stop(paste("update_fq_file failed:", detail))
    } else {
        return(httr::content(response, "parsed"))
    }
}

#' delete_fq_file_rs
#'
#' Delete fastq file from the ReadStore API
#' 
#' @param client ReadStore client
#' @param fq_file_id Fastq file ID
#' @return (json) list of projects (list objects)
delete_fq_file_rs <- function(client, fq_file_id) {

    fq_file_endpoint = file.path(client$endpoint_url, FQ_FILE_ENDPOINT, fsep = "")
    fq_file_endpoint = file.path(fq_file_endpoint, fq_file_id, '/', fsep = "")

    res <- httr::DELETE(fq_file_endpoint,
                        httr::authenticate(client$username, client$token),
                        encode = "json")

    if ((httr::status_code(res) == 400) | (httr::status_code(res) == 404)) {
        json = httr::content(res, "parsed")
        detail = json$detail
        
        if (detail == "FqFile not found") {
            stop("FqFile not found")
        } else {
            stop(paste("delete_fq_file", detail))
        }
    } else if (httr::status_code(res) == 403) {
        json = httr::content(res, "parsed")
        detail = json$detail
        stop(paste("FqFile Delete failed", detail))
    } else if (httr::status_code(res) %in% c(200,204)) {
        json = httr::content(res, "parsed")
        fq_file_id = as.integer(json$id)
        return(fq_file_id)
    } else {
        stop("delete_fq_file failed")
    }
}

#' list_fastq_datasets_rs
#'
#' Get list of fastq datasets from the ReadStore API
#' 
#' @param client ReadStore client
#' @param project_id Project ID to filter datasets for
#' @param project_name Project name to filter datasets for
#' @return (json) list of fastq datasets
list_fastq_datasets_rs <- function(client, project_id = NULL, project_name = NULL) {

    fq_dataset_endpoint = file.path(client$endpoint_url, FQ_DATASET_ENDPOINT, fsep = "")

    json_payload = list()

    if (!is.null(project_name)) {
        json_payload$project_name = project_name
    }
    if (!is.null(project_id)) {
        json_payload$project_id = project_id
    }

    res <- httr::GET(fq_dataset_endpoint,
                    query = json_payload,
                    httr::authenticate(client$username, client$token),
                    encode = "json")

    if (!(httr::status_code(res) %in% c(200,204))) {
        stop("list_fastq_datasets Failed")
    }

    json = httr::content(res, "parsed")

    return(json)
}

#' get_fastq_dataset_rs
#'
#' Get fastq dataset from the ReadStore API
#' ID or Name must be provided
#' 
#' @param client ReadStore client
#' @param dataset_id Dataset ID to return
#' @param dataset_name Dataset name to return
#' @return json object (list) with fastq dataset
get_fastq_dataset_rs <- function(client, dataset_id = NULL, dataset_name = NULL) {

    fq_dataset_endpoint = file.path(client$endpoint_url, FQ_DATASET_ENDPOINT, fsep = "")

    if (is.null(dataset_id) & is.null(dataset_name)) {
        stop("dataset_id or dataset_name required")
    }

    json_payload = list()

    if (!is.null(dataset_id)) {
        json_payload$id = dataset_id
    }
    if (!is.null(dataset_name)) {
        json_payload$name = dataset_name
    }

    res <- httr::GET(fq_dataset_endpoint,
                    query = json_payload,
                    httr::authenticate(client$username, client$token),
                    encode = "json")
    
    if (!(httr::status_code(res) %in% c(200,204))) {
        stop("get_fastq_dataset Failed")
    } else {
        json = httr::content(res, "parsed")

        if (length(json) == 0) {
            return(list())
        } else if (length(json) == 1) {
           return(json[[1]])
        } else {
            stop("get_fastq_dataset: Multiple fq datasets found for dataset id or name")
        }
    }
}

#' create_fastq_dataset_rs
#'
#' Create fastq dataset from the ReadStore API
#' ID or Name must be provided
#' 
#' @param client ReadStore client
#' @param name Dataset name
#' @param description Dataset description
#' @param qc_passed QC Pass
#' @param paired_end Paired End
#' @param index_read Index Read
#' @param project_ids List of project IDs
#' @param project_names List of project names
#' @param metadata Metadata (must be a named list in R)
#' @param fq_file_r1_id Fastq file R1 ID
#' @param fq_file_r2_id Fastq file R2 ID
#' @param fq_file_i1_id Fastq file I1 ID
#' @param fq_file_i2_id Fastq file I2 ID
#' @return A list containing the created Fastq dataset
create_fastq_dataset_rs <- function(client,
                                    name, 
                                    description, 
                                    qc_passed, 
                                    paired_end, 
                                    index_read, 
                                    project_ids, 
                                    project_names, 
                                    metadata, 
                                    fq_file_r1_id = NULL, 
                                    fq_file_r2_id = NULL, 
                                    fq_file_i1_id = NULL, 
                                    fq_file_i2_id = NULL) {
  
    fq_dataset_endpoint = file.path(client$endpoint_url, FQ_DATASET_ENDPOINT, fsep = "")

    if (name == "") {
        stop("create_fastq_dataset failed: Empty Name")
    }
    if (!(validate_charset(name))) {
        stop("create_fastq_dataset failed: Invalid Name")
    }
    
    validate_metadata(metadata)

    if (length(metadata) == 0) {
        metadata <- setNames(as.list(rep("", length(metadata))), metadata)
    }    

    # Define JSON for POST request
    json_data <- list(
        name = name,
        description = description,
        qc_passed = qc_passed,
        paired_end = paired_end,
        index_read = index_read,
        project_ids = project_ids,
        project_names = project_names,
        metadata = metadata,
        fq_file_r1 = fq_file_r1_id,
        fq_file_r2 = fq_file_r2_id,
        fq_file_i1 = fq_file_i1_id,
        fq_file_i2 = fq_file_i2_id
    )
    
    json_data_str <- jsonlite::toJSON(json_data, null="null", auto_unbox = TRUE)

    # Make the POST request
    response <- httr::POST(
        url = fq_dataset_endpoint,
        body = json_data_str,
        encode = "raw",
        config = httr::authenticate(client$username, client$token),
        content_type("application/json")
    )

    # Check the response
    if (httr::status_code(response) != 201) {
        tryCatch({
        detail <- httr::content(response, "parsed")
        }, error = function(e) {
        detail <- "No Message"
        })
        stop(paste("create_fastq_dataset failed:", detail))
    } else {
        return(httr::content(response, "parsed"))
    }
}


#' update_fastq_dataset_rs
#'
#' Get update fastq dataset from the ReadStore API
#' ID or Name must be provided
#' 
#' @param client ReadStore client
#' @param dataset_id Dataset ID
#' @param name Dataset name
#' @param description Dataset description
#' @param qc_passed QC Pass
#' @param paired_end Paired End
#' @param index_read Index Read
#' @param project_ids List of project IDs
#' @param project_names List of project names
#' @param metadata Metadata (must be a named list in R)
#' @param fq_file_r1_id Fastq file R1 ID
#' @param fq_file_r2_id Fastq file R2 ID
#' @param fq_file_i1_id Fastq file I1 ID
#' @param fq_file_i2_id Fastq file I2 ID
#' @return A list containing the created Fastq dataset
update_fastq_dataset_rs <- function(client,
                                    dataset_id,
                                    name, 
                                    description, 
                                    qc_passed, 
                                    paired_end, 
                                    index_read, 
                                    project_ids, 
                                    project_names, 
                                    metadata, 
                                    fq_file_r1_id = NULL, 
                                    fq_file_r2_id = NULL, 
                                    fq_file_i1_id = NULL, 
                                    fq_file_i2_id = NULL) {
  
    fq_dataset_endpoint = file.path(client$endpoint_url, FQ_DATASET_ENDPOINT, fsep = "")
    fq_dataset_endpoint = file.path(fq_dataset_endpoint, dataset_id, '/', fsep = "")

    if (name == "") {
        stop("create_fastq_dataset failed: Empty Name")
    }
    if (!(validate_charset(name))) {
        stop("create_fastq_dataset failed: Invalid Name")
    }
    
    validate_metadata(metadata)
    
    # Here is an error
    if (length(metadata) == 0) {
        metadata <- setNames(as.list(rep("", length(metadata))), metadata)
    }

    # Define JSON for POST request
    json_data <- list(
        name = name,
        description = description,
        qc_passed = qc_passed,
        paired_end = paired_end,
        index_read = index_read,
        project_ids = project_ids,
        project_names = project_names,
        metadata = metadata,
        fq_file_r1 = fq_file_r1_id,
        fq_file_r2 = fq_file_r2_id,
        fq_file_i1 = fq_file_i1_id,
        fq_file_i2 = fq_file_i2_id
    )
    
    # TODO: Still brings in an error

    json_data_str <- jsonlite::toJSON(json_data, null="null", auto_unbox = TRUE)

    # Make the POST request
    response <- httr::PUT(
        url = fq_dataset_endpoint,
        body = json_data_str,
        encode = "raw",
        config = httr::authenticate(client$username, client$token),
        content_type("application/json")
    )

    # Check the response
    if (httr::status_code(response) != 200) {
        tryCatch({
        detail <- httr::content(response, "parsed")
        }, error = function(e) {
        detail <- "No Message"
        })
        stop(paste("update_fastq_dataset failed:", detail))
    } else {
        return(httr::content(response, "parsed"))
    }
}

#' delete_fq_file_rs
#'
#' Delete fastq file from the ReadStore API
#' 
#' @param client ReadStore client
#' @param dataset_id Fastq file ID
#' @return (json) list of projects (list objects)
delete_fastq_dataset_rs <- function(client, dataset_id) {

    fq_dataset_endpoint = file.path(client$endpoint_url, FQ_DATASET_ENDPOINT, fsep = "")
    fq_dataset_endpoint = file.path(fq_dataset_endpoint, dataset_id, '/', fsep = "")

    res <- httr::DELETE(fq_dataset_endpoint,
                        httr::authenticate(client$username, client$token),
                        encode = "json")

    if (httr::status_code(res) == 400) {
        json = httr::content(res, "parsed")
        detail = json$detail
        if (detail == "FqDataset not found") {
            stop("FqDataset not found")
        } else {
            stop(paste("delete_fastq_dataset", detail))
        }
    } else if (httr::status_code(res) == 403) {
        json = httr::content(res, "parsed")
        detail = json$detail
        stop(paste("FqDataset Delete failed", detail))
    } else if (httr::status_code(res) %in% c(200,204)) {
        json = httr::content(res, "parsed")
        dataset_id = as.integer(json$id)
        return(dataset_id)
    } else {
        stop("delete_fastq_dataset failed")
    }
}


#' list_projects_rs
#'
#' Get list of projects from the ReadStore API
#' 
#' @param client ReadStore client
#' @return (json) list of projects (list objects)
list_projects_rs <- function(client) {

    project_endpoint = file.path(client$endpoint_url, PROJECT_ENDPOINT, fsep = "")

    res <- httr::GET(project_endpoint,
                    httr::authenticate(client$username, client$token),
                    encode = "json")

    if (!(httr::status_code(res) %in% c(200,204))) {
        stop("list_projects failed")
    }

    json = httr::content(res, "parsed")

    return(json)
}

#' get_project_rs
#'
#' Get project from the ReadStore API
#' ID or Name must be provided
#' 
#' @param client ReadStore client
#' @param project_id Project ID to return
#' @param project_name Project name to return
#' @return json object (list) with project
get_project_rs <- function(client, project_id = NULL, project_name = NULL) {

    project_endpoint = file.path(client$endpoint_url, PROJECT_ENDPOINT, fsep = "")

    if (is.null(project_id) & is.null(project_name)) {
        stop("project_id or project_name required")
    }

    json_payload = list()

    if (!is.null(project_id)) {
        json_payload$id = project_id
    }
    if (!is.null(project_name)) {
        json_payload$name = project_name
    }

    res <- httr::GET(project_endpoint,
                    query = json_payload,
                    httr::authenticate(client$username, client$token),
                    encode = "json")

    if (!(httr::status_code(res) %in% c(200,204))) {
        stop("get_project failed")
    } else {
        json = httr::content(res, "parsed")
        if (length(json) == 0) {
            return(list())
        } else if (length(json) == 1) {
           return(json[[1]])
        } else {
            stop("get_project: Multiple projects found for projects id or name")
        }
    }
}


#' create_project_rs
#'
#' Create project from the ReadStore API
#' 
#' @param client ReadStore client
#' @param name Project ID to return
#' @param description Project name to return
#' @param metadata List of project metadata
#' @param dataset_metadata_keys vector of dataset metadata keys
#' @return json object (list) with project
create_project_rs <- function(client,
                            name, 
                            description, 
                            metadata, 
                            dataset_metadata_keys) {
    # Create Project in ReadStore

    project_endpoint = file.path(client$endpoint_url, PROJECT_ENDPOINT, fsep = "")
    
    # Validate name
    if (name == "") {
        stop("create_project failed: Empty Name")
    }
    if (!(validate_charset(name))) {
        stop("create_project failed: Invalid Name")
    }

    validate_metadata(metadata)
    validate_metadata(dataset_metadata_keys)

    if (length(metadata) == 0) {    
        metadata <- setNames(as.list(rep("", length(metadata))), metadata)
    }
    if (length(dataset_metadata_keys) == 0) {
        dataset_metadata_keys <- setNames(as.list(rep("", length(dataset_metadata_keys))), dataset_metadata_keys)
    }

    # Prepare JSON for POST request
    json_data <- list(
        name = name,
        description = description,
        metadata = metadata,
        dataset_metadata_keys = dataset_metadata_keys
    )
    
    json_data_str <- jsonlite::toJSON(json_data, auto_unbox = TRUE, null = "null")

    # Make the POST request
    response <- httr::POST(
        url = project_endpoint,
        body = json_data_str,
        encode = "raw",
        config = httr::authenticate(client$username, client$token),
        httr::add_headers(`Content-Type` = "application/json")
    )

    # Check the response
    if (httr::status_code(response) != 201) {
        tryCatch({
        detail <- httr::content(response, "parsed")
        }, error = function(e) {
        detail <- "No Message"
        })
        stop(paste("create_project failed:", detail))
    } else {
        return(httr::content(response, "parsed"))
    }
}


#' update_project_rs
#'
#' Update project from the ReadStore API
#' 
#' @param client ReadStore client
#' @param name Project ID to return
#' @param description Project name to return
#' @param metadata List of project metadata
#' @param dataset_metadata_keys vector of dataset metadata keys
#' @return json object (list) with project
update_project_rs <- function(client,
                              project_id,
                              name,
                              description, 
                              metadata,
                              dataset_metadata_keys) {
    # Create Project in ReadStore

    project_endpoint = file.path(client$endpoint_url, PROJECT_ENDPOINT, fsep = "")
    project_endpoint = file.path(project_endpoint, project_id, '/', fsep = "")
    
    # Validate name
    if (name == "") {
        stop("create_project failed: Empty Name")
    }
    if (!(validate_charset(name))) {
        stop("create_project failed: Invalid Name")
    }

    validate_metadata(metadata)
    validate_metadata(dataset_metadata_keys)

    if (length(metadata) == 0) {    
        metadata <- setNames(as.list(rep("", length(metadata))), metadata)
    }
    if (length(dataset_metadata_keys) == 0) {
        dataset_metadata_keys <- setNames(as.list(rep("", length(dataset_metadata_keys))), dataset_metadata_keys)
    }

    # Prepare JSON for POST request
    json_data <- list(
        name = name,
        description = description,
        metadata = metadata,
        dataset_metadata_keys = dataset_metadata_keys
    )

    json_data_str <- jsonlite::toJSON(json_data, auto_unbox = TRUE, null = "null")

    # Make the POST request
    response <- httr::PUT(
        url = project_endpoint,
        body = json_data_str,
        encode = "raw",
        config = httr::authenticate(client$username, client$token),
        httr::add_headers(`Content-Type` = "application/json")
    )

    # Check the response
    if (httr::status_code(response) != 200) {
        tryCatch({
        detail <- httr::content(response, "parsed")
        }, error = function(e) {
        detail <- "No Message"
        })
        stop(paste("create_project failed:", detail))
    } else {
        return(httr::content(response, "parsed"))
    }
}



#' delete_project_rs
#'
#' Delete Project from the ReadStore API
#' 
#' @param client ReadStore client
#' @param project_id Project ID
#' @return (json) list of projects (list objects)
delete_project_rs <- function(client, project_id) {

    project_endpoint = file.path(client$endpoint_url, PROJECT_ENDPOINT, fsep = "")
    project_endpoint = file.path(project_endpoint, project_id, '/', fsep = "")

    res <- httr::DELETE(project_endpoint,
                        httr::authenticate(client$username, client$token),
                        encode = "json")

    if (httr::status_code(res) == 400) {
        json = httr::content(res, "parsed")
        detail = json$detail
        if (detail == "Project not found") {
            stop("Project not found")
        } else {
            stop(paste("delete_project", detail))
        }
    } else if (httr::status_code(res) == 403) {
        json = httr::content(res, "parsed")
        detail = json$detail
        stop(paste("Project Delete failed", detail))
    } else if (httr::status_code(res) %in% c(200,204)) {
        json = httr::content(res, "parsed")
        project_id = as.integer(json$id)
        return(project_id)
    } else {
        stop("delete_project failed")
    }
}



#' download_project_attachment_rs
#'
#' Download project attachment from the ReadStore API
#' Write the attachment to a file
#' ID or Name must be provided
#' 
#' @param client ReadStore client
#' @param attachment_name Attachment name to download
#' @param outpath Path to write the attachment
#' @param project_id Project ID to return
#' @param project_name Project name to return
download_project_attachment_rs <- function(client,
                                        attachment_name,
                                        outpath,
                                        project_id = NULL,
                                        project_name = NULL) {

    project_attachment_endpoint = file.path(client$endpoint_url, PROJECT_ATTACHMENT_ENDPOINT, fsep = "")
    
    if (is.null(project_id) & is.null(project_name)) {
        stop("project_id or project_name required")
    }

    json_payload = list(attachment_name = attachment_name)         

    if (!is.null(project_id)) {
        json_payload$project_id = project_id
    }
    if (!is.null(project_name)) {
        json_payload$project_name = project_name
    }

    res <- httr::GET(project_attachment_endpoint,
                    query = json_payload,
                    httr::authenticate(client$username, client$token),
                    encode = "json")

    if (!(httr::status_code(res) %in% c(200,204))) {
        stop("get_project failed")
    } else {
        json = httr::content(res, "parsed")
        if (length(json) == 0) {
            stop("download_project_attachment: Attachment not found")
        } else if (length(json) == 1) {
            attachment <- json[[1]]
            binary_data <- base64enc::base64decode(attachment$body)
            writeBin(binary_data, outpath)
        } else {
            stop("download_project_attachment: Multiple attachments found with name")
        }
    }
}

#' download_fq_dataset_attachment_rs
#'
#' Download dataset attachment from the ReadStore API
#' Write the attachment to a file
#' ID or Name must be provided
#' 
#' @param client ReadStore client
#' @param attachment_name Attachment name to download
#' @param outpath Path to write the attachment
#' @param dataset_id Project ID to return
#' @param dataset_name Project name to return
download_fq_dataset_attachment_rs <- function(client,
                                            attachment_name,
                                            outpath,
                                            dataset_id = NULL,
                                            dataset_name = NULL) {

    fq_attach_endpoint = file.path(client$endpoint_url, FQ_ATTACHMENT_ENDPOINT, fsep = "")
    
    if (is.null(dataset_id) & is.null(dataset_name)) {
        stop("dataset_id or dataset_name required")
    }

    json_payload = list(attachment_name = attachment_name)         

    if (!is.null(dataset_id)) {
        json_payload$dataset_id = dataset_id
    }
    if (!is.null(dataset_name)) {
        json_payload$dataset_name = dataset_name
    }

    res <- httr::GET(fq_attach_endpoint,
                    query = json_payload,
                    httr::authenticate(client$username, client$token),
                    encode = "json")

    if (!(httr::status_code(res) %in% c(200,204))) {
        stop("download_fq_dataset_attachment failed")
    } else {
        json = httr::content(res, "parsed")
        if (length(json) == 0) {
            stop("download_fq_dataset_attachment: Attachment not found")
        } else if (length(json) == 1) {
            attachment <- json[[1]]
            binary_data <- base64enc::base64decode(attachment$body)
            writeBin(binary_data, outpath)
        } else {
            stop("download_fq_dataset_attachment: Multiple attachments found with name")
        }
    }
}


#' upload_pro_data_rs
#'
#' Upload ProData entry using the ReadStore API
#' Dataset ID or Name must be provided to select the dataset to attach
#' 
#' Return 403 if the user does not have permission to upload ProData
#' 
#' @param client ReadStore client
#' @param name Name of ProData entry
#' @param pro_data_path Path to the ProData file
#' @param data_type Data type of the ProData entry
#' @param metadata Metadata key value list
#' @param description Description of the ProData entry
#' @param dataset_id Dataset ID to attach the ProData entry
#' @param dataset_name Dataset name to attach the ProData entry
upload_pro_data_rs <- function(client,
                                name,
                                pro_data_file,
                                data_type,
                                metadata = list(),
                                description = "",
                                dataset_id = NULL,
                                dataset_name = NULL) {

    pro_data_endpoint = file.path(client$endpoint_url, PRO_DATA_ENDPOINT, fsep = "")
    
    if (is.null(dataset_id) & is.null(dataset_name)) {
        stop("dataset_id or dataset_name required")
    }

    pro_data_file = normalizePath(pro_data_file)

    json_payload = list(name = name,
                        data_type = data_type,
                        upload_path = pro_data_file,
                        description = description)           

    if (length(metadata) > 0) {
        json_payload$metadata = metadata
    }

    if (!is.null(dataset_id)) {
        json_payload$dataset_id = dataset_id
    }
    if (!is.null(dataset_name)) {
        json_payload$dataset_name = dataset_name
    }

    res <- httr::POST(pro_data_endpoint,
                    body = json_payload,
                    config = httr::authenticate(client$username, client$token),
                    encode = "json")

    if (httr::status_code(res) == 403) {
        json = httr::content(res, "parsed")
        stop(paste("Upload ProData Failed:", json$detail))
    }
    else if (!(httr::status_code(res) %in% c(201,204))) {
        stop("upload_pro_data failed")
    }
}


#' list_pro_data_rs
#'
#' Get list of ProData entries from the ReadStore API
#' include_archived flag also shows archived entries 
#' 
#' @param client ReadStore client
#' @param project_id Project ID to filter
#' @param project_name Project name to filter
#' @param dataset_id Dataset ID to filter
#' @param dataset_name Dataset name to filter
#' @param name Name of ProData entry to filter
#' @param data_type Data type of ProData entry to filter
#' @param include_archived Include archived entries (bool)
#' @return (json) list of ProData entries
list_pro_data_rs <- function(client,
                            project_id = NULL,
                            project_name = NULL,
                            dataset_id = NULL,
                            dataset_name = NULL,
                            name = NULL,
                            data_type = NULL,
                            include_archived = FALSE) {
                            
    pro_data_endpoint = file.path(client$endpoint_url, PRO_DATA_ENDPOINT, fsep = "")
    
    json_payload = list(
        project_id = project_id,
        project_name = project_name,
        dataset_id = dataset_id,
        dataset_name = dataset_name,
        name = name,
        data_type = data_type
    )

    if (include_archived==FALSE) {
        json_payload$valid = TRUE
    }

    res <- httr::GET(pro_data_endpoint,
                    query = json_payload,
                    httr::authenticate(client$username, client$token),
                    encode = "json")

    if (!(httr::status_code(res) %in% c(200,204))) {
        stop("list_pro_data failed")
    }

    json = httr::content(res, "parsed")

    return(json)
}


#' get_pro_data_rs
#'
#' Get single ProData entry from the ReadStore API
#' Return by ProData Id or combination of Name and Dataset ID or Name 
#' If version is NULL, return the latest valid version
#' 
#' @param client ReadStore client
#' @param pro_data_id ProData ID to return
#' @param name Name of ProData entry to filter
#' @param version Version of ProData entry to filter
#' @param dataset_id Dataset ID to filter
#' @param dataset_name Dataset name to filter
#' @return (json) ProData entry
get_pro_data_rs <- function(client,
                            pro_data_id = NULL,
                            name = NULL,
                            version = NULL,
                            dataset_id = NULL,
                            dataset_name = NULL) {
                            
    pro_data_endpoint = file.path(client$endpoint_url, PRO_DATA_ENDPOINT, fsep = "")
    
    if (is.null(pro_data_id)) {
        if (is.null(name) | (is.null(dataset_id) & is.null(dataset_name))) {
            stop("name and dataset_id or dataset_name required")
        }
    }

    if (is.null(version)) {
        valid = 'true'
    } else {
        valid = 'false'
    } 
    
    json_payload = list(
        dataset_id = dataset_id,
        dataset_name = dataset_name,
        name = name,
        version = version,
        valid = valid,
        detail = 'true'
    )

    if (!is.null(pro_data_id)) {
        pro_data_endpoint = file.path(pro_data_endpoint, pro_data_id, '/', fsep = "")
        res = httr::GET(pro_data_endpoint,
                        httr::authenticate(client$username, client$token),
                        encode = "json")
    } else {
        res = httr::GET(pro_data_endpoint,
                        query = json_payload,
                        httr::authenticate(client$username, client$token),
                        encode = "json")
    }

    if (!(httr::status_code(res) %in% c(200,204))) {
        stop("list_pro_data failed")
    } else {
        json = httr::content(res, "parsed")

        if (length(json) == 0) {
            return(list())
        } else if (length(json) == 1) {
           return(json[[1]])
        } else {
            stop("get_pro_data: Multiple ProData found")
        }
    }
}


#' delete_pro_data_rs
#'
#' Delete ProData entry from ReadStore
#' Delete by ProData Id or combination of Name and Dataset ID or Name 
#' If version is NULL, delete the latest valid version
#' 
#' @param client ReadStore client
#' @param pro_data_id ProData ID to return
#' @param name Name of ProData entry to filter
#' @param dataset_id Dataset ID to filter
#' @param dataset_name Dataset name to filter
#' @param version Version of ProData entry to filter
delete_pro_data_rs <- function(client,
                                pro_data_id = NULL,
                                name = NULL,
                                dataset_id = NULL,
                                dataset_name = NULL,
                                version = NULL) {
    
    pro_data_endpoint = file.path(client$endpoint_url, PRO_DATA_ENDPOINT, fsep = "")
    
    if (is.null(pro_data_id)) {
        if (is.null(name) | (is.null(dataset_id) & is.null(dataset_name))) {
            stop("name and dataset_id or dataset_name required")
        }
    }

    json_payload = list(
        dataset_id = dataset_id,
        dataset_name = dataset_name,
        name = name,
        version = version
    )

    if (!is.null(pro_data_id)) {
        pro_data_endpoint = file.path(pro_data_endpoint, pro_data_id, '/', fsep = "")
        res <- httr::DELETE(pro_data_endpoint,
                        httr::authenticate(client$username, client$token),
                        encode = "json")
    } else {
        res <- httr::DELETE(pro_data_endpoint,
                            query = json_payload,
                            httr::authenticate(client$username, client$token),
                            encode = "json")
    }

    if (httr::status_code(res) == 400) {
        json = httr::content(res, "parsed")
        detail = json$detail
        if (detail == "ProData not found") {
            stop("ProData not found")
        } else {
            stop(paste("delete_pro_data", detail))
        }
    } else if (httr::status_code(res) == 403) {
        json = httr::content(res, "parsed")
        detail = json$detail
        stop(paste("ProData Delete failed", detail))
    } else if (httr::status_code(res) %in% c(200,204)) {
        json = httr::content(res, "parsed")
        pro_data_id = as.integer(json$id)
        return(pro_data_id)
    } else {
        stop("delete_pro_data failed")
    }
}
