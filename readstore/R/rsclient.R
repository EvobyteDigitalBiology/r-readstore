library(httr)
library(base64enc)

REST_API_VERSION = "api_v1/"
USER_AUTH_TOKEN_ENDPOINT = "user/auth_token/"
FASTQ_UPLOAD_ENDPOINT = "fq_file_upload/"
FQ_DATASET_ENDPOINT = "fq_dataset/token/"
FQ_FILE_ENDPOINT = "fq_file/token/"
FQ_ATTACHMENT_ENDPOINT = "fq_attachment/token/"
PROJECT_ENDPOINT = "project/token/"
PROJECT_ATTACHMENT_ENDPOINT = "project_attachment/token/"

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

auth_user_token <- function(endpoint, username, token) {

    tryCatch({
    # Construct the authentication endpoint URL
    auth_endpoint <- file.path(endpoint, USER_AUTH_TOKEN_ENDPOINT, fsep = "")
    
    # Prepare the payload
    payload <- list(username = username, token = token)
    
    # Send a POST request
    res <- httr::POST(auth_endpoint, body = payload, encode = "json")
    
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


get_rs_client <- function(username, token, endpoint_url) {

    endpoint_url = file.path(endpoint_url, REST_API_VERSION, fsep = "")

    if (!test_server_connection(endpoint_url)) {
        stop("Connection Error")
    }
    if (!auth_user_token(endpoint_url, username, token)) {
        stop("Authentication Error")
    }

    return(list(username = username, token = token, endpoint_url = endpoint_url))
}

upload_fastq_rs <- function(client, fastq_files) {

    fq_upload_endpoint = file.path(client$endpoint_url, FASTQ_UPLOAD_ENDPOINT, fsep = "")
    
    for(fq in fastq_files) {
        
        fq_path = normalizePath(fq)

        if (!(file.exists(fq_path))) {
            stop(paste("fastq file not found for path", fq_path))
        }
        if (file.access(fq_path,mode=4) != 0) {
            stop(paste("No read permissions for path", fq_path))
        }
        
        json_payload = list(username = client$username,
                            token = client$token,
                            fq_file_path = fq_path)

        res <- httr::POST(fq_upload_endpoint, body = json_payload, encode = "json")

        if (!(httr::status_code(res) %in% c(200, 204))) {
            stop("upload_fastq failed")
        }
    }
}


get_fq_file_rs <- function(client, fq_file_id) {
    # Get a fastq file from the ReadStore API
    #
    # Args:
    #     client: ReadStore client
    #     fq_file_id: Fastq file ID
    #
    # Returns:
    #     Fastq file object

    fq_file_endpoint = file.path(client$endpoint_url, FQ_FILE_ENDPOINT, fsep = "")

    json_payload = list(username = client$username, token = client$token, fq_file_id = fq_file_id)

    res <- httr::POST(fq_file_endpoint, body = json_payload, encode = "json")

    if (!(httr::status_code(res) %in% c(200, 204))) {
        stop("get_fq_file Failed")
    }

    json = httr::content(res, "parsed")[[1]]

    return(json)
}


get_fq_file_upload_path_rs <- function(client, fq_file_id) {
    # Get the upload path for a fastq file from the ReadStore API
    #
    # Args:
    #     client: ReadStore client
    #     fq_file_id: Fastq file ID
    #
    # Returns:
    #     Upload path

    fq_file = get_fq_file(client, fq_file_id)

    if (!("upload_path" %in% names(fq_file))) {
        stop("Upload Path Not Found")
    }

    return(fq_file$upload_path)
}


list_fastq_datasets_rs <- function(client, project_id = NULL, project_name = NULL) {

    fq_dataset_endpoint = file.path(client$endpoint_url, FQ_DATASET_ENDPOINT, fsep = "")

    json_payload = list(username = client$username, token = client$token)

    if (!is.null(project_name)) {
        json_payload$project_name = project_name
    }
    if (!is.null(project_id)) {
        json_payload$project_id = project_id
    }

    res <- httr::POST(fq_dataset_endpoint, body = json_payload, encode = "json")

    if (!(httr::status_code(res) %in% c(200,204))) {
        stop("list_fastq_datasets Failed")
    }

    json = httr::content(res, "parsed")

    return(json)
}


get_fastq_dataset_rs <- function(client, dataset_id = NULL, dataset_name = NULL) {

    fq_dataset_endpoint = file.path(client$endpoint_url, FQ_DATASET_ENDPOINT, fsep = "")

    if (is.null(dataset_id) & is.null(dataset_name)) {
        stop("dataset_id or dataset_name required")
    }

    json_payload = list(username = client$username, token = client$token)

    if (!is.null(dataset_id)) {
        json_payload$dataset_id = dataset_id
    }
    if (!is.null(dataset_name)) {
        json_payload$dataset_name = dataset_name
    }

    res <- httr::POST(fq_dataset_endpoint, body = json_payload, encode = "json")
    
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


list_projects_rs <- function(client) {

    project_endpoint = file.path(client$endpoint_url, PROJECT_ENDPOINT, fsep = "")

    json_payload = list(username = client$username, token = client$token)

    res <- httr::POST(project_endpoint, body = json_payload, encode = "json")

    if (!(httr::status_code(res) %in% c(200,204))) {
        stop("list_projects failed")
    }

    json = httr::content(res, "parsed")

    return(json)
}


get_project_rs <- function(client, project_id = NULL, project_name = NULL) {

    project_endpoint = file.path(client$endpoint_url, PROJECT_ENDPOINT, fsep = "")

    if (is.null(project_id) & is.null(project_name)) {
        stop("project_id or project_name required")
    }

    json_payload = list(username = client$username, token = client$token)

    if (!is.null(project_id)) {
        json_payload$project_id = project_id
    }
    if (!is.null(project_name)) {
        json_payload$project_name = project_name
    }

    res <- httr::POST(project_endpoint, body = json_payload, encode = "json")

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

download_project_attachment_rs <- function(client,
                                        attachment_name,
                                        outpath,
                                        project_id = NULL,
                                        project_name = NULL) {

    project_attachment_endpoint = file.path(client$endpoint_url, PROJECT_ATTACHMENT_ENDPOINT, fsep = "")
    
    if (is.null(project_id) & is.null(project_name)) {
        stop("project_id or project_name required")
    }

    json_payload = list(username = client$username,
                        token = client$token,
                        attachment_name = attachment_name)         

    if (!is.null(project_id)) {
        json_payload$project_id = project_id
    }
    if (!is.null(project_name)) {
        json_payload$project_name = project_name
    }

    res <- httr::POST(project_attachment_endpoint, body = json_payload, encode = "json")

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


download_fq_dataset_attachment_rs <- function(client,
                                        attachment_name,
                                        outpath,
                                        dataset_id = NULL,
                                        dataset_name = NULL) {

    fq_attach_endpoint = file.path(client$endpoint_url, FQ_ATTACHMENT_ENDPOINT, fsep = "")
    
    if (is.null(dataset_id) & is.null(dataset_name)) {
        stop("dataset_id or dataset_name required")
    }

    json_payload = list(username = client$username,
                        token = client$token,
                        attachment_name = attachment_name)         

    if (!is.null(dataset_id)) {
        json_payload$dataset_id = dataset_id
    }
    if (!is.null(dataset_name)) {
        json_payload$dataset_name = dataset_name
    }

    res <- httr::POST(fq_attach_endpoint, body = json_payload, encode = "json")

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