library(dplyr)


RETURN_TYPES = c('json')

check_return_type <- function(return_type) {
    if (!(return_type %in% RETURN_TYPES)) {
        stop('Invalid return_type')
    }
}

get_client <- function(config_dir = '~/.readstore',
                        username = NULL,
                        token = NULL,
                        host = 'http://localhost',
                        return_type = 'json',
                        port = 8000,
                        fastq_extensions = c('.fastq','.fastq.gz','.fq','.fq.gz')) {

    endpoint <- paste0(paste(host, port, sep=':'),'/')

    check_return_type(return_type)

    if (xor(!(is.null(username)),(!(is.null(username))))) {
        stop('Both Username and Token must be provided')

    } else if (file.exists(config_dir)) {
        
        if (grepl('~', config_dir, fixed=TRUE)) {
            config <- path.expand(config_dir)
        } else {
            config <- normalizePath(config_dir)
        }

        config_path <- file.path(config, 'config')

        if (!(file.exists(config_path))) {
            stop(paste('Config file not found at ', config_path))
        }

        rs_config <- load_rs_config(config_path)
        username <- rs_config$username
        token <- rs_config$token
        endpoint <- rs_config$endpoint_url
        fastq_extensions <- rs_config$fastq_extensions
    }

    # Check for env variables taking precedence
    if (!(Sys.getenv("READSTORE_USERNAME") == "")) {
        username <- Sys.getenv("READSTORE_USERNAME")
    }
    if (!(Sys.getenv("READSTORE_TOKEN") == "")) {
        token <- Sys.getenv("READSTORE_TOKEN")
    }
    if (!(Sys.getenv("READSTORE_ENDPOINT_URL") == "")) {
        endpoint <- Sys.getenv("READSTORE_ENDPOINT_URL")
    }
    if (!(Sys.getenv("READSTORE_FASTQ_EXTENSIONS") == "")) {
        fastq_extensions <- strsplit(Sys.getenv("READSTORE_FASTQ_EXTENSIONS"), ",")[[1]]
    }

    if (is.null(username)) {
        stop('No username found in configuration')
    }
    if (is.null(token)) {
        stop('No token found in configuration')
    }

    # Validate connection through rs_client
    rs_client <- get_rs_client(
        username = username,
        token = token,
        endpoint = endpoint
    )

    rs_client$fastq_extensions = fastq_extensions
    rs_client$return_type = return_type

    return(rs_client)
}

list_datasets <- function(client,
                        project_id = NULL,
                        project_name = NULL,
                        return_type = NULL) {

    fq_datasets = list_fastq_datasets_rs(client,project_id,project_name)

    return(fq_datasets)
}

get_dataset <- function(client,
                        dataset_id = NULL,
                        dataset_name = NULL,
                        return_type = NULL) {

    if (is.null(dataset_id) & is.null(dataset_name)) {
        stop("dataset_id or dataset_name required")
    }

    fq_dataset <- get_fastq_dataset_rs(client, dataset_id, dataset_name)

    return(fq_dataset)

    }

get_fastq <- function(client,
                        dataset_id = NULL,
                        dataset_name = NULL,
                        return_type = NULL) {
        
    if (is.null(dataset_id) & is.null(dataset_name)) {
        stop("dataset_id or dataset_name required")
    }

    fq_dataset <- get_fastq_dataset_rs(client, dataset_id, dataset_name)

    if (length(fq_dataset) == 0) {
        return(list())
    } else {
        fq_file_ids = c(fq_dataset$fq_file_r1,
                        fq_dataset$fq_file_r2,
                        fq_dataset$fq_file_i1,
                        fq_dataset$fq_file_i2
                        )


        fq_files = lapply(fq_file_ids, function(x) get_fq_file(client, fq_file_id=x))
        return(fq_files)
    }
}

download_dataset_attachment <- function(
    client,
    attachment_name,
    dataset_id = NULL,
    dataset_name = NULL,
    outpath = NULL) {

    if (is.null(dataset_id) & is.null(dataset_name)) {
        stop("dataset_id or dataset_name required")
    }

    fq_dataset <- get_fastq_dataset_rs(client, dataset_id, dataset_name)

    if (length(fq_dataset) == 0) {
        stop('Dataset not found')
    }

    attachments = fq_dataset$attachments

    if (!(attachment_name %in% attachments)) {
        stop('Attachment not found')
    } else {
        
        if (is.null(outpath)) {
            outpath <- getwd()
            outpath <- file.path(outpath, attachment_name)
        }

        output_dirname <- dirname(outpath)

        if ((output_dirname != '') & (!file.exists(output_dirname))) {
            stop(paste0('Output directory does not exists ', output_dirname))
        }

        download_fq_dataset_attachment(client, attachment_name, outpath, dataset_id, dataset_name)
    }
}


list_projects <- function(client,
                        return_type = NULL) {

    projects <- list_projects_rs(client)

    return(projects)
}


get_project <- function(client,
                        project_id = NULL,
                        project_name = NULL,
                        return_type = NULL) {

    if (is.null(project_id) & is.null(project_name)) {
        stop("project_id or project_name required")
    }

    project <- get_project_rs(client, project_id, project_name)

    return(project)
}


download_project_attachment <- function(
    client,
    attachment_name,
    project_id = NULL,
    project_name = NULL,
    outpath = NULL) {

    if (is.null(project_id) & is.null(project_name)) {
        stop("project_id or project_name required")
    }

    project = get_project_rs(client, project_id, project_name)

    if (length(project) == 0) {
        stop('Project not found')
    }

    attachments = project$attachments

    if (!(attachment_name %in% attachments)) {
        stop('Attachment not found')
    } else {
        
        if (is.null(outpath)) {
            outpath <- getwd()
            outpath <- file.path(outpath, attachment_name)
        }

        output_dirname = dirname(outpath)

        if ((output_dirname != '') & (!file.exists(output_dirname))) {
            stop(paste0('Output directory does not exists ', output_dirname))
        }

        download_project_attachment_rs(client, attachment_name, outpath, project_id, project_name)
    }
}

upload_fastq <- function(client, fastq) {

    upload_files = c()
    for (fq in fastq) {
        if (!(file.exists(fq))) {
            stop(paste('fastq file does not exists', fq))
        }
        
        if (any(sapply(client$fastq_extensions, function(x) endsWith(fq, x)))) {
            stop(paste('fastq file is not a valid FASTQ files ', fq))
        } 

        upload_files = c(upload_files, fq)
    }

    upload_fastq_rs(upload_files)
}