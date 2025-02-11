RETURN_TYPES = c('list', 'data.frame')

#' check_return_type
#' 
#' Stop if return_type is not valid
#' 
#' @param return_type The return type to check 
check_return_type <- function(return_type) {
    if (!(return_type %in% RETURN_TYPES)) {
        stop('Invalid return_type')
    }
}


#' get_metadata_df
#' 
#' Return metadata entries from a list and cast to a data frame
#' Metadata is key value pairs and key values are used to 
#' construct the data.frame names
#' 
#' @param input_list List of metadata entries
#' @return data.frame with metadata entries
get_metadata_df <- function(input_list) {
  
  # Extract metadata as a list
  metadata <- lapply(input_list, function(x) x$metadata)
  all_keys <- unique(unlist(lapply(metadata, names)))
  
  # Convert the metadata to a data frame
  metadata_df <- do.call(rbind, lapply(metadata, function(x) {
    # Create a named vector with NA for missing keys
    row <- setNames(rep(NA, length(all_keys)), all_keys)
    # Populate the row with existing values
    row[names(x)] <- x
    return(as.data.frame(as.list(row), stringsAsFactors = FALSE))
  }))  
  
  return(metadata_df)
}


#' get_client
#'
#' Get ReadStore client to interact with ReadStore API.
#' Providing username and token will take precedence over configuration file.
#' Environment variables will take precedence over configuration file.
#' 
#' @param config_dir Path to the configuration directory to load the configuration from
#' @param username Username to use for authentication
#' @param token Token to use for authentication
#' @param host Host to connect to
#' @param port Port to connect to
#' @param fastq_extensions FASTQ file extensions to use
#' @return ReadStore client (list)
#' @export
get_client <- function(config_dir = '~/.readstore',
                        username = NULL,
                        token = NULL,
                        host = 'http://localhost',
                        port = 8000,
                        fastq_extensions = c('.fastq','.fastq.gz','.fq','.fq.gz')) {

    # TODO: username+token should take precedence over config_dir if provided

    if (xor(!(is.null(username)),(!(is.null(token))))) {
        stop('Both Username and Token must be provided')
    # If username and token are provided, they take precedence over config file
    } else if (!(is.null(username) | is.null(token))) {
        
        username <- username
        token <- token
        fastq_extensions <- fastq_extensions
        endpoint <- paste0(paste(host, port, sep=':'),'/')

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

    # Check if username and token are provided,
    # if true they take precedence over config file


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

    # Make sure endpoint end with trailing slash
    if (!(grepl('/$', endpoint))) {
        endpoint <- paste0(endpoint, '/')
    }

    # Validate connection through rs_client
    rs_client <- get_rs_client(
        username = username,
        token = token,
        endpoint_url = endpoint
    )

    rs_client$fastq_extensions = fastq_extensions

    return(rs_client)
}


#' list_datasets
#'
#' Get list of fastq datasets from the ReadStore API
#' 
#' @param client ReadStore client
#' @param project_id Project ID to filter datasets for
#' @param project_name Project name to filter datasets for
#' @param return_type The return type (list | data.frame)
#' @return data.frame or list of datasets
#' @export
list_datasets <- function(client,
                        project_id = NULL,
                        project_name = NULL,
                        return_type = 'data.frame') {
    
    check_return_type(return_type)

    fq_datasets = list_fastq_datasets_rs(client, project_id, project_name)
    
    if (return_type == 'data.frame') {
        
        # Return empty data.frame if no datasets found
        if (length(fq_datasets) == 0) {
            return(data.frame())
        } else {
            # Extract project metadata
            metadata <- lapply(fq_datasets, function(x) x$metadata) 
            # Extract project attachments
            attachments <- lapply(fq_datasets, function(x) x$attachments) 
            # Extract project
            project_names <- lapply(fq_datasets, function(x) as.character(x$project_names)) 
            # Extract project
            project_ids <- lapply(fq_datasets, function(x) as.integer(x$project_ids))

            # Flatten the dataset list
            # Extract id, name, description, paired_end, index_read metrics
            flattened_df <- do.call(rbind, lapply(fq_datasets, function(x) {
                data.frame(
                    id = x$id,
                    name = x$name,
                    description = x$qc_passed,
                    paired_end = x$paired_end,
                    index_read = x$index_read,
                    stringsAsFactors = FALSE
                    )
                }))

            # Add metadata to the data.frame
            flattened_df$project_ids <- I(project_ids)
            flattened_df$project_names <- I(project_names)

            flattened_df$metadata <- I(metadata)
            flattened_df$attachments <- I(attachments)

            return(flattened_df)
        }
    } else {
        return(fq_datasets)
    }
}


#' list_datasets_metadata
#'
#' Return metadata for datasets from the ReadStore API
#' Order of returned metadata is the same as the dataset list
#' 
#' Return a list of metadata entries for each dataset
#' 
#' or return a data.frame with metadata entries
#' Here metadata keys will be cast to column names (wide format)
#' 
#' @param client ReadStore client
#' @param project_id Subset metadata for a project ID
#' @param project_name Subset metadata for a project name
#' @param return_type The return type (data.frame | list)
#' @return data.frame or list of datasets of metadata
#' @export
list_datasets_metadata <- function(client,
                                    project_id = NULL,
                                    project_name = NULL,
                                    return_type = 'data.frame') {

    check_return_type(return_type)                

    fq_datasets = list_fastq_datasets_rs(client, project_id, project_name)
    
    if (return_type == 'data.frame') {
        if (length(fq_datasets) == 0) {
            return(data.frame())
        } else {
            metadata_df = get_metadata_df(fq_datasets)
            return(metadata_df)
        } 
    } else {
        metadata <- lapply(fq_datasets, function(x) x$metadata) 
        return(metadata)
    }
}

#' get_dataset
#'
#' Get fastq dataset from the ReadStore API
#' ID or Name must be provided
#' 
#' @param client ReadStore client
#' @param dataset_id Dataset ID to return
#' @param dataset_name Dataset name to return
#' @param return_type The return type (currently only list)
#' @return json object (list) with fastq dataset
#' @export
get_dataset <- function(client,
                        dataset_id = NULL,
                        dataset_name = NULL,
                        return_type = 'list') {
    
    check_return_type(return_type)

    if (is.null(dataset_id) & is.null(dataset_name)) {
        stop("dataset_id or dataset_name required")
    }

    fq_dataset <- get_fastq_dataset_rs(client, dataset_id, dataset_name)

    return(fq_dataset)
    }


#' create_dataset
#'
#' Create a new dataset by calling create_fastq_dataset_rs
#'
#' @param client ReadStore client
#' @param name Dataset name
#' @param description Dataset description
#' @param qc_passed QC Pass
#' @param paired_end Paired End
#' @param index_read Index read
#' @param project_ids List of project IDs
#' @param project_names List of project names
#' @param metadata Named list of metadata
#' @return Created dataset object
#' @export
create_dataset <- function(client,
                           dataset_name,
                           description = '',
                           project_ids = list(),
                           project_names = list(),
                           metadata = list()) {
    
    # Check if a dataset with the same name already exists
    dataset_check <- get_dataset(client, dataset_name = dataset_name)

    if (length(dataset_check) > 0) {
        stop('Dataset with the same name already exists')
    }

    # Validate project_ids and project_names
    if (!is.null(project_ids)) {

        if (!is.list(project_ids)) {
            stop('project_ids must be a list')
        }

        existing_project_ids <- sapply(project_ids, function(id) {
            project <- get_project_rs(client, project_id = id)
            return(!is.null(project))
        })
        if (!all(existing_project_ids)) {
            stop('One or more project_ids do not exist in the database')
        } else {
            project_ids <- as.list(project_ids)
        }
    } else {
        project_ids <- list()
    }

    if (!is.null(project_names)) {

        if (!is.list(project_names)) {
            stop('project_names must be a vector')
        }

        existing_project_names <- sapply(project_names, function(name) {
            project <- get_project_rs(client, project_name = name)
            return(!is.null(project))
        })
        if (!all(existing_project_names)) {
            stop('One or more project_names do not exist in the database')
        } else {
            project_names <- as.list(project_names)
        }
    } else {
        project_names <- list()
    }

    res <- create_fastq_dataset_rs(
                        client = client,
                        name = dataset_name,
                        description = description,
                        qc_passed = FALSE,
                        paired_end = FALSE,
                        index_read = FALSE,
                        project_ids = project_ids,
                        project_names = project_names,
                        metadata = metadata)

    return(res)
}


#' update_dataset
#'
#' Update a dataset defined by dataset_id
#' NULL values in the arguments will not be updated
#' 
#' @param client ReadStore client
#' @param dataset_id Dataset ID
#' @param dataset_name Dataset name (NULL or string)
#' @param description Dataset description (NULL or string)
#' @param project_ids List of project IDs (NULL or vector)
#' @param project_names List of project names (NULL or vector)
#' @param metadata Named list of metadata (NULL or list)
#' @return Created dataset object
#' @export
update_dataset <- function(client,
                           dataset_id,
                           dataset_name = NULL,
                           description = NULL,
                           project_ids = NULL,
                           project_names = NULL,
                           metadata = NULL) {
    
    # Check if a dataset with the same name already exists
    dataset_check <- get_dataset(client, dataset_id = dataset_id)

    if (length(dataset_check) == 0) {
        stop('Dataset not found')
    }

    project_ids_new <- NULL
    project_names_new <- NULL
    
    if (is.null(project_ids) & is.null(project_names)) {
        
        # Case both project_ids and project_names are not specified
        # Take over values from the existing dataset
        
        project_ids_new <- dataset_check[['project_ids']]
        project_names_new <- dataset_check[['project_names']]
    } 
    
    if (!is.null(project_ids)) {
        
        # Case project_ids are specified, check if exists

        existing_project_ids <- sapply(project_ids, function(id) {
            project <- get_project_rs(client, project_id = id)
            return(!is.null(project))
        })
        
        if (!all(existing_project_ids)) {
            stop('One or more project_ids do not exist in the database')
        } else {
            project_ids_new <- project_ids
        }
    }

    if (!is.null(project_names)) {
        
        # Check if projects defined by name exist
        existing_project_names <- sapply(project_names, function(name) {
            project <- get_project_rs(client, project_name = name)
            return(!is.null(project))
        })
        if (!all(existing_project_names)) {
            stop('One or more project_names do not exist in the database')
        } else {
            project_names_new <- project_names
        }
    }

    # Set project_ids and project_names to empty list if any is left NULL
    # This is the case if only project_ids OR project_names were specified
    if (is.null(project_ids_new)) {
        project_ids_new <- list()
    }

    if (is.null(project_names_new)) {
        project_names_new <- list()
    }


    # Define 
    dataset_name_new <- ifelse(is.null(dataset_name), dataset_check$name, dataset_name)
    description_new <- ifelse(is.null(description), dataset_check$description, description)
    
    if (is.null(metadata)) {
        metadata_new <- dataset_check[['metadata']]
    } else {
        metadata_new <- metadata
    }

    update_fastq_dataset_rs(client = client,
                            dataset_id = dataset_id,
                            name = dataset_name_new,
                            description = description_new,
                            qc_passed = dataset_check$qc_passed,
                            paired_end = dataset_check$paired_end,
                            index_read = dataset_check$index_read,
                            project_ids = project_ids_new,
                            project_names = project_names_new,
                            metadata = metadata_new)
}

#' delete_dataset
#'
#' Delete a dataset by checking its existence and calling delete_fastq_dataset_rs
#'
#' @param client ReadStore client
#' @param dataset_id Dataset ID
#' @param dataset_name Dataset name
#' @return ID of the deleted dataset
#' @export
delete_dataset <- function(client, dataset_id = NULL, dataset_name = NULL) {
    
    if (is.null(dataset_id) && is.null(dataset_name)) {
        stop("dataset_id or dataset_name is required")
    }
    ds <- get_fastq_dataset_rs(client, dataset_id = dataset_id, dataset_name = dataset_name)
    if (length(ds) == 0 || is.null(ds$id)) {
        stop("No dataset found with the provided dataset_id or dataset_name")
    }

    # Returns integer
    res <- delete_fastq_dataset_rs(client, ds$id)

    return(res)
}





#' get_fastq
#'
#' Get fastq files from the ReadStore API for a dataset
#' ID or Name of dataset must be provided
#' 
#' @param client ReadStore client
#' @param dataset_id Dataset ID to return
#' @param dataset_name Dataset name to return
#' @param return_type The return type (currently only json)
#' @return json object (list) with fastq dataset
#' @export
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


        fq_files = lapply(fq_file_ids, function(x) get_fq_file_rs(client, fq_file_id=x))
        return(fq_files)
    }
}

#' download_dataset_attachment
#'
#' Download attachment from a dataset
#' ID or Name of dataset must be provided
#' Also provide the attachment name to download
#' 
#' @param client ReadStore client
#' @param attachment_name Attachment name to download
#' @param dataset_id Dataset ID to return
#' @param dataset_name Dataset name to return
#' @param outpath Path to save the attachment
#' @export
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

        download_fq_dataset_attachment_rs(client, attachment_name, outpath, dataset_id, dataset_name)
    }
}

#' list_projects
#'
#' Get list of projects from the ReadStore API
#' 
#' @param client ReadStore client
#' @param return_type The return type  (list | data.frame)
#' @return data.frame or list of projects
#' @export
list_projects <- function(client,
                        return_type = 'data.frame') {
    
    check_return_type(return_type)

    projects <- list_projects_rs(client)

    if (return_type == 'data.frame') {
        
        if (length(projects) == 0) {
            return(data.frame())
        } else {

            # Extract project metadata
            metadata <- lapply(projects, function(x) x$metadata) 
            # Extract project attachments
            attachments <- lapply(projects, function(x) x$attachments) 

            # Flatten the dataset list
            # Extract id, name, description, paired_end, index_read metrics
            flattened_df <- do.call(rbind, lapply(projects, function(x) {
                data.frame(
                    id = x$id,
                    name = x$name,
                    stringsAsFactors = FALSE
                    )
                }))

            # Add metadata to the data.frame
            flattened_df$metadata <- I(metadata)
            flattened_df$attachments <- I(attachments)

            return(flattened_df)
        }
    } else {
        return(projects)
    }
}


#' list_projects_metadata
#'
#' Return metadata for projects from the ReadStore API
#' Order of returned projects is the same as the dataset list
#' 
#' Return a list of metadata entries for each project
#' 
#' or return a data.frame with metadata entries
#' Here metadata keys will be cast to column names (wide format)
#' 
#' @param client ReadStore client
#' @param return_type The return type (data.frame | list)
#' @return data.frame or list of datasets of metadata
#' @export
#' 
list_projects_metadata <- function(client,
                                    return_type = 'data.frame') {

    check_return_type(return_type)                

    projects <- list_projects_rs(client)
    
    if (return_type == 'data.frame') {
        if (length(projects) == 0) {
            return(data.frame())
        } else {    
            metadata_df = get_metadata_df(projects)
            return(metadata_df)
        } 
    } else {
        metadata <- lapply(projects, function(x) x$metadata) 
        return(metadata)
    }
}


#' get_projects
#'
#' Get project from the ReadStore API
#' ID or Name of project must be provided
#' 
#' @param client ReadStore client
#' @param project_id Project ID to return
#' @param project_name Project name to return
#' @param return_type The return type (currently only json)
#' @return json object (list) with project
#' @export
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


#' create_project
#'
#' Create a new project by calling create_project_rs
#'
#' @param client ReadStore client
#' @param name Project name
#' @param description Project description
#' @param metadata List of project metadata
#' @param dataset_metadata_keys vector of dataset metadata keys
#' @return The created project object
#' @export
create_project <- function(client,
                           project_name,
                           description = '',
                           metadata = list(),
                           dataset_metadata_keys = c()) {
    
    project_check <- get_project(client, project_name = project_name)
    if (length(project_check) > 0) {
        stop("Project with the same name already exists")
    }

    # Validate dataset_metadata_keys is a vector
    if ((length(dataset_metadata_keys) > 0) & (!is.vector(dataset_metadata_keys))) {
        stop('dataset_metadata_keys must be a vector')
    }

    # Convert dataset_metadata_keys to a list
    dataset_metadata_keys <- as.list(setNames(rep('', length(dataset_metadata_keys)), dataset_metadata_keys))

    result <- create_project_rs(
        client = client,
        name = project_name,
        description = description,
        metadata = metadata,
        dataset_metadata_keys = dataset_metadata_keys
    )

    return(result)
}


#' update_project
#'
#' Update a dataset defined by dataset_id
#' NULL values in the arguments will not be updated
#' 
#' @param client ReadStore client
#' @param project_id Dataset ID
#' @param project_name Dataset name (NULL or string)
#' @param description Dataset description (NULL or string)
#' @param metadata Named list of metadata (NULL or list)
#' @param dataset_metadata_keys Named list of metadata (NULL or vector)
#' @return Created dataset object
#' @export
update_project <- function(client,
                           project_id,
                           project_name = NULL,
                           description = NULL,
                           metadata = NULL,
                           dataset_metadata_keys = NULL) {
    
    # Check if a dataset with the same name already exists
    project_check <- get_project(client, project_id = project_id)

    if (length(project_check) == 0) {
        stop('Project Not Found')
    }

    # Define 
    project_name_new <- ifelse(is.null(project_name), project_check$name, project_name)
    description_new <- ifelse(is.null(description), project_check$description, description)
    
    if (is.null(metadata)) {
        metadata_new <- project_check[['metadata']]
    } else {
        metadata_new <- metadata
    }
    
    if (is.null(dataset_metadata_keys)) {
        dataset_metadata_keys_new <- project_check[['dataset_metadata_keys']]
    } else {
        dataset_metadata_keys_new <- dataset_metadata_keys
        dataset_metadata_keys_new <- as.list(setNames(rep('', length(dataset_metadata_keys_new)), dataset_metadata_keys_new))
    }

    update_project_rs(client = client,
                        project_id = project_id,
                        name = project_name_new,
                        description = description_new,
                        metadata = metadata_new,
                        dataset_metadata_keys = dataset_metadata_keys_new)
}



#' delete_project
#'
#' Delete an existing project by checking its existence and calling delete_project_rs
#'
#' @param client ReadStore client
#' @param project_id Project ID to delete
#' @param project_name Project name to delete
#' @return ID of the deleted project
#' @export
delete_project <- function(client,
    project_id = NULL,
    project_name = NULL) {
    
    if (is.null(project_id) && is.null(project_name)) {
        stop("project_id or project_name is required")
    }
    project <- get_project_rs(client, project_id = project_id, project_name = project_name)
    if (length(project) == 0 || is.null(project$id)) {
        stop("No project found with the provided project_id or project_name")
    }
    result <- delete_project_rs(client, project$id)
    return(result)
}


#' download_project_attachment
#'
#' Download attachment from a project
#' ID or Name of project be provided
#' Also provide the attachment name to download
#' 
#' @param client ReadStore client
#' @param attachment_name Attachment name to download
#' @param project_id Project ID to return
#' @param project_name Project name to return
#' @param outpath Path to save the attachment
#' @export
download_project_attachment <- function(
    client,
    attachment_name,
    project_id = NULL,
    project_name = NULL,
    outpath = NULL) {

    if (is.null(project_id) & is.null(project_name)) {
        stop("project_id or project_name required")
    }

    project <- get_project_rs(client, project_id, project_name)

    if (length(project) == 0) {
        stop('Project not found')
    }

    attachments <- project$attachments

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

        download_project_attachment_rs(client, attachment_name, outpath, project_id, project_name)
    }
}


#' upload_pro_data
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
#' @export
upload_pro_data <- function(client,
                            name,
                            pro_data_file,
                            data_type,
                            metadata = list(),
                            description = "",
                            dataset_id = NULL,
                            dataset_name = NULL) {
    
    if (is.null(dataset_id) & is.null(dataset_name)) {
        stop("dataset_id or dataset_name required")
    }

    if (!(file.exists(pro_data_file))) {
        stop('pro_data_file does not exists')
    }

    upload_pro_data_rs(client,
                        name,
                        pro_data_file,
                        data_type,
                        metadata,
                        description,
                        dataset_id,
                        dataset_name)
}

#' list_pro_data
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
#' @param return_type The return type (list | data.frame)
#' @return data.frame or list of datasets
#' @export
list_pro_data <- function(client,
                            project_id = NULL,
                            project_name = NULL,
                            dataset_id = NULL,
                            dataset_name = NULL,
                            name = NULL,
                            data_type = NULL,
                            include_archived = FALSE,
                            return_type = 'data.frame') {
    
    check_return_type(return_type)

    pro_data <- list_pro_data_rs(client,
                            project_id,
                            project_name,
                            dataset_id,
                            dataset_name,
                            name,
                            data_type,
                            include_archived)
    
    if (return_type == 'data.frame') {
        
        # If ProData is empty return empty data.frame
        if (length(pro_data) == 0) {
            return(data.frame())
        } else {

            metadata <- lapply(pro_data, function(x) x$metadata) 

            flattened_df <- do.call(rbind, lapply(pro_data, function(x) {
                data.frame(
                    id = x$id,
                    name = x$name,
                    data_type = x$data_type,
                    version = x$version,
                    dataset_id = x$dataset_id,
                    dataset_name = x$dataset_name,
                    upload_path = x$upload_path,
                    stringsAsFactors = FALSE
                )
            }))

            flattened_df$metadata <- I(metadata)

            return(flattened_df)
        }
    } else {
        return(pro_data)
    }   
    return(pro_data)
}

#' list_pro_data_metadata
#'
#' Return metadata for processed data (pro_data) from the ReadStore API
#' Order of the returned metadata is the same as the pro_data list
#' 
#' Return a list of pro_data entries for each dataset
#' 
#' or return a data.frame with metadata entries
#' Here metadata keys will be cast to column names (wide format)
#' 
#' @param client ReadStore client
#' @param project_id Project ID to filter
#' @param project_name Project name to filter
#' @param dataset_id Dataset ID to filter
#' @param dataset_name Dataset name to filter
#' @param name Name of ProData entry to filter
#' @param data_type Data type of ProData entry to filter
#' @param include_archived Include archived entries (bool)
#' @param return_type The return type (list | data.frame)
#' @return data.frame or list of datasets
#' @export
list_pro_data_metadata <- function(client,
                                    project_id = NULL,
                                    project_name = NULL,
                                    dataset_id = NULL,
                                    dataset_name = NULL,
                                    name = NULL,
                                    data_type = NULL,
                                    include_archived = FALSE,
                                    return_type = 'data.frame') {

    check_return_type(return_type)                

    pro_data <- list_pro_data_rs(client,
                            project_id,
                            project_name,
                            dataset_id,
                            dataset_name,
                            name,
                            data_type,
                            include_archived)
    
    if (return_type == 'data.frame') {
        if (length(pro_data) == 0) {
            return(data.frame())
        } else {
            metadata_df = get_metadata_df(pro_data)
            return(metadata_df)
        } 
    } else {
        metadata <- lapply(pro_data, function(x) x$metadata) 
        return(metadata)
    }
}


#' get_pro_data
#'
#' Get single Processed Data (ProData) entry from the ReadStore API
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
#' @export
get_pro_data <- function(client,
                        pro_data_id = NULL,
                        name = NULL,
                        version = NULL,
                        dataset_id = NULL,
                        dataset_name = NULL) {
    
    if (is.null(pro_data_id)) {
        if (is.null(name) | (is.null(dataset_id) & is.null(dataset_name))) {
            stop("name and dataset_id or dataset_name required")
        }
    }

    pro_data <- get_pro_data_rs(client,
                        pro_data_id,
                        name,
                        version,
                        dataset_id,
                        dataset_name)

    return(pro_data)
}


#' delete_pro_data
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
#' @export
delete_pro_data <- function(client,
                            pro_data_id = NULL,
                            name = NULL,
                            dataset_id = NULL,
                            dataset_name = NULL,
                            version = NULL) {
    
    if (is.null(pro_data_id)) {
        if (is.null(name) | (is.null(dataset_id) & is.null(dataset_name))) {
            stop("name and dataset_id or dataset_name required")
        }
    }

    delete_pro_data_rs(client,
                        pro_data_id,
                        name,
                        version,
                        dataset_id,
                        dataset_name)
}

#' upload_fastq
#'
#' Upload fastq files (fastq) to the ReadStore API
#' 
#' If fastq_name is provided, it must have the same length as fastq
#' If read_type is provided, it must have the same length as fastq
#'
#' Arguments can be provided as strings or vectors 
#' 
#' @param client ReadStore client
#' @param fastq string or vector of fastq files to upload
#' @param fastq_name optional string or vector of fastq names to upload
#' @param read_type optional string or vector of read types to upload
#' @export
upload_fastq <- function(client, fastq, fastq_name = NULL, read_type = NULL) {

    # Check if fastq is a string, if so cast to vector
    if (is.character(fastq)) {
        fastq = c(fastq)
    }
    if (is.character(fastq_name)) {
        fastq_name = c(fastq_name)
    }
    if (is.character(read_type)) {
        read_type = c(read_type)
    }

    if (!(is.null(fastq_name))) {
        if (length(fastq) != length(fastq_name)) {
            stop('fastq and fastq_name must have the same length')
        }
    }
    if (!(is.null(read_type))) {
        if (length(fastq) != length(read_type)) {
            stop('fastq and read_type must have the same length')
        }
    }

    ix = 1
    upload_files = c()
    upload_fq_names = c()
    upload_read_types = c()

    for (fq in fastq) {
        if (!(file.exists(fq))) {
            stop(paste('fastq file does not exists', fq))
        }
        
        if (!(any(sapply(client$fastq_extensions, function(x) endsWith(fq, x))))) {          
            stop(paste('fastq file is not a valid FASTQ files ', fq))
        }

        if (!(is.null(fastq_name))) {
            upload_fq_names = c(upload_fq_names, fastq_name[ix])
        } else {
            upload_fq_names = c(upload_fq_names, NULL)
        }

        if (!(is.null(read_type))) {
            upload_read_types = c(upload_read_types, read_type[ix])
        } else {
            upload_read_types = c(upload_read_types, NULL)
        }
        
        upload_files = c(upload_files, fq)
        ix = ix + 1
    }

    # Perform the upload
    ix = 1
    for (fq in upload_files) {
        upload_fastq_rs(client, fq, upload_fq_names[ix], upload_read_types[ix])
        ix = ix + 1
    }
}