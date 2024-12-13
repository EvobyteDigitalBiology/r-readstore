# r-readstore SDK

This README describes r-readstore, the R client (SDK) for the ReadStore API. 

r-readstore can be used to access projects, datasets, metadata and attachment files in the ReadStore Database from  Python code. 
The package enables you to automate your bioinformatics pipelines, Python scripts and notebooks.

Check the [ReadStore Github repository](https://github.com/EvobyteDigitalBiology/readstore) for more information on how to get started with ReadStore and setting up your server.

More infos on the [ReadStore website](https://evo-byte.com/readstore/)

Tutorials and Intro Videos: https://www.youtube.com/@evobytedigitalbio

Blog posts and How-Tos: https://evo-byte.com/blog/

For general questions reach out to info@evo-byte.com

Happy analysis :)

## Table of Contents
- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Credits and Acknowledgments](#acknowledgments)

## The Lean Solution for Managing NGS and Omics Data

ReadStore is a platform for storing, managing, and integrating omics data. It speeds up analysis and offers a simple way of managing and sharing NGS omics datasets, metadata and processed data (**Pro**cessed **Data**).
Built-in project and metadata management structures your workflows and a collaborative user interface enhances teamwork — so you can focus on generating insights.

The integrated Webservice (API) enables your to directly retrieve data from ReadStore via the terminal [Command-Line-Interface (CLI)](https://github.com/EvobyteDigitalBiology/readstore-cli) or [Python](https://github.com/EvobyteDigitalBiology/pyreadstore) / [R](https://github.com/EvobyteDigitalBiology/r-readstore) SDKs.

The ReadStore Basic version provides a local webserver with a simple user management. If you need an organization-wide deployment, advanced user and group management or cloud integration please check the ReadStore Advanced versions and reach out to info@evo-byte.com.

## Description

r-readstore is a R client (SDK) that lets you easily connect to your ReadStore server and interact with the ReadStore API.
By importing the readstore package in R, you can quickly retrieve data from a ReadStore server.

This tool provides streamlined and standardized access to NGS datasets and metadata, helping you run analyses more efficiently and with fewer errors.
You can easily scale your pipelines, and if you need to migrate or move NGS data, updating the ReadStore database ensures all your workflows stay up-to-date.

## Security and Permissions<a id="backup"></a>

**PLEASE READ AND FOLLOW THESE INSTRUCTIONS CAREFULLY!**

### User Accounts and Token<a id="token"></a>

Using **r-readstore** requires an active user account and a token (and a running ReadStore server). 

You should **never enter your user account password** when working with r-readstore.

To retrieve your token

1. Login to the ReadStore app via your browser
2. Navigate to `Settings` page and click on `Token`
3. You can regenerate your token anytime (`Reset`). This will invalidate the previous token

For uploading FASTQ files your user account needs to have `Staging Permission`.
You can check this in the `Settings` page of your account.
If you not have `Staging Permission`, ask your ReadStore server admin to grant you permission.

### Setting Your Credentials

You need to provide the r-readstore client with valid ReadStore credentials.

There are different options

1. Load credentials from the ReadStore `config` file. 
The file is generated by the [ReadStore CLI](https://github.com/EvobyteDigitalBiology/readstore-cli),
by default in your home directory (`~/.readstore/`). Make sure to keep read permissions to the file restrictive

2. Directly enter your username and token when instantiating a r-readstore client within your R code

3. Set username and token via environment variables (`READSTORE_USERNAME`, `READSTORE_TOKEN`). This is useful in container or cloud environments.


## Installation

In your R environment you can directly install the readstore client from the r-readstore
[GitHub repository](https://github.com/EvobyteDigitalBiology/r-readstore)

```r
# with remotes library
library(remotes)
remotes::install_github('https://github.com/EvobyteDigitalBiology/r-readstore', subdir='readstore')

# or with devtools library
library(devtools)
devtools::install_github('https://github.com/EvobyteDigitalBiology/r-readstore', subdir='readstore')
```

Validate the successful install by running

```r
library(readstore)
```

## Usage

Detailed tutorials, videos and explanations are found on [YouTube](https://www.youtube.com/playlist?list=PLk-WMGySW9ySUfZU25NyA5YgzmHQ7yquv) or on the [**EVO**BYTE blog](https://evo-byte.com/blog).

### Quickstart

Let's access some dataset and project data from the ReadStore database!

Make sure a ReadStore server is running and reachable (by default under `127.0.0.1:8000`).
You can enter (`http://127.0.0.1:8000/api_x_v1/`) in your browser and should get a response from the API.

We assume you ran `readstore configure` before to create a config file for your user.
If not, consult the [ReadStore CLI](https://github.com/EvobyteDigitalBiology/readstore-cli) README on how to set this up.

We first will create the R client and perform some operations to retrieve data from the ReadStore database.
More information on all available methods can be found below.

```r
library(readstore)

client <- get_client()  # Create an instance of the ReadStore client

# Manage Datasets 

datasets <- list_datasets(client)          # List all datasets and return json-style list of list

datasets_project_1 <- list_datasets(client,          # List all datasets for project 1
                                    project_id = 1)  # return data.frame
                                                    
datasets_metadata <- list_datasets_metadata(client)  # Get metadata for datasets cast to data.frame
                                                     # metadata keys as column names

datasets_id_25 <- get_dataset(client,                # Get detailed data for dataset 25
                             dataset_id = 25)        # return json-style R list

fastq_files_dataset_25 <- get_fastq(client,          # Get individual fastq files for dataset 25
                                   dataset_id = 25)  # return json-style R nested list

download_dataset_attachment(client,                 # Download file attached to dataset 25
                            dataset_id = 25,        
                            attachment_name = 'gene_table.tsv') 

# Manage Projects

projects <- list_projects(client)          # List all projects and return data.frame

projects_metadata <- list_projects_metadata(client) # Get metadata for projects cast to data.frame
                                                    # metadata keys as column names

projects <- get_project(client,                         # Get details for MyProject
                       project_name = 'MyProject')      # return json-style list

download_project_attachment(client,                     # Download file attached to project
                            project_name = 'MyProject'
                            attachment_name = 'project_plan.pptx')

# Upload FASTQ datasets

upload_fastq(client,                                                 # Upload FASTQ files
            fastq = c('path/to/fastq_R1.fq', 'path/to/fastq_R2.fq'), # Define path to FASTQ  
            fastq_name = c('sample_R1', 'sample_R2'),                # Set names of FASTQ files
            read_type = c('R1', 'R2'))                               # Define type of Read

# Manage ProData

upload_pro_data(client,                                         # Upload Processed Data files
                name = 'sample_1_count_matrix',                 # Set name of count matrix
                pro_data_path = 'path/to/sample_1_counts.h5',   # Set file path
                data_type = 'count_matrix',                     # Set type to 'count_matrix'
                dataset_id = 25)                                # Attach ProData to dataset_id 25

pro_data_list <- list_pro_data(client,                # List Processed Data (ProData)
                            dataset_id = 25)          # Get ProData for dataset_id 25

pro_data_metadata <- list_pro_data(client,                   # List ProData metadata as data.frame
                                   dataset_id = 25)          # For dataset_id 25

pro_data <- get_pro_data(client,                         # Get individual ProData entry
                        name = 'sample_1_count_matrix',  # Get ProData with name 'sample_1_count_matrix'
                        dataset_id = 25)                 # Get ProData for dataset_id 25

pro_data <- delete_pro_data(client,                         # Delete ProData
                           name = 'sample_1_count_matrix',  # Get ProData with name 'sample_1_count_matrix'
                           dataset_id = 25)                 # Get ProData for dataset_id 25
```


### Configure the R Client

The Client is the central object and provides authentication against the ReadStore API.
By default, the client will try to read the `~/.readstore/config` credentials file.
You can change the directory if your config file is located in another folder.

If you set the `username` and `token` arguments, the client will use these credentials instead.

If your ReadStore server is not running under localhost (`127.0.0.1`) port `8000`, you can adapt the default settings.

```r 
client <- get_client(config_dir = '~/.readstore',  # Directory containing ReadStore credentials
                    username = NULL,               # Username
                    token = NULL,                  # Token
                    host = 'http://localhost',     # Hostname / IP of ReadStore server
                    port = 8000,                   # Server Port Number
                    fastq_extensions= c('.fastq','.fastq.gz','.fq','.fq.gz')) 
                    # Accepted FASTQ file extensions for upload validation 
```

Is is possible to set userame, token, server endpoint and fastq extensions using the listed environment variables. 
The enironment variables precede over other client configurations.

- `READSTORE_USERNAME` (username)
- `READSTORE_TOKEN`    (token)
- `READSTORE_ENDPOINT_URL` (`http://host:post`, e.g. `http://localhost:8000`)
- `READSTORE_FASTQ_EXTENSIONS` (fastq_extensions, `'.fastq',.fastq.gz,.fq,.fq.gz'`)

**Possible errors**

    - Connection Error:     If no ReadStore server was found at the provided endpoint
    - Authentication Error: If provided username or token are not found
    - No Permission to Upload/Delete FASTQ/ProData: User has no [Staging Permissions]

### Access Datasets

```r
# List ReadStore Datasets

# Option to subset by project_id OR project_name
# Option to return either a data.frame or list 

list_datasets(client,
              project_id = NULL,            # Filter datasets for project with id `project_id`
              project_name = NULL           # Filter datasets for project with name `project_name`
              return_type = 'data.frame')   # return_type (data.frame | list)

# List ReadStore Dataset Metadata

# Returns formatted metadata for each dataset
# Metadata keys are cast into data.frame columns (wide-format)

# Option to subset by project_id OR project_name
# Option to return either a data.frame or list 

list_datasets_metadata(client,
                    project_id = NULL,            # Filter metadata for project with id `project_id`
                    project_name = NULL           # Filter metadata for project with name `project_name`
                    return_type = 'data.frame')   # return_type (data.frame | list)

# Get ReadStore Dataset Details

# MUST provide dataset_id OR dataset_name
# Option to return either a data.frame or list 

get_dataset(client,
            dataset_id = NULL,
            dataset_name = NULL) # Return json style list


# Get FASTQ file(s)for a specific dataset

# Must provide dataset_id OR dataset_name

get_fastq(client,
          dataset_id = NULL,       # Get fastq data for dataset with id `dataset_id`
          dataset_name = NULL)     # Get fastq data for dataset `dataset_name`
                                   # Return  json style list of lists
```


### Access Projects

```r
# List ReadStore Projects

list_projects(client,
              return_type = 'data.frame') # Return type (data.frame | list)


# Get Metadata for Projects from ReadStore

# Return metadata for projects from the ReadStore API
# Order of returned projects is the same as the dataset list

# Return a list of metadata entries for each project

# or return a data.frame with metadata entries
# Here metadata keys will be cast to column names (wide format)

list_projects_metadata(client,
                       return_type = 'data.frame')  # The return type (data.frame | list)


# Get ReadStore Project Details

# Must provide project_id OR project_name

get_project(client,
            project_id = NULL,    # Get project with id `project_id`
            project_name = NULL)  # Filter project with name `project_name`
```

### Access ProData

```r
# Upload ProData

upload_pro_data(client,
                name,                # name of ProData entry
                pro_data_path,       # Path to Processed Data file
                data_type,           # Type of dataset (e.g. count_matrix)
                metadata = list(),   # Metadata key values list for ProData
                description = "",    # Set description
                dataset_id = NULL,   # Dataset ID to attach ProData to 
                dataset_name = NULL) # Dataset Name to attach ProData to

# List ProData

list_pro_data(client,
                project_id = NULL,       # Filter by Project ID
                project_name = NULL,     # Filter by Project Name
                dataset_id = NULL,       # Filter by Dataset ID
                dataset_name = NULL,     # Filter by Dataset Name
                name = NULL,             # Filter by Name
                data_type = NULL,        # Filter by Data Type
                include_archived = FALSE, # Return archived ProData
                return_type = 'data.frame') # The return type (data.frame | list)


# List ProData Metadata

list_pro_data_metadata(client,
                        project_id = NULL,       # Filter by Project ID
                        project_name = NULL,     # Filter by Project Name
                        dataset_id = NULL,       # Filter by Dataset ID
                        dataset_name = NULL,     # Filter by Dataset Name
                        name = NULL,             # Filter by Name
                        data_type = NULL,        # Filter by Data Type
                        include_archived = FALSE,# Return archived ProData
                        return_type = 'data.frame') # The return type (data.frame | list)

# This function returns metadata for processed data (ProData) from the ReadStore API. The order of the returned metadata is the same as the ProData list. You can choose to return a list of metadata entries for each dataset or a data.frame with metadata entries where metadata keys will be cast to column names (wide format).


# Get ProData

get_pro_data(client,
            pro_data_id = NULL,     # Get ProData by ID
            name = NULL,            # Get ProData by Name
            version = NULL,         # Get Specific Version
            dataset_id = NULL,      # Get ProData by dataset_id
            dataset_name = NULL)    # Get ProData by dataset_name

# Must provide pro_data_id OR (dataset_name/dataset_id and name)

# Delete ProData

delete_pro_data(client,                 
                pro_data_id = NULL,  # Get ProData by ID
                name = NULL,         # Get ProData by Name    
                dataset_id = NULL,   # Get ProData by dataset_id
                dataset_name = NULL, # Get ProData by dataset_name
                version = NULL)      # Set specific version to delete

# Must provide pro_data_id OR (dataset_name/dataset_id and name)
```


### Download Attachments

```r 
# Download project attachment file from ReadStore Database

# Must provide project_id OR project_name

download_project_attachment(client,
                            attachment_name,          # name of attachment file
                            project_id = NULL,        # project id with attachment
                            project_name = NULL,      # project name with attachment
                            outpath = NULL)           # Path to download file to
                                                      # default NULL download to working dir

# Download dataset attachment file from ReadStore Database 

# Must provide dataset_id OR dataset_name

download_dataset_attachment(client,
                            attachment_name,             # name of attachment file
                            dataset_id = NULL,           # datatset id with attachment
                            dataset_name = NULL,         # datatset name with attachment
                            outpath = NULL)              # Path to download file to
```


### Upload FASTQ files

Upload FASTQ files to ReadStore server. The methods checks if the FASTQ files exist and end with valid FASTQ ending.

```r 
# Upload FASTQ files to ReadStore 

upload_fastq(client,                                         
             fastq,             # Path to FASTQ file (string/vector)
             fastq_name = NULL, # Names of FASTQ files (string/vector)
             read_type = NULL)  # Read Types (string/vector)
```

## Contributing

Contributions make this project better! Whether you want to report a bug, improve documentation, or add new features, any help is welcomed!

### How You Can Help
- Report Bugs
- Suggest Features
- Improve Documentation
- Code Contributions

### Contribution Workflow
1. Fork the repository and create a new branch for each contribution.
2. Write clear, concise commit messages.
3. Submit a pull request and wait for review.

Thank you for helping make this project better!

## License

The r-readstore is licensed under an GLP3 Open Source License.
See the LICENSE file for more information.

## Credits and Acknowledgments<a id="acknowledgments"></a>

r-readstore is built upon the following open-source python packages and would like to thank all contributing authors, developers and partners.

- R (https://www.r-project.org/)
- httr (https://httr.r-lib.org/)
- base64enc (https://github.com/s-u/base64enc)
- ini (https://cran.r-project.org/web/packages/ini/index.html)
- jsonlite (https://cran.r-project.org/web/packages/jsonlite/index.html)