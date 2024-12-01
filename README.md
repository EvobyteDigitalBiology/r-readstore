# r-readstore SDK

This README describes r-readstore, the R client (SDK) for the ReadStore API. 

r-readstore can be used to access projects, datasets, metadata and attachment files in the ReadStore Database from  Python code. 
The package enables you to automate your bioinformatics pipelines, Python scripts and notebooks.

Check the [ReadStore Github repository](https://github.com/EvobyteDigitalBiology/readstore) for more information on how to get started with ReadStore and setting up your server.

More infos on the [ReadStore website](#https://evo-byte.com/readstore/)

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

## The Lean Solution for Managing FASTQ and NGS Data

ReadStore is a platform for storing, managing, and integrating genomic data. It speeds up analysis and offers a simple way of managing and sharing FASTQ and NGS datasets.
Built-in project and metadata management structures your workflows and a collaborative user interface enhances teamwork — so you can focus on generating insights.

The integrated Webservice enables your to directly retrieve data from ReadStore via the terminal Command-Line-Interface (CLI) or Python/R SDKs.

The ReadStore Basic version provides a local webserver with a simple user management. If you need an organization-wide deployment, advanced user and group management or cloud integration please check the ReadStore Advanced versions and reach out to info@evo-byte.com.

## Description

r-readstore is a R client (SDK) that lets you easily connect to your ReadStore server and interact with the ReadStore API.
By importing the readstore package in R, you can quickly retrieve data from a ReadStore server.

This tool provides streamlined and standardized access to NGS datasets and metadata, helping you run analyses more efficiently and with fewer errors.
You can easily scale your pipelines, and if you need to migrate or move NGS data, updating the ReadStore database ensures all your workflows stay up-to-date.


## Security and Permissions<a id="backup"></a>

**PLEASE READ AND FOLLOW THESE INSTRUCTIONS CAREFULLY!**

### User Accounts and Token<a id="token"></a>

Using r-readstore requires an active user account and a token (and a running ReadStore server). 

You should **never enter your user account password** when working with r-readstore.

To retrieve your token:

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
[GitHub repository](#https://github.com/EvobyteDigitalBiology/r-readstore)

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
You can enter (`http://127.0.0.1:8000/api_v1/`) in your browser and should get a response from the API.

We assume you ran `readstore configure` before to create a config file for your user.
If not, consult the [ReadStore CLI](https://github.com/EvobyteDigitalBiology/readstore-cli) README on how to set this up.

We first will create the R client and perform some operations to retrieve data from the ReadStore database.
More information on all available methods can be found below.


```r
library(readstore)

client <- get_client()  # Create an instance of the ReadStore client

datasets = list_datasets(client)          # List all datasets and return json-style list of list

datasets_project_1 = list_datasets(client,          # List all datasets for project 1
                                    project_id = 1) # return json-style list of list
                                                    

datasets_id_25 = get_dataset(client,                # Get detailed data for dataset 25
                             dataset_id = 25)       # return json-style list

fastq_files_dataset_25 = get_fastq(client,          # Get individual fastq files for dataset 25
                                   dataset_id = 25) # return json-style list of list

projects = list_projects(client)          # List all projects and return json-style list of list

projects = get_project(client,                      # Get details for MyProject
                       project_name = 'MyProject')  # return json-style list

download_dataset_attachment(client,                 # Download file attached to dataset 25
                            dataset_id = 25,        
                            attachment_name = 'gene_table.tsv') 

download_project_attachment(client,                 # Download file attached to project
                            project_name = 'MyProject'
                            attachment_name = 'project_plan.pptx')

upload_fastq(client, fastq = c('path/to/fastq_R1.fq', 'path/to/fastq_R2.fq'))

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

### Access Datasets

```r
# List ReadStore Datasets

list_datasets(client,
              project_id = NULL,   # Filter datasets for project with id `project_id`
              project_name = NULL) # Filter datasets for project with name `project_name`
                                   # Return json style list of lists

# Get ReadStore Dataset Details
# Must provide dataset_id OR dataset_name

get_dataset(client,
            dataset_id = NULL,
            dataset_name = NULL) # Return json style list

# Get FASTQ file data for a dataset
# Must provide dataset_id OR dataset_name

get_fastq(client,
          dataset_id = NULL,       # Get fastq data for dataset with id `dataset_id`
          dataset_name = NULL)     # Get fastq data for dataset `dataset_name`
                                   # Return  json style list of lists
```


### Access Projects

```r
# List ReadStore Projects

list_projects(client) # Return  json style list of lists

# Get ReadStore Project Details
# Must provide project_id OR project_name

get_project(client,
            project_id = NULL,    # Get project with id `project_id`
            project_name = NULL)  # Filter project with name `project_name`
```

### Download Attachmeents

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

upload_fastq(client, fastq = c('myfastq_r1.fq', 'myfastq_r2.fq'))  # Vector of FASTQ file path to upload
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