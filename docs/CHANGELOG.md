# Changelog

## [1.4.0] - 2025-02-10

**Features**
- Add method update_dataset
- Add method update_project

**Bugfixes**
- Provided user and token in client definition take precedence over config file

## [1.3.0] - 2025-01-02

**Features**
- Add method create_dataset
- Add method delete_dataset
- Add method create_project
- Add method delete_project

**Updates**
- Update validation of metadata in backend

## [1.2.0] - 2024-12-13

**Features**
- Add method list_datasets_metadata
- Add method list_projects_metadata
- Add method list_pro_data_metadata

**Updates**
- list_datasets: add option to specify data.frame or list as return type
- list_projects: add option to specify data.frame or list as return type
- list_pro_data: add option to specify data.frame or list as return type

**Bugfixes**
Refactor: pro_data_path to pro_data_file

## [1.1.0] - 2024-12-01

**Features**
- Add method upload_pro_data
- Add method list_pro_data
- Add method get_pro_data
- Add method delete_pro_data

**Updates**

- upload_fastq: Add fastq_name and read_type arguments to enable definition of FASTQ names and read types

## [1.0.0] - 2024-11-06

Initial Version