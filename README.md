# Colocator API

This repository contains the backend for the website Colocator and is generated with django 4.2. <br>
Colocator allows for colocalization analysis to compare genetic association datasets (e.g. GWAS and eQTL). 
With Colocator you are able to upload your own GWAS summary statistics for colocalization analysis.
It returns summary files and locus visualization plots to allow for detailed review of the results

## Data

This project requires two datasets:
- A preprocessed version of the eQTLGen phase 1 dataset in HDF5 format
- A preprocessed version of the Gencode Human Release 43 (GRCh37) file in HDF5 format

For access to these datasets, contact one of the authors.

## Installation

### Software Requirements
To run the code in this repository, the following software is required:
- [Python](https://www.python.org/downloads/) (tested under version 3.8, 3.10, and 3.11)
- [R](https://www.r-project.org/) (tested under version 4.2.3)
- [PostgreSQL](https://www.postgresql.org/download/) (tested under version 15)
- [Redis](https://redis.io/download/) (tested under version 7.0.11)

When all the software above is installed, follow the steps below to install the project.

### Getting started
#### Clone the repository
First, clone the Github repository by running the following command:

```bash
git clone https://github.com/joost011/ColocAPI.git
```

When you have cloned the repository, switch into the newly created directory:

```bash
cd ColocAPI
```

#### Install python packages
After switching to the project directory, run the command below to install the required python packages. It is recommended to do this in a virtual environment.

```bash
pip install -r requirements.txt
```

#### Install R packages
Install the following R packages:
- jsonlite (version 1.8.4)
- coloc (version 5.2.1)
- rjson (version 0.2.21)

#### Create database
After installing the required python and R packages, you have to create a PostgreSQL database. See the [PostgreSQL documentation](https://www.postgresql.org/docs/current/sql-createdatabase.html). This can be accomplished with the command line tool of PostgreSQL called `psql`. Mind that it might be possible that you have to add it manually to you Windows Path when using Windows. The default super user of PostgreSQL is called ```postgres```. To create a database, run the following commands in a new terminal that will create a database called ```coloc```:

```bash
> psql -U postgres
> CREATE DATABASE coloc;
```

After creating the database, add a record to the ```pg_hba.conf``` file in the installation directory of PostgreSQL. See [the documentation](https://www.postgresql.org/docs/current/auth-pg-hba-conf.html) for further explanation.

#### Set environment variables
Next, you have to set the environment variables. The repository contains a file ```.env.example```. Copy this file and remove ```.example``` from the file name, what will result in a new file called ```.env```. Fill in the variable fields. If during installation of PostgreSQL you did not change the port number, the default port is ```5432```.

#### Run migrations
After creating a database and setting the environment variables, run the migrations with the following command (inside the project directory):

```bash
python manage.py migrate
```

The above command will create the necessary tables in the database.

#### Create storage directory
Finally, you have to create the storage directory of the project. In the root of the project, create a directory called ```storage```, containing four subdirectories: ```in_files```, ```out_files```, ```processed_files```, and ```static_files```. Place the two datasets mentioned earlier in this readme in the ```static_files``` directory.

After following the above steps, the installation is completed.

## Usage
### Run the development server
To run the development server, run the following command in the project directory:

```bash
python manage.py runserver
```

### Run the RQ workers
To run the Redis Queue workers, run the following command in a different terminal in the project directory:

```bash
python manage.py rqworker low default high
```

Don't forget to start the Redis server:

```bash
sudo service redis-server start
```

Redis only works on Unix systems. When using Windows, run the two commands in a WSL terminal.

## Authors
Contributers to this repository:
- [ExquizeNyx](https://github.com/ExquizeNyx)
- [JobMaathuis](https://github.com/JobMaathuis)
- [joost011](https://github.com/joost011)

## Acknowledgements 
We want to thank the Hanze university of Applied Sciences and the University Medical Centre Groningen (UMCG) for making this project possible. In particular, we want to thank Lude Franke and Harm-Jan Westra from the UMCG for their expertise and guidancethroughout this project.

## License
Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3. See LICENSE.md for more information.

## Built with
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![DjangoREST](https://img.shields.io/badge/DJANGO-REST-ff1709?style=for-the-badge&logo=django&logoColor=white&color=ff1709&labelColor=gray)
![Postgres](https://img.shields.io/badge/postgres-%23316192.svg?style=for-the-badge&logo=postgresql&logoColor=white)
![Redis](https://img.shields.io/badge/redis-%23DD0031.svg?style=for-the-badge&logo=redis&logoColor=white)


