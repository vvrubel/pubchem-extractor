# molecad

API for retrieving data from Pubchem databases and performing substructure search through local mongo database.

## Table of content

* [Confluence documentation](#confluence-documentation)
* [Installation](#installation)
* [Code formatting](#code-formatting)
* [Environment variables](#environment-variables)
* [How to](#how-to)
  * [How to set environmental variable](#how-to-set-environmental-variable)
  * [How to add a new frontend page](#how-to-add-a-new-frontend-page)
  * [How to deploy](#how-to-deploy)
* [Project structure](#project-structure)

## Confluence documentation - TODO
* [API]() 
* [Database schema]()

## Installation
1. Create conda virtual environment
```console
$ conda env create -f conda.yml
$ conda activate molecad
```
2. Install dependencies
```console
$ poetry install
```

## Code formatting

For constant code style we use multiple tools to enforce it (mypy, black, isort, flake8). To avoid running these tools manually all the time we use pre-commit hooks. To install and run it run:

```shell
poetry run pre-commit install
poetry run pre-commit run --all-files
```

## Environment variables

The application needs a set of environment variables to run. There are 3 types of such variables:
1. **Necessary**. You must set these variables before running the application. If some of these variables are missing, an error will be thrown
2. **Important**. You should think about these variables in production. But during the development, you probably can live with defaults
3. **Optional**. It's good to know about the existence of these variables. But their default values are probably good enough

| Variable | Type | Description |
| -------- | ---- | ----------- |
| ONEQ_TOKEN | Necessary | Your One-q token. You can not run One-q tasks without it. You can get your token [here]( https://one-q-api.net.biocad.ru/login) |
| SERVICE_URL | Necessary | URL of the running application. You can run the app without it, but then you can not catch any webhooks. So even for the local development, it is recommended to set your API from BIOCAD network (e.g. `http://10.249.233.40:5000/`). Don't forget to set port (i. e. 5000) and protocol (http or https) |
| MONGO_HOST | Important | A host for mongo DB. For example, `mongo.net.biocad.ru/lechatelier-test`. If not set, a localhost will be used |
| MONGO_USER | Important | Mongo DB username. Not necessary for the local development: if not set, a standard user will be used |
| MONGO_PASSWORD | Important | Mongo DB password |
| MONGO_DB_NAME | Important | Database name. The default value is "chem-mol-props" |
| MONGO_DB_COLLECTION | Important | Collection name. The default value is "predictions" |
| MONGO_PORT | Optional | Port for the database. The default value is 27017 |
| SECRET_KEY | Optional | Random string that is used by Flask for cryptographical [protection](https://flask.palletsprojects.com/en/1.1.x/config/#SECRET_KEY). You can leave default value for the local development. For the production, please generate random key at least once. Keep this key, well, secret |
| FLASK_ENV | Optional | Can be "production" or "development". It influences how logs look like. It is recommended to check locally, how the application works when this variable set to  "production", before the deployment|
| APP_VERSION | Optional | Prefix for all the URLs. The default value is /v0/" |
| ONEQ_CPU_PER_TASK | Optional | How much CPU is required by each One-Q task. The default value is 1 |
| ONEQ_GPU_PER_TASK | Optional | How much GPU is required by each One-Q task. The default value is 0 |
| ONEQ_URL | Optional | The default value is https://one-q-api.prod.net.biocad.ru |
| UPLOAD_FOLDER | Optional | Folder to upload files. The default value is "." |
## Environment variables
### Setting env variables:

Copy `.env.sample` file  into `.env`. Fill missing values. Here are suggestions for some of them. For local development you need to run mongodb database on your machine. [Read here](https://docs.mongodb.com/manual/tutorial/install-mongodb-on-os-x/) on how to do it.

1. Save your IP-address from BIOCAD network to the `SERVICE_URL` environment variable. Set the port 5000:
```console
$ export SERVICE_URL=<your_IP>:5000/
```
For example:
```console
$ export SERVICE_URL=http://10.249.233.40:5000/
```
4. Export your One-q token to the `ONEQ_TOKEN` environment variable. You can get your token [here]( https://one-q-api.net.biocad.ru/login)
```console
$ export ONEQ_TOKEN=<your token>
```

The table of all environment variables for the project can be found in this document below.

5. Run local build:
```console
ENV=development poetry run uvicorn src:app --port 5000 --reload --host 0.0.0.0
```

Now you can use API and develop locally. API documentation can be found on [confluence](https://confluence.biocad.ru/display/CHEMB/API)
## How to
Few useful guides, which might help you to develop a project

* [How to set environmental variable](#how-to-set-environmental-variable)
* [How to add a new frontend page](#how-to-add-a-new-frontend-page)
* [How to deploy](#how-to-deploy)

### How to set environmental variable

It is very convenient to set these variables once in `.env` file. Create such file in the project root, and fill it with any variables that you need. You can find the example in `.env.sample` file.

If you want to add a new environment variable to the project, you should do the following steps:
1. Set in in the console, or in `.env` file. Example of setting the variable in bash console:
```console
$ export GLORIOUS_VAR="beautiful_value"
$ echo $GLORIOUS_VAR # Checking that the variable is set
beautiful_value
```

2. Read it in the `chem_mol_props/setup_.py` `Settings` class. For example, if you added variable "GLORIOUS_VAR", add this code to `Settings class`:
```python
glorious_var: str = Field("shiny_default_value", env="GLORIOUS_VAR")
```
GLORIOUS_VAR will be set to "shiny default value", if it is not set in the environment or in the `.env` file 

3. The above steps will allow you to use new variable in local development. For the production, you also need to add information about this variable in `iac/template-nomad.hcl` file. Do the following:
* Find section `job "${CI_PROJECT_NAME}-${NOMAD_DC}"` (it is probably the very first line of the file)
* Find section `group "${CI_PROJECT_NAME}"`
* Find section `env`. Inside it, set your new variable. For example:
```.hcl
job "${CI_PROJECT_NAME}-${NOMAD_DC}" {
...
    group "${CI_PROJECT_NAME}" {
    ...
        env {
            "ONEQ_TOKEN"   = "${ONEQ_TOKEN}"
            "SERVICE_URL"  = "${SERVICE_URL}"
            "GLORIOUS_VAR" = "${GLORIOUS_VAR}"
         }
    ...
    }
...
}
```

### How to deploy

1. Make some changes in the project

2. Push them to the repository. It will automatically run [lint](https://github.com/biocad/molprop-back/blob/master/iac/ci/lint.yaml) stage of the CI pipeline with some basic checks (e.g. code style)

3. If you want to build a new docker-image, create a new release. It will automatically run [build](https://github.com/biocad/molprop-back/blob/master/iac/ci/build.yaml) and [deploy-test](https://github.com/biocad/molprop-back/blob/master/iac/ci/deploy-test.yaml) stages of the pipeline

4. Go to [Gitlab pipelines](https://gitlab.math.bio/biocad/molprop-back/-/pipelines) of the project. There, you can manually run [deploy-dev](https://github.com/biocad/molprop-back/blob/master/iac/ci/deploy-dev.yaml) and  [deploy-prod](https://github.com/biocad/molprop-back/blob/master/iac/ci/deploy-prod.yaml) stages of the CI pipeline

5. That's it! Now you can find the project on a server. For example, for the "dev" version of the application, go to https://molprop-back.dev.net.biocad.ru/

If something went wrong without an obvious reason, contact DevOps team
```

## Project structure:
├── docker                          <- Files related to deployment via docker
|   ├── Dockerfile
|   ├── flask_cfg.py
|   └── gunicorn_cfg.py
|
├── iac                             <- Files related to deployment via docker
|   ├── ci
|   |   ├── build.yaml              <- Describes how and on what event to build project
|   |   ├── deploy-prod.yaml        <- Describes how and on what event to deploy project to pord
|   |   └── deploy-test.yaml        <- Describes how and on what event to deploy project to test
|   └── template-nomad.hcl          <- Get and set env vars. Set other variables of app 
├── tests          
|   └── test_api.py                 <- Api tests
|
├── molecad 
    ^--------------------------------- Name of service in snake case
|   ├── __init__.py                 <- Says that folder is module
|   ├── error_handlers.py           <- Caught errors handlers 
|   ├── errors.py                   <- Custom errors implementation
|   ├── log.py                      <- Loguru sink, formatter and logger patch for correlation_id
|   ├── main.py                     <- Main file with configuration
|   ├── py.typed                    <- Says that module is typed
|   ├── routes.py                   <- Application routes
|   └── validator.py                <- Implement your data validator here 
│
├── .env.sample                     <- Excample of .env file. You can use it to set environment variables.
├── .gitignore                      <- Git will ignore matching files and folders
├── .gitlab-ci.yml                  <- Combines instructions from iac folder
├── pyproject.toml                  <- Dependencies. execute in console `poetry install` to install them. `poetry add LIB_NAME` to add new dependency.
├── README.md                       <- file you are reading
└── setup.cfg                       <- linter and isort settings
```

Repository created from [bcd-flask-template](https://github.com/biocad/bcd-flask-template) using template version [1b3fa43431d12f520e406712812fbee23beb7b49](https://github.com/biocad/bcd-flask-template/commit/1b3fa43431d12f520e406712812fbee23beb7b49).