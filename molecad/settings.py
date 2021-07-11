from pydantic import (
    BaseSettings,
    Field,
)


class Settings(BaseSettings):
    mongo_host = Field("127.0.0.1", env="MONGO_HOST")
    mongo_port = Field(27017, env="MONGO_PORT")
    mongo_user = Field("", env="MONGO_USER")
    mongo_password = Field("", env="MONGO_PASSWORD")
    mongo_db_name = Field("molecad", env="MONGO_DB_NAME")
    mongo_db_collection = Field("mol_props", env="MONGO_DB_COLLECTION")

    project_dir = Field(".", env="PROJ_DIR")
    fetch_dir = Field("./molecad/downloads/fetch", env="FETCH_DIR")
    split_file = Field("./molecad/downloads/split", env="SPLIT_DIR")
    saved_json = Field("", env="JSON")

    component_name = "molecad"
    api_version: str = "/v0"

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


setup = Settings()
