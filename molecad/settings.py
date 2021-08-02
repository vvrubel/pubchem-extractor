from pathlib import Path

from pydantic import BaseSettings, Field, HttpUrl


class Settings(BaseSettings):
    env: str = Field("PROD", env="ENV")
    project_dir: Path = Field(".", env="PROJ_DIR")
    fetch_dir: Path = Field("./data/fetch", env="FETCH_DIR")
    split_dir: Path = Field("./data/split", env="SPLIT_DIR")

    api_url: HttpUrl = Field("http://127.0.0.1:8000", env="API_URL")
    api_version: str = Field("", env="API_VERSION")
    app_url: HttpUrl = Field("http://127.0.0.1:8050", env="APP_URL")

    mongo_host: str = Field("127.0.0.1", env="MONGO_HOST")
    mongo_port: int = Field(27017, env="MONGO_PORT")
    mongo_user: str = Field("", env="MONGO_USER")
    mongo_password: str = Field("", env="MONGO_PASSWORD")
    mongo_auth_source: str = Field("admin", env="MONGO_AUTH_SOURCE")

    db_name: str = Field("", env="MONGO_DB_NAME")
    properties: str = Field("properties", env="MONGO_PROPERTIES_COLLECTION")
    molecules: str = Field("molecules", env="MONGO_MOLECULES_COLLECTION")
    mfp_counts: str = Field("mfp_counts", env="MONGO_MFP_COUNTS_COLLECTION")

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


settings = Settings()
