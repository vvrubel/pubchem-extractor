from pathlib import Path

from pydantic import BaseSettings, Field


class Settings(BaseSettings):
    env: str = Field("dev", env="ENV")
    project_dir: Path = Field(".", env="PROJ_DIR")
    fetch_dir: Path = Field("./data/fetch", env="FETCH_DIR")
    split_dir: Path = Field("./data/split", env="SPLIT_DIR")

    mongo_host: str = Field("127.0.0.1", env="MONGO_HOST")
    mongo_port: int = Field(27017, env="MONGO_PORT")
    mongo_user: str = Field("", env="MONGO_USER")
    mongo_password: str = Field("", env="MONGO_PASSWORD")
    mongo_auth_source: str = Field("admin", env="MONGO_AUTH_SOURCE")

    db_name: str = Field("compound", env="MONGO_DB_NAME")
    prop_collection: str = Field("properties", env="MONGO_PROPERTIES_COLLECTION")
    mol_collection: str = Field("molecules", env="MONGO_MOLECULES_COLLECTION")
    mfp_collection: str = Field("mfp_counts", env="MONGO_MFP_COUNTS_COLLECTION")

    api_host: str = Field("127.0.0.1", env="API_HOST")
    api_port: int = Field(8000, env="API_PORT")
    api_version: str = Field("", env="API_VERSION")

    component_name: str = "molecad"

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"

    @property
    def mongo_uri(self):
        return (
            f"mongodb://{self.mongo_user}:{self.mongo_password}@{self.mongo_host}:"
            f"{self.mongo_port}/{self.mongo_auth_source}"
        )


settings = Settings()
