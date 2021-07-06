from pydantic import (
    BaseSettings,
    Field,
)


class Settings(BaseSettings):
    env: str = Field("production", env="ENV")
    api_version: str = "/v0"
    service_url: str = Field("http://localhost:5000", env="SERVICE_URL")

    mongo_host = Field("127.0.0.1", env="MONGO_HOST")
    mongo_port = Field(27017, env="MONGO_PORT")
    mongo_user = Field("", env="MONGO_USER")
    mongo_password = Field("", env="MONGO_PASSWORD")
    mongo_db_name = Field("molecad", env="MONGO_DB_NAME")
    mongo_db_collection = Field("properties", env="MONGO_DB_COLLECTION")
    mongo_auth_source = Field("admin", env="MONGO_AUTH_SOURCE")

    data_file = Field("", env="FILE")
    data_dir = Field("", env="OUT_DIR")

    component_name = "molecad"

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


setup = Settings()
