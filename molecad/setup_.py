import logging
import sys

from fastapi import FastAPI
from loguru import logger
from pydantic import (
    BaseSettings,
    Field,
)
from pymongo import MongoClient

from molecad.log import DevelopFormatter, JsonSink


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

    component_name = "molprop-back"

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"

    @property
    def webhook_uri(self):
        return self.service_url + self.api_version + "/webhooks"

    @property
    def solubility_api_url(self):
        return "https://chem-solubility-api.prod.net.biocad.ru/v0"

    @property
    def oneq(self) -> OneQAsync:
        if self.oneq_ is None:
            self.oneq_ = OneQAsync(OneQConfig(url=self.oneq_url, token=self.oneq_token))
        return self.oneq_

    @property
    def version(self):
        import importlib.metadata as importlib_metadata

        return importlib_metadata.version("molprop_back")

    def init_app(self):
        from starlette_context import plugins
        from starlette_context.middleware import ContextMiddleware

        app = FastAPI(title=self.component_name, version=self.version)
        plugins = [plugins.RequestIdPlugin(), plugins.CorrelationIdPlugin()]
        app.add_middleware(ContextMiddleware, plugins=plugins)
        return app

    def setup_logging(self):
        class InterceptHandler(logging.Handler):
            def emit(self, record):
                # Get corresponding Loguru level if it exists
                try:
                    level = logger.level(record.levelname).name
                except ValueError:
                    level = record.levelno

                # Find caller from where originated the logged message
                frame, depth = logging.currentframe(), 2
                while frame.f_code.co_filename == logging.__file__:
                    frame = frame.f_back
                    depth += 1

                logger.opt(depth=depth, exception=record.exc_info).log(level, record.getMessage())

        logging.getLogger().handlers = [InterceptHandler()]
        logger.remove()
        if self.env == "development":
            develop_fmt = DevelopFormatter(self.component_name)
            logger.add(sys.stdout, format=develop_fmt, level="DEBUG", backtrace=True, diagnose=True)
        else:
            json_sink = JsonSink(self.component_name)
            logger.add(json_sink)

        logging.getLogger("uvicorn.access").handlers = [InterceptHandler()]

    @property
    def db_collection(self):
        if self.collection_ is None:
            if self.mongo_user:
                mongo_client = MongoClient(
                    f"mongodb://{self.mongo_user}:{self.mongo_password}@{self.mongo_host}"
                )
            else:  # For some reason this syntax doesn't work with a real database
                mongo_client = MongoClient(
                    host=self.mongo_host,
                    port=self.mongo_port,
                    username=self.mongo_user,
                    password=self.mongo_password,
                    authSource=self.mongo_db_name,
                )

            db = mongo_client[self.mongo_db_name]
            self.collection_ = db[self.mongo_db_collection]

        return self.collection_


setup = Settings()


def setup_routing(app: FastAPI):
    from .routes import route
    from .webhooks import webhooks_route

    app.include_router(route, prefix=setup.api_version + "/api")
    app.include_router(webhooks_route, prefix=setup.api_version + "/webhooks")