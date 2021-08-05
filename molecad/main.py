import uvicorn
from fastapi import FastAPI
from pymongo import MongoClient

from .error_handler import app_error_handler, exception
from .errors import BaseAppException
from .routes import router
from .settings import settings

app = FastAPI(title=settings.component_name, version=settings.version)


@app.on_event("startup")
def startup_db_client():
    app.mongodb_client = MongoClient(settings.mongo_uri)
    app.mongodb = app.mongodb_client[settings.db_name]


@app.on_event("shutdown")
def shutdown_db_client():
    app.mongodb_client.close()


settings.setup_logging()
app.include_router(router, tags=["compound"], prefix=f"/{settings.api_version}")
app.add_exception_handler(BaseAppException, app_error_handler)
app.add_exception_handler(Exception, exception)


if __name__ == "__main__":
    uvicorn.run("molecad.main:api", reload=True, host=settings.api_host, port=settings.api_port)
