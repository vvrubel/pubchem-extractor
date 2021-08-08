import uvicorn
from fastapi import FastAPI

# TODO: change to motor
from pymongo import MongoClient

from .routes import router
from .settings import settings

app = FastAPI(title=settings.component_name, version=settings.version)

app.include_router(router, prefix=settings.api_version)


@app.on_event("startup")
def startup_db_client():
    app.mongodb_client = MongoClient(settings.mongo_uri)
    app.mongodb = app.mongodb_client[settings.db_name]


@app.on_event("shutdown")
async def shutdown_db_client():
    app.mongodb_client.close()


if __name__ == "__main__":
    uvicorn.run("molecad.main:app", reload=True, host=settings.api_host, port=settings.api_port)
