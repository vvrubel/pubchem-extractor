import uvicorn
from fastapi import FastAPI
from motor.motor_asyncio import AsyncIOMotorClient

from .error_handler import app_error_handler, exception
from .errors import BaseAppException
from .routes import route
from .settings import settings

app = FastAPI(title=settings.component_name, version=settings.version)


@app.on_event("startup")
async def startup_db_client():
    app.mongodb_client = AsyncIOMotorClient(settings.mongo_uri)
    app.mongodb = app.mongodb_client[settings.db_name]


@app.on_event("shutdown")
async def shutdown_db_client():
    app.mongodb_client.close()


settings.setup_logging()
app.include_router(route, tags=["compound"], prefix=f"/{settings.api_version}/compound")
app.add_exception_handler(BaseAppException, app_error_handler)
app.add_exception_handler(Exception, exception)


if __name__ == "__main__":
    uvicorn.run("molecad.main:api", reload=True, host=settings.api_host, port=settings.api_port)
