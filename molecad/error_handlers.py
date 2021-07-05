from loguru import logger
from starlette.responses import JSONResponse

from molecad.errors import BaseAppException


def app_error_handler(request, e: BaseAppException):
    return JSONResponse(content=e.to_dict(), status_code=e.error_code)


def exception(request, e: Exception):
    logger.exception("Unexpected error occurred. {}", e)
    return JSONResponse(
        content={
            "error": "unexpected",
            "params": {
                "message": str(e).replace("\n", " ") + str(getattr(e, "original_exception", ""))
            },
        },
        status_code=500,
    )
