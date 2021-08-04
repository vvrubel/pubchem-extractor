from loguru import logger


class BadDomainError(Exception):
    pass


class BadNamespaceError(Exception):
    pass


class BadOperationError(Exception):
    pass


class EmptySmilesError(Exception):
    pass


class BaseAppException(Exception):
    def __init__(self, message: str = "") -> None:
        super().__init__(message)
        self.message = message
        logger.trace("Exception '{}' was raised. {}.", self.__class__.__name__, self.message)

    @property
    def error_code(self) -> int:
        raise NotImplementedError("Error code should be implemented.")

    def to_dict(self) -> dict:
        raise NotImplementedError("`to_dict` is not implemented.")


class NoDatabaseRecordError(BaseAppException):
    def __init__(self, message: str = "", function_name: str = "") -> None:
        super().__init__(message=message)
        self.function_name = function_name

    @property
    def error_code(self) -> int:
        return 404

    def to_dict(self) -> dict:
        return {
            "error": "database.record.missing",
            "params": {"function": self.function_name},
            "message": str(self),
        }

    def __str__(self) -> str:
        return "Couldn't find record in the database"


class BadRequestError(BaseAppException):
    def to_dict(self) -> dict:
        return {"error": "Bad request"}

    @property
    def error_code(self) -> int:
        return 400


class ResultUnexpectedError(BadRequestError):
    def to_dict(self) -> dict:
        return {"error": "result.unexpected", "message": str(self)}

    def __str__(self):
        return "Unexpected error, when trying to get result"


class UnknownPipelineError(BadRequestError):
    def to_dict(self) -> dict:
        return {"error": "pipeline.unknown", "message": str(self)}

    def __str__(self) -> str:
        return "Pipeline is unknown or unimplemented, can't run the task"
