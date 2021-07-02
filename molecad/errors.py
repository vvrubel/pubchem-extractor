class BaseAppException(Exception):
    def __init__(self, message: str = "") -> None:
        super().__init__(message)
        self.message = message
        logger.trace("Exception '{}' was raised. {}.", self.__class__.__name__, self.message)

    @property
    def error_code(self) -> int:
        raise NotImplementedError(
            "Error code should be implemented. Choose one from "
        )

    def to_dict(self) -> dict:
        raise NotImplementedError("`to_dict` is not implemented.")


class BadRequestError(BaseAppException):
    def to_dict(self) -> dict:
        return {"error": "Bad request"}

    @property
    def error_code(self) -> int:
        return 400


class BadDomainError(BadRequestError):
    def to_dict(self) -> dict:
        return {"error": "Bad domain", "message": str(self)}

    def __str__(self):
        return 'В данной версии сервиса поиск возможен только по базе данных "Compound"'


class BadNamespaceError(BadRequestError):
    def to_dict(self) -> dict:
        return {"error": "Bad namespace", "message": str(self)}

    def __str__(self):
        return 'Пространство имен для поиска по базе данных "Compound" задано некорректно'


class OperationError(BadRequestError):
    def to_dict(self) -> dict:
        return {"error": "Bad namespace", "message": str(self)}

    def __str__(self):
        return 'Ошибка при составлении операции'


class DirNotFoundError(BadRequestError, FileNotFoundError, FileExistsError):
    def to_dict(self) -> dict:
        return {"error": "Bad directory path", "message": str(self)}

    def __str__(self):
        return "Указанной директории не существует"

#
# class FileMissingError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "file.missing", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "Input file is missing"
#
#
# class WrongExtensionError(BadRequestError):
#     def __init__(self, message: str = "", allowed_extensions: tuple = (".smi",)) -> None:
#         super().__init__(message)
#         self.allowed_extensions = allowed_extensions
#
#     def to_dict(self) -> dict:
#         return {
#             "error": "file.wrong_extension",
#             "params": {"allowed_extensions": self.allowed_extensions},
#             "message": str(self),
#         }
#
#     def __str__(self) -> str:
#         extensions_string: str = ", ".join(self.allowed_extensions)
#         return f"Wrong file extension. Allowed extensions: {extensions_string}"
#
#
# class EmptyFileError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "file.empty", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "File is empty"
#
#
# class ContainsEmptyStringsError(BadRequestError):
#     def __init__(self, message: str = "", line: int = 0) -> None:
#         super().__init__(message)
#         self.line = line
#
#     def to_dict(self) -> dict:
#         return {
#             "error": "file.has_empty_strings",
#             "params": {"line": self.line},
#             "message": str(self),
#         }
#
#     def __str__(self) -> str:
#         return f"File contains an empty string at the line {self.line}"
#
#
# class InvalidMoleculeError(BadRequestError):
#     def __init__(self, message: str = "", wrong_molecule: str = "") -> None:
#         super().__init__(message)
#         self.wrong_molecule = wrong_molecule
#
#     def to_dict(self) -> dict:
#         return {
#             "error": "file.molecules.invalid",
#             "params": {"wrong_molecule": self.wrong_molecule},
#             "message": str(self),
#         }
#
#     def __str__(self) -> str:
#         return f"Can't read molecule '{self.wrong_molecule}'. Please provide valid SMILES formula."
#
#
# class LigPrepError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "ligprep.error", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "There was an error with LigPrep run"
#
#
# class NoLigPrepOutputError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "ligprep.output.missing", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "LigPrep didn't return an output"
#
#
# class QikpropError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "qikprop.error", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "There was an error with QikProp run"
#
#
# class QikpropReturnedNoCSVError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "qikprop.csv.missing", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "QikProp didn't return a CSV output"
#
#
# class UnknownPipelineError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "pipeline.unknown", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "Pipeline is unknown or unimplemented, can't run the task"
#
#
# class NoDatabaseRecordError(BaseAppException):
#     def __init__(self, message: str = "", function_name: str = "") -> None:
#         super().__init__(message=message)
#         self.function_name = function_name
#
#     @property
#     def error_code(self) -> int:
#         return 404
#
#     def to_dict(self) -> dict:
#         return {
#             "error": "database.record.missing",
#             "params": {"function": self.function_name},
#             "message": str(self),
#         }
#
#     def __str__(self) -> str:
#         return "Couldn't find record in the database"
#
#
# class EpikError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "epik.error", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "There was an error with Epik run"
#
#
# class EpikReturnedNoOutputError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "epik.output.missing", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "Epik didn't return any output"
#
#
# class TaskIdNotFoundError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "task_id.not_found", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "Task ID not found"
#
#
# class UnknownResultTypeError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "result.type.unknown", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "Can't understand the result type of the task"
#
#
# class UnprocessedWebhookError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "webhook.catch.failed", "message": str(self)}
#
#     def __str__(self):
#         return 'Task result is "Completed", but we did not receive a webhook'
#
#
# class ResultUnexpectedError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "result.unexpected", "message": str(self)}
#
#     def __str__(self):
#         return "Unexpected error, when trying to get result"
#
#
# class IncorrectTaskIDFormatError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "uuid.format.incorrect", "message": str(self)}
#
#     def __str__(self) -> str:
#         return "Incorrect format of task ID"
#
#
# class AmbiguousTaskStatusError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "task.status.ambiguous", "message": str(self)}
#
#     def __str__(self):
#         return (
#             "Task's status is completed, but we didn't receive a webhook. Please reload the "
#             "page. If this error remains, contact us, and provide the details about the task"
#         )
#
#
# class PipelineMissingTaskIDError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "pipeline.task_id.missing", "message": str(self)}
#
#     def __str__(self):
#         return "A step of the pipeline has no task ID"
#
#
# class PipelineMissingError(BadRequestError):
#     def to_dict(self) -> dict:
#         return {"error": "pipeline.missing", "message": str(self)}
#
#     def __str__(self):
#         return "The task doesn't have pipeline"
