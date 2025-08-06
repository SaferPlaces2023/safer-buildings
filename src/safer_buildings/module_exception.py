import os
import traceback

from . import _consts
from .module_log import Logger, is_debug_mode


# DOC: Exception mapping for status code in exeptions
_EXCEPTION__OUTPUT_CODE__MAP = {
    # DOC: 4xx Client Errors
    FileNotFoundError: 404,                 # Not Found

    # DOC: 5xx Server Errors
    ConnectionError: 502,                   # Bad Gateway
    TimeoutError: 504,                      # Gateway Timeout
}


class CustomException(Exception):
    """
    Custom exception class for safer_buildings module.
    """
    def __init__(self, message, status_code=None, ex_type=None, ex_parent=None, *args):
        super().__init__(message, *args)
        self.message = message
        self.status_code = status_code or self.__class__.status_code
        self.status_code = _EXCEPTION__OUTPUT_CODE__MAP.get(ex_type, self.status_code)
        self.ex_type = ex_type or self.__class__
        self.ex_parent = ex_parent

    @classmethod
    def from_exception(cls, exception: Exception):
        """
        Create a CustomException from an existing native exception.
        """
        return cls(
            message = str(exception),
            status_code = cls.status_code,
            ex_type = type(exception),
            ex_parent = exception
        )

    def to_dict(self):
        """
        Convert the exception to a dictionary representation.
        """
        additional_info = dict()

        if is_debug_mode():
            if self.ex_parent:
                additional_info['traceback'] = ''.join(traceback.format_exception(self.ex_parent))

        ex_dict = {
            "exception_source": type(self).__name__,
            "message": self.message,
            _consts._ERROR_CODE_KEY: self.status_code,
            "exception_type": self.ex_type.__name__,
            ** additional_info
        }

        return ex_dict


# DOC: Custom exception for module args
class ArgsException(CustomException):
    """
    Custom exception for argument validation errors.
    """
    status_code = 400


# DOC: Custom exception for module retriever
class RetrieverException(CustomException):
    """
    Custom exception for module retriever errors.
    """
    status_code = 500


# DOC: cosumo exception for module flood
class FloodException(CustomException):
    """
    Custom exception for module flood errors.
    """
    status_code = 500


# DOC: Custom exception for module stats
class StatsException(CustomException):
    """
    Custom exception for module stats errors.
    """
    status_code = 500


# DOC: Custom exception for module add_ops
class AddOpsException(CustomException):
    """
    Custom exception for module add_ops errors.
    """
    status_code = 500


# DOC: Custom exception for module outputs
class OutputsException(CustomException):
    """
    Custom exception for module outputs errors.
    """
    status_code = 500